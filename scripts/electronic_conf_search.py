import sys
import os
import pandas as pd
import copy
import re

import concurrent.futures
from concurrent.futures import ThreadPoolExecutor, TimeoutError, as_completed

import time
from rdkit import Chem
from rdkit.Chem import AllChem, rdchem
# Add the directory of the script to the Python path
sys.path.append(os.path.dirname(os.path.realpath(__file__)))

from db_utils import update_data, retrieve_data
from settings import QMC_PATH
sys.path.append(QMC_PATH)

from QMC.conformers.create_conformers import RotatableBonds
from QMC.calculator.xtb import xTB
from QMC.qmconf import QMConf
from QMC.qmmol import QMMol
from utils import execute_shell_command, get_total_energy_xtb, get_statistics, format_time

def find_atom_mapping_reactant_to_product(reactant, product):
    """Find and return new atom order in product based on reactant"""
    # Copy molecules to prevent mutation of original molecules
    reactant_copy = copy.deepcopy(reactant)
    product_copy = copy.deepcopy(product)

    # Kekulize molecules to make aromaticity explicit
    Chem.Kekulize(reactant_copy, clearAromaticFlags=True)
    Chem.Kekulize(product_copy, clearAromaticFlags=True)

    # Normalize molecules to only compare connectivity
    reactant_norm = normalize_molecule(reactant_copy)
    product_norm = normalize_molecule(product_copy)

    # Attempt to match structure of product to reactant by breaking a bond
    smarts_bond = Chem.MolFromSmarts('[CX4;H0;R]-[CX4;H1;R]')
    atom_indices = list(reactant_norm.GetSubstructMatch(smarts_bond))

    if len(atom_indices) != 0:
        bond = reactant_norm.GetBondBetweenAtoms(
            atom_indices[0], atom_indices[1])
        broken_bond_reactant = Chem.FragmentOnBonds(
            reactant_norm, [bond.GetIdx()], addDummies=False)

        # find new atom order for product
        product_order = product_norm.GetSubstructMatch(broken_bond_reactant)
    else:
        product_order = product_norm.GetSubstructMatch(reactant_norm)

    return product_order


def normalize_molecule(mol):
    """Normalize molecule to compare only connectivity"""
    # Change all bond types to single
    for bond in mol.GetBonds():
        if bond.GetBondType() != Chem.BondType.SINGLE:
            bond.SetBondType(Chem.BondType.SINGLE)
    # Remove formal charges
    for atom in mol.GetAtoms():
        if atom.GetFormalCharge() != 0:
            atom.SetFormalCharge(0)
    return mol


def reorder_product_to_match_reactant(reactant, product):
    """Change atom order of product to match reactant"""
    new_product_order = find_atom_mapping_reactant_to_product(
        reactant, product)
    reordered_product = Chem.RenumberAtoms(product, new_product_order)
    return reordered_product


def generate_product_from_reactant(reactant_smi):
    """Generate product from reactant SMILES string"""
    # Azobenzene
    if "/N=N/" in reactant_smi:
        # Replace /N=N/ with /N=N\ in the reactant to create the product
        product_smi = reactant_smi.replace("/N=N/", "/N=N\\")
    elif "N=N" in reactant_smi:
        # Replace N=N with /N=N/ in the reactant first
        reactant_smi_new = reactant_smi.replace("N=N", "/N=N/")
        # Then, replace /N=N/ with /N=N\ to create the product
        product_smi = reactant_smi_new.replace("/N=N/", "/N=N\\")
        reactant_smi = reactant_smi_new
    else:
        # If neither /N=N/ nor N=N are present, return None or some default value
        product_smi = None

    # NBD
    # nbd_reaction_smarts = "[C:1]1[C:2]2[C:3]=[C:4][C:5]1[C:6]=[C:7]2>>[C:1]1[C:2]2[C:3]3[C:4]4[C:5]1[C:6]4[C:7]23"

    # Reaction Smarts
    # __rxn__ = AllChem.ReactionFromSmarts(nbd_reaction_smarts)
    # reactant_mol =  Chem.MolFromSmiles(reactant_smi)
    #
    # product_mol = __rxn__.RunReactants((reactant_mol,))[0][0]
    # product_smi = Chem.MolToSmiles(product_mol)

    return reactant_smi, product_smi

def get_rmsd_values(ref_mol, optimized_confs):
    rmsd_values = []
    for conformer in optimized_confs:
        target_mol = conformer.get_rdkit_mol()
        rmsd = AllChem.GetBestRMS(ref_mol, target_mol)
        rmsd_values.append(rmsd)
    return rmsd_values
    
def create_conformers(name, rdkit_conf, charge, multiplicity, num_cpus):
    qmmol = QMMol()
    qmmol.add_conformer(rdkit_conf, fmt='rdkit', label=name, charged_fragments=True, set_initial=True)
    
    molecule = qmmol.get_rdkit_mol()
    triple_bond_smarts = '[c]:[c]-[CH0]#[CH0]'
    num_extra_rotatable_bonds = int(len(molecule.GetSubstructMatches(Chem.MolFromSmarts(triple_bond_smarts))) / 2)
    
    num_rotatable_bonds = len(RotatableBonds(qmmol.initial_conformer.get_rdkit_mol())) + num_extra_rotatable_bonds
    num_conformers = 5 + 5 * num_rotatable_bonds
    print(f"Number of conformers: {num_conformers}")
    
    qmmol.create_random_conformers(threads=num_cpus, num_confs=num_conformers)
    
    return qmmol.conformers

def optimize_conformer(conformer):

    conformer_name = conformer.label
    xyz_file = f"{conformer_name}.xyz"
    xyz_result = "xtbopt.xyz"
    conformer.write_xyz(to_file=True)
    command_input = f"xtb {xyz_file} --gfn 1 --opt crude"
    output = execute_shell_command(command_input, shell=False, timeout=300)
    
    if output is None:
        raise TimeoutError(f"Optimization for {conformer_name} timed out") 
    if output == 999:
        raise Exception("Calculation Error")
    else:       
        energy = get_total_energy_xtb(output)
        optimized_conformer = QMConf(xyz_result, fmt='xyz', label=conformer_name)
        optimized_conformer.results['energy'] = energy
        return optimized_conformer


def remove_chirality(mol):
    for atom in mol.GetAtoms():
        atom.SetChiralTag(rdchem.ChiralType.CHI_UNSPECIFIED)
    return mol

def find_and_validate_lowest_energy_conformer(initial_smiles, optimized_confs):

    try:
        # Filter out conformers where energy is None
        valid_conformers = [conf for conf in optimized_confs if conf.results['energy'] is not None]
    except Exception as exc:
        print("Error! No optimized conformers available: {exc}")
        return None
    
    lowest_energy_conformer = min(valid_conformers, key=lambda x: x.results['energy'])
    conf_smiles = Chem.MolToSmiles(Chem.RemoveHs(lowest_energy_conformer.get_rdkit_mol()), canonical=True)
    initial_smiles = Chem.MolToSmiles(Chem.RemoveHs(Chem.MolFromSmiles(initial_smiles)), canonical=True)
   
    while initial_smiles != conf_smiles:
        # Remove chirality before comparison
        initial_mol_no_chirality = remove_chirality(Chem.MolFromSmiles(initial_smiles))
        conf_mol_no_chirality = remove_chirality(lowest_energy_conformer.get_rdkit_mol())

        # Generate SMILES without chirality
        initial_smiles_no_chirality = Chem.MolToSmiles(Chem.RemoveHs(initial_mol_no_chirality), canonical=True)
        conf_smiles_no_chirality = Chem.MolToSmiles(Chem.RemoveHs(conf_mol_no_chirality), canonical=True)

        if initial_smiles_no_chirality == conf_smiles_no_chirality:
            break
        valid_conformers.remove(lowest_energy_conformer)
        
        lowest_energy_conformer = min(valid_conformers, key=lambda x: x.results['energy'])
        conf_smiles = Chem.MolToSmiles(Chem.RemoveHs(lowest_energy_conformer.get_rdkit_mol()), canonical=True)

        if len(valid_conformers) == 0:
            print("Error: no conformers match the initial input")
            return None
    
    return lowest_energy_conformer

def find_ground_state_conformers(name, rdkit_conf, charge, multiplicity, num_cpus, smi):
    all_conformers = create_conformers(name, rdkit_conf, charge, multiplicity, num_cpus)
    optimized_confs = []  # List to gather all optimized conformers
    time_outs = 0
    # Dictionary to keep track of the number of optimization attempts for each conformer
    optimization_attempts = {conformer: 0 for conformer in all_conformers}

    # Initial submission
    with ThreadPoolExecutor(max_workers=num_cpus) as executor:
        future_to_conformer = {executor.submit(optimize_conformer, conformer): conformer for conformer in all_conformers}
    
        while future_to_conformer:
            for future in as_completed(future_to_conformer):
                conformer = future_to_conformer.pop(future)

                try:
                    optimized_conformer = future.result()
                    optimized_confs.append(optimized_conformer)
                except TimeoutError:
                    time_outs += 1

                    if conformer not in optimization_attempts:
                        print(f"Conformer {conformer.label} not found in optimization_attempts!")
                        optimization_attempts[conformer] = 0  # Initialize it

                    optimization_attempts[conformer] += 1

                    # Check if we should retry
                    if optimization_attempts[conformer] < 2 or len(optimized_confs) < 2:
                        future_to_conformer[executor.submit(optimize_conformer, conformer)] = conformer
                except Exception as exc:
                    print(f"Optimizing {conformer.label} raised an exception: {exc}")
                    
    print(f"{time_outs} conformers timed out")
    initial_smiles = smi
    
    lowest_energy_conformer = find_and_validate_lowest_energy_conformer(initial_smiles, optimized_confs)

    if lowest_energy_conformer != None:

        energies = [conf.results['energy'] for conf in optimized_confs]
        energy_stats = get_statistics(energies)

        ref_mol = lowest_energy_conformer.get_rdkit_mol()
        rmsd_values = get_rmsd_values(ref_mol, optimized_confs)
        rmsd_stats = get_statistics(rmsd_values)
        print(f"Energy: {lowest_energy_conformer.results['energy']}")
        return lowest_energy_conformer, energy_stats, rmsd_stats
    else:
        return None


def perform_ground_state_search(name, smi, charge, multiplicity, num_cpus):
    """Perform ground state search given a SMILES string"""
    reactant_smi, product_smi = generate_product_from_reactant(smi)

    reactant_mol = Chem.AddHs(Chem.MolFromSmiles(reactant_smi))

    # Attempt to reorder product to match reactant
    try:
        product_mol = reorder_product_to_match_reactant(
            reactant_mol, Chem.AddHs(Chem.MolFromSmiles(product_smi)))
    except:
        print("Product reorder failed, using original order")
        product_mol = Chem.AddHs(Chem.MolFromSmiles(product_smi))

    reactant_qmconf = None
    product_qmconf = None
    stat_dictionary = {}

    for molecule, suffix in [(reactant_mol, "_r"), (product_mol, "_p")]:
        # Skip product optimization if reactant ground state is None
        if suffix == '_p' and reactant_qmconf is None:
            print("Skipping product optimization as reactant ground state is None.")
            continue
        
        AllChem.EmbedMolecule(molecule)
        rdkit_conf = molecule.GetConformer()
        molecule_name = name + suffix
        initial_smiles = Chem.MolToSmiles(Chem.RemoveHs(molecule), canonical=True)
        
        if suffix == '_r':
            molecule_type = "Reactant"
        else:
            molecule_type = "Product"
        
        try:
            print(f"{molecule_type} ground state search")
            qmconf, energy_stats, rmsd_stats = find_ground_state_conformers(
                molecule_name, rdkit_conf, charge, multiplicity, num_cpus, initial_smiles)
            stat_dictionary[f"{molecule_type}_Energy"] = energy_stats
            stat_dictionary[f"{molecule_type}_RMSD"] = rmsd_stats
            
            if molecule_type == "Reactant":
                reactant_qmconf = qmconf
            else:
                product_qmconf = qmconf
                
                
        except Exception as e:
            print(f"Unable to find ground state conformer raised exception: {e}")
            retry_count = 0
            qmconf = None
            while qmconf is None and retry_count < 3:
                try:
                    print(f"Retrying {molecule_type} ground state search")
                    initial_smiles = Chem.MolToSmiles(Chem.RemoveHs(molecule), canonical=True)    
                    qmconf, energy_stats, rmsd_stats = retry_ground_state_search(
                        molecule_name, charge, multiplicity, num_cpus, initial_smiles)
                    stat_dictionary[f"{molecule_type}_Energy"] = energy_stats
                    stat_dictionary[f"{molecule_type}_RMSD"] = rmsd_stats
                    
                    if suffix == '_r':
                        reactant_qmconf = qmconf
                    else:
                        product_qmconf = qmconf
                except Exception as e:
                    retry_count += 1
                    print(f"Retry failed: {retry_count}, Exception: {e}")

    # Calculate and store energy difference between product and reactant
    if reactant_qmconf is not None and product_qmconf is not None:
        energy_diff = product_qmconf.results['energy'] - reactant_qmconf.results['energy']
    else:
        energy_diff = None  # Or some other default value

    return reactant_qmconf, product_qmconf, energy_diff, stat_dictionary


def retry_ground_state_search(molecule_name, charge, multiplicity, num_cpus, smi):
    """Retry ground state search if initial attempt fails"""
    try:
        molecule_mol = Chem.MolFromSmiles(smi)
        new_mol = Chem.AddHs(molecule_mol)
        
        # Sanitize the molecule after adding hydrogens
        Chem.SanitizeMol(new_mol)
        
        # Generate 3D coordinates
        AllChem.EmbedMolecule(new_mol)
        
        # Get the conformer
        rdkit_conf = new_mol.GetConformer()
        
        # Find ground state conformers
        return find_ground_state_conformers(molecule_name, rdkit_conf, charge, multiplicity, num_cpus, smi)
    except:
        print(f"An error occurred while processing {molecule_name} for retry.")
        return None

def main(batch_id, num_cpus):
    # Load the data
    filters = {'BatchID': batch_id, 'CalculationStage': 'storage'}
    data = retrieve_data(filters)
    print(data)
    start_time = time.time()
    # Find energy difference and product for each molecule in the dataset
    results = []
    for index, compound in data.iterrows():
      print(f"""
Index: {index+1}/{len(data)}, Compound: {compound.HashedName}
""")
      
      # Get reactant conformer, product conformer, and energy difference
      reactant_qmconf, product_qmconf, energy_diff, stat_dictionary = perform_ground_state_search(
          str(compound.HashedName),
          compound.Smiles,
          compound.Charge,
          compound.Multiplicity,
          num_cpus
      )
      
      if reactant_qmconf is None or product_qmconf is None:
          results.append({
              'HashedName': str(compound.HashedName),
              'CalculationStage': 'storage'
          })
      else:
          results.append({
              'HashedName': str(compound.HashedName),
              'ReactantObject': reactant_qmconf,
              'ProductObject': product_qmconf,
              'StorageEnergy': energy_diff,
              'CalculationStage': 'storage_completed',
              'GroundStateStats': stat_dictionary
          })
          
    end_time = time.time()
    elapsed_time = end_time - start_time
    formatted_elapsed_time = format_time(elapsed_time)
    # Gather results
    results_df = pd.DataFrame(results)
    print(f"""
    /)＿/)☆
 ／(๑^᎑^๑)っ ＼
|￣∪￣  ￣ |＼／
|＿＿_＿＿|／
{batch_id} finished:
{results_df}

Time: {formatted_elapsed_time}
          """)

    update_data(results_df)

if __name__ == '__main__':
    batch_id = sys.argv[1]
    num_cpus = 8
    main(batch_id, num_cpus)
