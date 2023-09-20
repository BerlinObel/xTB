import sys
import os
import pandas as pd
import copy
import re

import concurrent.futures
from concurrent.futures import ThreadPoolExecutor, TimeoutError

import time
from rdkit import Chem
from rdkit.Chem import AllChem

# Add the directory of the script to the Python path
sys.path.append(os.path.dirname(os.path.realpath(__file__)))

from db_utils import update_data, retrieve_data
from settings import QMC_PATH
sys.path.append(QMC_PATH)

from QMC.conformers.create_conformers import RotatableBonds
from QMC.calculator.xtb import xTB
from QMC.qmconf import QMConf
from QMC.qmmol import QMMol
from utils import execute_shell_command, get_total_energy_xtb, get_statistics

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
    product_smi = reactant_smi.replace("/N=N/", "/N=N\\")

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
    output = execute_shell_command(command_input, shell=False, timeout=100)
    
    energy = get_total_energy_xtb(output)
    
    optimized_conformer = QMConf(xyz_result, fmt='xyz', label=conformer_name)
    optimized_conformer.results['energy'] = energy
     
    return optimized_conformer

def find_and_validate_lowest_energy_conformer(initial_smiles, optimized_confs):

    lowest_energy_conformer = min(optimized_confs, key=lambda x: x.results['energy'])
    conf_smiles = Chem.MolToSmiles(Chem.RemoveHs(lowest_energy_conformer.get_rdkit_mol()))
    
    index = 1
    while initial_smiles != conf_smiles:
        optimized_confs.remove(lowest_energy_conformer)
        new_lowest_energy_conformer = min(optimized_confs, key=lambda x: x.results['energy'])
        conf_smiles = Chem.MolToSmiles(Chem.RemoveHs(new_lowest_energy_conformer.get_rdkit_mol()))
        
        index += 1
        if len(optimized_confs) < index:
            sys.exit('Error: no conformers match the initial input')
            
    return lowest_energy_conformer

def find_ground_state_conformers(name, rdkit_conf, charge, multiplicity, num_cpus):
    all_conformers = create_conformers(name, rdkit_conf, charge, multiplicity, num_cpus)
    optimized_confs = []  # List to gather all optimized conformers

    with ThreadPoolExecutor(max_workers=num_cpus) as executor:
        future_to_conformer = {executor.submit(optimize_conformer, conformer): conformer for conformer in all_conformers}
        
        for future in concurrent.futures.as_completed(future_to_conformer):
            conformer = future_to_conformer[future]

            try:
                optimized_conformer = future.result()
                # Append the successfully computed optimized conformer to the list
                optimized_confs.append(optimized_conformer)
            except TimeoutError:
                print(f"The task for {conformer.label} did not complete within the timeout")
            except Exception as exc:
                print(f"The task for {conformer.label} raised an exception: {exc}")

    initial_smiles = Chem.MolToSmiles(Chem.RemoveHs(all_conformers[0].get_rdkit_mol()))
    
    lowest_energy_conformer = find_and_validate_lowest_energy_conformer(initial_smiles, optimized_confs)

    energies = [conf.results['energy'] for conf in optimized_confs]
    energy_stats = get_statistics(energies)
    print("Energy stats")
    print(energy_stats)

    
    ref_mol = lowest_energy_conformer.get_rdkit_mol()
    rmsd_values = get_rmsd_values(ref_mol, optimized_confs)
    rmsd_stats = get_statistics(rmsd_values)
    print("RMSD stats")
    print(rmsd_stats)
    
    print(f"Lowest energy conf: {lowest_energy_conformer.label}")
    print(f"Energy: {lowest_energy_conformer.results['energy']}")
    
    return lowest_energy_conformer, energy_stats, rmsd_stats


def perform_ground_state_search(name, smi, charge, multiplicity, num_cpus):
    """Perform ground state search given a SMILES string"""
    print(f"Create product for {name}:")
    reactant_smi, product_smi = generate_product_from_reactant(smi)

    reactant_mol = Chem.AddHs(Chem.MolFromSmiles(reactant_smi))
    # Attempt to reorder product to match reactant
    try:
        print("Attempting to reorder product...")
        product_mol = reorder_product_to_match_reactant(
            reactant_mol, Chem.AddHs(Chem.MolFromSmiles(product_smi)))
    except:
        print("Product reorder failed, using original order...")
        product_mol = Chem.AddHs(Chem.MolFromSmiles(product_smi))

    for molecule, suffix in [(reactant_mol, "_r"), (product_mol, "_p")]:
        AllChem.EmbedMolecule(molecule)
        rdkit_conf = molecule.GetConformer()
        molecule_name = name + suffix
        stat_dictionary = {}
        # Conduct ground state conformer search
        if suffix == '_r':
            try:
                print(f"Starting reactant conformer search for {molecule_name}")
                reactant_qmconf, energy_stats, rmsd_stats = find_ground_state_conformers(
                    molecule_name, rdkit_conf, charge, multiplicity, num_cpus)
                stat_dictionary["Reactant_Energy"] = energy_stats
                stat_dictionary["Reactant_RMSD"] = rmsd_stats
                
            except:
                print("Reactant conformer search failed, retrying without reordering...")
                retry_count = 0
                reactant_qmconf = None
                while reactant_qmconf is None and retry_count < 40:
                    try:
                        reactant_qmconf, energy_stats, rmsd_stats = retry_ground_state_search(
                        molecule_name, rdkit_conf, charge, multiplicity, num_cpus)
                        stat_dictionary["Reactant_Energy"] = energy_stats
                        stat_dictionary["Reactant_RMSD"] = rmsd_stats
                    except:
                        retry_count += 1
                        print(f"Retry failed: {retry_count}")
        elif suffix == '_p':
            try:
                print(f"Starting product conformer search for {molecule_name}")
                product_qmconf, energy_stats, rmsd_stats = find_ground_state_conformers(
                    molecule_name, rdkit_conf, charge, multiplicity, num_cpus)
                stat_dictionary["Product_Energy"] = energy_stats
                stat_dictionary["Product_RMSD"] = rmsd_stats
            except:
                print("Product conformer search failed, retrying without reordering...")
                retry_count = 0
                product_qmconf = None
                while product_qmconf is None and retry_count < 40:
                    try:
                        product_qmconf, energy_stats, rmsd_stats = retry_ground_state_search(
                        molecule_name, rdkit_conf, charge, multiplicity, num_cpus)
                        stat_dictionary["Product_Energy"] = energy_stats
                        stat_dictionary["Product_RMSD"] = rmsd_stats
                    except:
                        retry_count += 1
                        print(f"Retry failed: {retry_count}")

    # Calculate and store energy difference between product and reactant
    energy_diff = product_qmconf.results['energy'] - \
        reactant_qmconf.results['energy']
    return reactant_qmconf, product_qmconf, energy_diff, stat_dictionary


def retry_ground_state_search(charge, multiplicity, num_cpus, name, smi, suffix):
    """Retry ground state search if initial attempt fails"""
    molecule = Chem.AddHs(Chem.MolFromSmiles(smi))
    AllChem.EmbedMolecule(molecule)
    rdkit_conf = molecule.GetConformer()
    molecule_name = name + suffix
    return find_ground_state_conformers(molecule_name, rdkit_conf, charge, multiplicity, num_cpus)

def main(batch_id, num_cpus):
    # Load the data
    filters = {'BatchID': batch_id, 'CalculationStage': 'storage'}
    data = retrieve_data(filters)
    print(data)
    # Find energy difference and product for each molecule in the dataset
    results = []
    for index, compound in data.iterrows():
        print(f"Index: {index+1}/{len(data)}, Compound: {compound.HashedName}")
        
        # Get reactant conformer, product conformer, and energy difference
        reactant_qmconf, product_qmconf, energy_diff, stat_dictionary = perform_ground_state_search(
            str(compound.HashedName),
            compound.Smiles,
            compound.Charge,
            compound.Multiplicity,
            num_cpus
        )
        results.append({
            'HashedName': str(compound.HashedName),
            'ReactantObject': reactant_qmconf,
            'ProductObject': product_qmconf,
            'StorageEnergy': energy_diff,
            'CalculationStage': 'storage_completed',
            'GroundStateStats': stat_dictionary
        })
    
    # Gather results
    results_df = pd.DataFrame(results)
    print("Finished Batch:")
    print(results_df)
    
    # Update the database with the results
    update_data(results_df)

if __name__ == '__main__':
    batch_id = sys.argv[1]
    num_cpus = 8
    main(batch_id, num_cpus)
