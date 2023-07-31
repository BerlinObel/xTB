import sys
import os
import pandas as pd
import copy

from rdkit import Chem
from rdkit.Chem import AllChem

script_dir = os.path.dirname(os.path.abspath(__file__))
sys.path.append(script_dir)
from settings import QMC_PATH
sys.path.append(QMC_PATH)

from QMC.qmmol import QMMol
from QMC.qmconf import QMConf
from QMC.calculator.xtb import xTB

from QMC.conformers.create_conformers import RotatableBonds

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
        bond = reactant_norm.GetBondBetweenAtoms(atom_indices[0], atom_indices[1])
        broken_bond_reactant = Chem.FragmentOnBonds(reactant_norm, [bond.GetIdx()], addDummies=False)
        
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
    new_product_order = find_atom_mapping_reactant_to_product(reactant, product)
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

def find_ground_state_conformers(name, rdkit_conf, charge, multiplicity, num_cpus):
    """Find ground state conformers"""
    # Create QMMol object and add conformer
    qmmol = QMMol()
    qmmol.add_conformer(rdkit_conf, fmt='rdkit', label=name, charged_fragments=True, set_initial=True)

    # Find number of rotatable bonds, accounting for triple bonds
    molecule = qmmol.get_rdkit_mol()
    triple_bond_smarts = '[c]:[c]-[CH0]#[CH0]'
    num_extra_rotatable_bonds = int(len(molecule.GetSubstructMatches(Chem.MolFromSmarts(triple_bond_smarts))) / 2)

    # Compute number of conformers to find
    num_rotatable_bonds = len(RotatableBonds(qmmol.initial_conformer.get_rdkit_mol())) + num_extra_rotatable_bonds
    num_conformers = 5 + 5 * num_rotatable_bonds
    qmmol.create_random_conformers(threads=num_cpus, num_confs=num_conformers)

    # Set up and run optimization
    xtb_parameters = {'method': 'gfn2', 'opt': 'tight', 'cpus': 1}
    qmmol.calc = xTB(parameters=xtb_parameters)
    qmmol.optimize(num_procs=num_cpus, keep_files=False)

    # Get most stable conformer. If most stable conformer is not identical to initial conformer, try next lowest.
    initial_smiles = Chem.MolToSmiles(Chem.RemoveHs(qmmol.initial_conformer.get_rdkit_mol()))
    lowest_energy_conformer = qmmol.nlowest(1)[0]
    conf_smiles = Chem.MolToSmiles(Chem.RemoveHs(lowest_energy_conformer.get_rdkit_mol()))

    index = 1
    while initial_smiles != conf_smiles:
        lowest_energy_conformer = qmmol.nlowest(index+1)[-1]
        conf_smiles = Chem.MolToSmiles(Chem.RemoveHs(lowest_energy_conformer.get_rdkit_mol()))
        index += 1
        if len(qmmol.conformers) < index:
            sys.exit('Error: no conformers match the initial input')
    return lowest_energy_conformer

def perform_ground_state_search(name, smi, charge, multiplicity, num_cpus):
    """Perform ground state search given a SMILES string"""
    reactant_smi, product_smi = generate_product_from_reactant(smi)

    reactant_mol = Chem.AddHs(Chem.MolFromSmiles(reactant_smi))
    # Attempt to reorder product to match reactant
    try:
        print("Attempting to reorder product...")
        product_mol = reorder_product_to_match_reactant(reactant_mol, Chem.AddHs(Chem.MolFromSmiles(product_smi)))
    except:
        print("Product reorder failed, using original order...")
        product_mol = Chem.AddHs(Chem.MolFromSmiles(product_smi))

    for molecule, suffix in [(reactant_mol, "_r"), (product_mol, "_p")]:
        AllChem.EmbedMolecule(molecule)
        rdkit_conf = molecule.GetConformer()
        molecule_name = name + suffix
        # Conduct ground state conformer search
        if suffix == '_r':
            try:
                print("Starting reactant conformer search...")
                reactant_qmconf = find_ground_state_conformers(molecule_name, rdkit_conf, charge, multiplicity, num_cpus)
            except:
                print("Reactant conformer search failed, retrying without reordering...")
                retry_count = 0
                reactant_qmconf = None
                while reactant_qmconf is None and retry_count < 40:
                    try:
                        reactant_qmconf = retry_ground_state_search(charge, multiplicity, num_cpus, name, reactant_smi, suffix)
                    except:
                        retry_count += 1
                        print(f"Retry failed: {retry_count}")
        elif suffix == '_p':
            try:
                print("Starting product conformer search...")
                product_qmconf = find_ground_state_conformers(molecule_name, rdkit_conf, charge, multiplicity, num_cpus)
            except:
                print("Product conformer search failed, retrying without reordering...")
                retry_count = 0
                product_qmconf = None
                while product_qmconf is None and retry_count < 40:
                    try:
                        product_qmconf = retry_ground_state_search(charge, multiplicity, num_cpus, name, product_smi, suffix)
                    except:
                        retry_count += 1
                        print(f"Retry failed: {retry_count}")

    # Calculate and store energy difference between product and reactant
    energy_diff = product_qmconf.results['energy'] - reactant_qmconf.results['energy']
    return reactant_qmconf, product_qmconf, energy_diff


def retry_ground_state_search(charge, multiplicity, num_cpus, name, smi, suffix):
    """Retry ground state search if initial attempt fails"""
    molecule = Chem.AddHs(Chem.MolFromSmiles(smi))
    AllChem.EmbedMolecule(molecule)
    rdkit_conf = molecule.GetConformer()
    molecule_name = name + suffix
    return find_ground_state_conformers(molecule_name, rdkit_conf, charge, multiplicity, num_cpus)


if __name__ == '__main__':
    num_cpus = 2
    dataset = pd.read_csv(sys.argv[1])
    # Find energy difference and product for each molecule in the dataset
    results = []
    for index, compound in dataset.iterrows():
        # Get reactant conformer, product conformer, and energy difference
        reactant_qmconf, product_qmconf, energy_diff = perform_ground_state_search(
            str(compound.comp_name),
            compound.smiles,
            compound.charge,
            compound.multiplicity,
            num_cpus
        )
        results.append({
            'rep': str(compound.comp_name),
            'reac': reactant_qmconf,
            'prod': product_qmconf,
            'storage': energy_diff
        })

    # Gather results
    results_df = pd.DataFrame(results)
    pickle_filename = sys.argv[1].split('/')[0].split('.')[0] + '.pkl'
    print(f"pickle filename: {pickle_filename}")
    print(results_df)
    results_df.to_pickle(pickle_filename)

