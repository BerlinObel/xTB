import os
import sys
import pandas as pd
sys.path.append("/groups/kemi/brq616/speciale/opt/xTB/tQMC/QMC")

from qmmol import QMMol
from qmconf import QMConf

from rdkit.Chem import AllChem
from rdkit.Chem import Draw


from rdkit.Chem import rdchem
from tqdm import tqdm
import time
from rdkit import Chem, rdBase
import itertools
import hashlib
from rdkit.Chem import rdMolDescriptors as rdMD
from rdkit.Chem import Descriptors
import logging



# Set up logging
logging.basicConfig(level=logging.INFO)

def generate_molecules(base_mol, substituents, indices, max_subs):
    new_mols = set()
    failed = 0

    for num_subs in range(1, min(max_subs, len(indices), len(substituents)) + 1):
        for subs in itertools.combinations_with_replacement(substituents, num_subs):
            for inds in itertools.combinations(indices, num_subs):
                new_mol = Chem.RWMol(base_mol)
                for (sub, sub_index, bond_type), base_index in zip(subs, inds):
                    num_atoms_new = new_mol.GetNumAtoms()
                    for atom in sub.GetAtoms():
                        new_mol.AddAtom(atom)
                    for bond in sub.GetBonds():
                        new_mol.AddBond(bond.GetBeginAtomIdx() + num_atoms_new, bond.GetEndAtomIdx() + num_atoms_new, bond.GetBondType())
                    new_mol.AddBond(base_index, num_atoms_new + sub_index, bond_type)
                try:
                    Chem.SanitizeMol(new_mol)
                    new_mols.add(new_mol)
                except:
                    failed += 1

    logging.debug(f'Failed sanitizations: {failed}')                
    return new_mols

def generate_and_print_molecules(base_molecule, substituents, positions, max_subs, ring_side, position_type):
    all_new_molecules = set()
    logging.debug(f"Generating {ring_side} ring {position_type} positions with up to {max_subs} substitutions")
    for num_subs in range(0, max_subs+1):
        new_molecules = generate_molecules(base_molecule, substituents, positions, num_subs)
        all_new_molecules.update(new_molecules)
    return all_new_molecules

def add_atom_numbers(molecule):
    # Add atom index property to each atom
    for atom in molecule.GetAtoms():
        atom.SetProp('molAtomMapNumber', str(atom.GetIdx()))
        
def validate_molecules(molecule_tuples):
    invalid_molecules = []
    for (mol, index, bond) in molecule_tuples:
        if mol is None:
            invalid_molecules.append((mol, index, bond))
    return invalid_molecules

def get_hashed_label(smiles, length=10):
    if smiles is not None:
        hashed_label = hashlib.sha1(smiles.encode()).hexdigest()[:length]
        return hashed_label
    else:
        return None
# Generate SMILES for substituents
# each element in the list is a tuple: (substituent, connecting atom index, bond type)

electron_withdrawing_substituents = [
    (Chem.MolFromSmiles('FC(F)F'), 1, rdchem.BondType.SINGLE),  # Trifluoromethyl
    (Chem.MolFromSmiles('Cl'), 0, rdchem.BondType.SINGLE),  # Chloro
    # (Chem.MolFromSmiles('F'), 0, rdchem.BondType.SINGLE),  # Fluoro
    # (Chem.MolFromSmiles('I'), 0, rdchem.BondType.SINGLE),  # Iodo
    # (Chem.MolFromSmiles('Br'), 0, rdchem.BondType.SINGLE),  # Bromo
    # (Chem.MolFromSmiles('C=O'), 0, rdchem.BondType.SINGLE),  # Carbonyl
    # (Chem.MolFromSmiles('OC(=O)'), 1, rdchem.BondType.SINGLE),  # Carboxyl
    # (Chem.MolFromSmiles('NC(=O)'), 1, rdchem.BondType.SINGLE),  # Amide
    (Chem.MolFromSmiles('C#N'), 0, rdchem.BondType.TRIPLE),  # Cyano
    (Chem.MolFromSmiles('S(=O)(=O)O'), 0, rdchem.BondType.SINGLE),  # Sulfonic acid
    (Chem.MolFromSmiles('[N+](=O)[O-]'), 0, rdchem.BondType.SINGLE),  # Nitro
]

electron_donating_substituents = [
    (Chem.MolFromSmiles('N'), 0, rdchem.BondType.SINGLE),  # Amino
    # (Chem.MolFromSmiles('CO'), 1, rdchem.BondType.SINGLE),  # Methoxy
    # (Chem.MolFromSmiles('O'), 0, rdchem.BondType.SINGLE),  # Hydroxy
    # (Chem.MolFromSmiles('C'), 0, rdchem.BondType.SINGLE),  # Methyl
    (Chem.MolFromSmiles('CN'), 0, rdchem.BondType.SINGLE),  # Methylamine
]

electron_donating_rings = [
    (Chem.MolFromSmiles('C1=CC=C2C(=C1)C=CS2'), 0, rdchem.BondType.SINGLE),  # Benzothiophene
    # (Chem.MolFromSmiles('C1=COC=C1'), 0, rdchem.BondType.SINGLE),  # Furan
    # (Chem.MolFromSmiles('C1=CSC=C1'), 0, rdchem.BondType.SINGLE),  # Thiophene
    # (Chem.MolFromSmiles('C1=CNC=C1'), 0, rdchem.BondType.SINGLE),  # Pyrrole
    # (Chem.MolFromSmiles('C1=CN=CN1'), 0, rdchem.BondType.SINGLE),  # Imidazole
    # (Chem.MolFromSmiles('C1=CON=C1'), 0, rdchem.BondType.SINGLE),  # Isoxazole
    # (Chem.MolFromSmiles('COc1cccs1'), 5, rdchem.BondType.SINGLE),  # Methoxy-substituted thiophene
    # (Chem.MolFromSmiles('NCc1cccs1'), 5, rdchem.BondType.SINGLE),  # Amino-substituted thiophene
    # (Chem.MolFromSmiles('c1cc[nH]c1'), 2, rdchem.BondType.SINGLE),  # Pyrrole
    # (Chem.MolFromSmiles('c1cc[nH]c1'), 1, rdchem.BondType.SINGLE),  # Pyrrole
    # (Chem.MolFromSmiles('Cc1cc(C)[nH]n1'), 2, rdchem.BondType.SINGLE),  # Methyl-substituted pyrazole
    (Chem.MolFromSmiles('c1cc(C)on1'), 2, rdchem.BondType.SINGLE),  # Methyl-substituted isoxazole
]

electron_withdrawing_rings = [
    (Chem.MolFromSmiles('C1=CC=NC=C1'), 0, rdchem.BondType.SINGLE),  # Pyridine
    # (Chem.MolFromSmiles('C1=CNN=C1'), 0, rdchem.BondType.SINGLE),  # Pyrazole
]

steric_blocking_substituents = [
    # (Chem.MolFromSmiles('C'), 0, rdchem.BondType.SINGLE),  # Methyl
    # (Chem.MolFromSmiles('CC'), 0, rdchem.BondType.SINGLE),  # Ethyl
    # (Chem.MolFromSmiles('C(C)C'), 0, rdchem.BondType.SINGLE),  # Isopropyl
    (Chem.MolFromSmiles('C(C)(C)C'), 0, rdchem.BondType.SINGLE)  # Tert-butyl
]
# # Validate substituent molecules
# invalid_substituent_EWG = validate_molecules(electron_withdraw_substituents)
# print("Invalid EWG substituents:", invalid_substituent_EWG)

# # Validate substituent molecules
# invalid_substituent_EDG = validate_molecules(electron_donating_substituents)
# print("Invalid EDG substituents:", invalid_substituent_EDG)

# # Validate ring substituent molecules
# invalid_ring_sub_EDG = validate_molecules(electron_donating_rings)
# print("Invalid EDG rings:", invalid_ring_sub_EDG)

# # Validate ring substituent molecules
# invalid_ring_sub_EWG = validate_molecules(electron_withdrawing_rings)
# print("Invalid EWG rings:", invalid_ring_sub_EWG)

# Generate azobenzene variations
azobenzene = Chem.MolFromSmiles('C1=CC=C(C=C1)N=NC2=CC=CC=C2')
all_positions = [0, 1, 2, 4, 5, 9, 10, 11, 12, 13]

ortho_positions = [2, 4, 9, 13]
meta_positions = [1, 5, 10, 12]
para_positions = [0, 11]

left_ring = [0, 1, 2, 4, 5]
left_ring_meta = [1, 5]
left_ring_ortho = [2, 4]
left_ring_para = [0]

right_ring = [9, 10, 11, 12, 13]
right_ring_meta = [10, 12]
right_ring_ortho = [9, 13]
right_ring_para = [11]

all_new_molecules = set()
max_subs_para = 1
max_subs_ortho = 2
max_subs_meta = 2

logging.info('Generating push-pull azobenzene variations')

# Generate molecules with electron-donating substituents on the left ring
new_azobenzenes_meta = generate_and_print_molecules(azobenzene, electron_donating_substituents, left_ring_meta, max_subs_meta, 'left', 'meta')
new_azobenzenes_para = generate_and_print_molecules(azobenzene, electron_donating_substituents, left_ring_para, max_subs_para, 'left', 'para')

# Generate molecules with electron-withdrawing substituents on the right ring
for new_azobenzene in new_azobenzenes_meta:
    all_new_molecules.update(generate_and_print_molecules(new_azobenzene, electron_withdrawing_substituents, right_ring_meta, max_subs_meta, 'right', 'meta'))
    all_new_molecules.update(generate_and_print_molecules(new_azobenzene, electron_withdrawing_substituents, right_ring_para, max_subs_para, 'right', 'para'))

for new_azobenzene in new_azobenzenes_para:
    all_new_molecules.update(generate_and_print_molecules(new_azobenzene, electron_withdrawing_substituents, right_ring_meta, max_subs_meta, 'right', 'meta'))
    all_new_molecules.update(generate_and_print_molecules(new_azobenzene, electron_withdrawing_substituents, right_ring_para, max_subs_para, 'right', 'para'))

logging.info('Generating azobenzene variations that hinder ts')

# Add substituents that destabilize the z isomer 
steric_azobenzenes = generate_and_print_molecules(azobenzene, steric_blocking_substituents, ortho_positions, max_subs_ortho, 'both', 'ortho')

for new_steric_azobenzene_right in steric_azobenzenes:
    more_steric_azobenzenes_right = generate_molecules(new_steric_azobenzene_right, electron_withdrawing_substituents, right_ring_para, 1)
    for num_subs in range(0,2):      
        for new_azobenzene_left in more_steric_azobenzenes_right:
            more_steric_azobenzenes = generate_molecules(new_azobenzene_left, electron_donating_substituents, left_ring_para, num_subs)
            all_new_molecules.update(more_steric_azobenzenes)

# Generate half-azobenzene variations with the ring substituents
half_azobenzene = Chem.MolFromSmiles('C1=CC=C(C=C1)N=N')
ring_sub_placement = [7]
ring_subs = electron_donating_rings + electron_withdrawing_rings

logging.info('Generating half-azobenzene variations with different rings')

half_azobenzenes_ring_subs = generate_molecules(half_azobenzene, ring_subs, ring_sub_placement, 1)

# Expand each of the molecules with the ring substituent by adding up to three additional substituents at the sub placements
ring_positions = [0, 1, 2, 4, 5]

logging.info("Adding subs to the half-azobenzene variations")

for max_subs in range(0, 4):
    # Generate molecules with electron-donating substituents
    for mol in half_azobenzenes_ring_subs:
        all_new_molecules.update(generate_and_print_molecules(mol, electron_donating_substituents, left_ring_meta, 2, 'half-azobenzene', 'electron-donating meta'))
        new_ortho = generate_and_print_molecules(mol, electron_donating_substituents, left_ring_ortho, 2, 'half-azobenzene', 'electron-donating ortho')
        all_new_molecules.update(new_ortho)
        for ortho_mol in new_ortho:
            added_para = generate_and_print_molecules(mol, electron_donating_substituents, left_ring_para, 1, 'half-azobenzene', 'electron-donating para')
            all_new_molecules.update(added_para)
            
    # Generate molecules with electron-withdrawing substituents
    for mol in half_azobenzenes_ring_subs:
        all_new_molecules.update(generate_and_print_molecules(mol, electron_withdrawing_substituents, left_ring_meta, 2, 'half-azobenzene', 'electron-withdrawing meta'))
        new_ortho = generate_and_print_molecules(mol, electron_withdrawing_substituents, left_ring_ortho, 2, 'half-azobenzene', 'electron-withdrawing ortho')
        all_new_molecules.update(new_ortho)
        for ortho_mol in new_ortho:
            added_para = generate_and_print_molecules(mol, electron_withdrawing_substituents, left_ring_para, 1, 'half-azobenzene', 'electron-withdrawing para')
            all_new_molecules.update(added_para)
# Convert all new molecules to SMILES strings
# all_new_molecules_smiles = {Chem.MolToSmiles(mol, canonical=True) for mol in all_new_molecules}
            
#
 #           all_new_molecules_smiles = {Chem.MolToSmiles(Chem.MolFromSmiles('Oc2ccc(/N=N/c1ccc(C(F)(F)F)cc1)cc2'), canonical=True),Chem.MolToSmiles(Chem.MolFromSmiles('N#Cc2ccc(/N=N/c1ccc(N)cc1)cc2'),canonical=True)}


easy_mol=['c1ccccc1N=Nc2occc2','c1ccccc1N=Nc2sccc2','c1ccccc1N=Nc2[nH]ccc2','c1ccccc1N=Nc2[nH]ncc2','c1ccccc1N=Nc2ncccc2','c1ccccc1N=Nc2nnccc2']

all_new_molecules_smiles = {Chem.MolToSmiles(Chem.MolFromSmiles(mol), canonical=True) for mol in easy_mol}

logging.info(f"Generated a total of {len(all_new_molecules_smiles)} new molecules.")


input_df = pd.DataFrame(columns=['comp_name', 'smiles', 'charge', 'multiplicity'])

for smiles in tqdm(all_new_molecules_smiles):
    mol = Chem.MolFromSmiles(smiles)
    comp_hash = get_hashed_label(smiles)
    comp_name = f"azo_{comp_hash}"
    charge = 0
    multiplicity = 1 

    input_df = pd.concat([input_df, pd.DataFrame([{'comp_name': comp_name, 'smiles': smiles, 'charge': charge, 'multiplicity': multiplicity}])], ignore_index=True)

input_df['molecular_weight'] = input_df['smiles'].apply(lambda x: Descriptors.MolWt(Chem.MolFromSmiles(x)))


# save to csv
input_df.to_csv('molecules.csv', index=False)
