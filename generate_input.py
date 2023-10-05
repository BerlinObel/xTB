import os
import sys
import pandas as pd
sys.path.append("/groups/kemi/brq616/speciale/opt/xTB/tQMC/QMC")

from qmmol import QMMol
from qmconf import QMConf

from tqdm import tqdm
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors as rdMD
from rdkit.Chem import Descriptors
import logging



easy_mol=['c1ccccc1N=Nc2occc2','c1ccccc1N=Nc2sccc2','c1ccccc1N=Nc2[nH]ccc2','c1ccccc1N=Nc2[nH]ncc2','c1ccccc1N=Nc2ncccc2','CCCCc1ccccc1/N=N/c2ccccc2']

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
