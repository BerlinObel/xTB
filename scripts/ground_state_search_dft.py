import pandas as pd
import numpy as np
import sys
import os
import re
import math
import time
import shutil
import glob

from rdkit import Chem
from rdkit.Chem import Descriptors
from rdkit.Chem import AllChem
from rdkit.Chem import rdmolfiles

# Add the directory of the script to the Python path
sys.path.append(os.path.dirname(os.path.realpath(__file__)))
from db_utils import retrieve_data, insert_data, create_orca_table
from settings import QMC_PATH
sys.path.append(QMC_PATH)

from utils import format_time, get_final_energy_orca, run_orca, reorder_product_to_match_reactant, orca_terminated_normally

# from qmmol import QMMol
# from qmconf import QMConf
from QMC.qmmol import QMMol
from QMC.qmconf import QMConf


def write_orca_geo_opt_input(xyz_file, input_file): 
    current_directory = os.getcwd()
    file_path = os.path.join(current_directory, xyz_file)
    
    input_string = f"""! PAL2
! RKS cc-pVDZ M062X TightSCF TightOPT
! Opt Freq
%method
Grid 7
end
%maxcore 4800
* xyzfile 0 1 {file_path} 
    """
    print(f"input orca: {input_string}")
    with open(input_file, 'w') as file:
            file.write(input_string) 
            
def write_xyz(mol, filename):
    xyz_data = rdmolfiles.MolToXYZBlock(mol)
    with open(filename, 'w') as f:
        f.write(xyz_data) 
        
def set_up_reactant_and_product(reactant, product):
    reactant_mol =  (reactant.get_rdkit_mol())
    product_mol = Chem.AddHs(product.get_rdkit_mol())
    product_mol = reorder_product_to_match_reactant(reactant_mol, product_mol)
    return reactant_mol, product_mol

def run_geometry_optimization(reactant, product):

    try:
        reactant_mol, product_mol = set_up_reactant_and_product(reactant, product)
        
        reactant_name = reactant.label
        reactant_xyz_file = f"{reactant_name}.xyz"
        product_name = product.label
        product_xyz_file = f"{product_name}.xyz"
        write_xyz(reactant_mol, reactant_xyz_file)
        write_xyz(product_mol, product_xyz_file)
    except Exception as e:
        print(f"Error in setting up reactant and product: {e}")
        
    try:
        orca_input_reactant = f"{reactant_name}_opt.inp" 
        write_orca_geo_opt_input(reactant_xyz_file, orca_input_reactant)
        
        start_time_orca = time.time()
        orca_output, orca_error = run_orca(orca_input_reactant)
        normal_termination_reactant = orca_terminated_normally(orca_output)
        end_time_orca = time.time()
        elapsed_time_orca = end_time_orca - start_time_orca
        reactant_enthalpy = get_final_energy_orca(orca_output)
        print(f"Orca Time: {format_time(elapsed_time_orca)}")
    except Exception as e:
        print(f"Error running Geometry Optimization for reactant: {e}")
        
    try:
        orca_input_product = f"{product_name}_opt.inp"
        write_orca_geo_opt_input(product_xyz_file, orca_input_product)
    
        start_time_orca = time.time()
        orca_output, orca_error = run_orca(orca_input_product)
        normal_termination_product = orca_terminated_normally(orca_output)
        end_time_orca = time.time()
        elapsed_time_orca = end_time_orca - start_time_orca
        product_enthalpy = get_final_energy_orca(orca_output)
        print(f"Orca Time: {format_time(elapsed_time_orca)}")
    except Exception as e:
        print(f"Error running Geometry Optimization for product: {e}")

    storage_energy = product_enthalpy-reactant_enthalpy
    return normal_termination_reactant, normal_termination_product, storage_energy
        

def main(batch_id):
    # Step 1: Load the data from the database
    filters = {'BatchID': batch_id, 'CalculationStage': 'storage_completed'}
    data = retrieve_data(filters)
    print(data)
    start_time = time.time()
    
    # Step 2: Perform the calculations and store the new data
    for index, compound in data.iterrows():
        print(f"""
Index: {index+1}/{len(data)}, Compound: {compound.HashedName}
""")
        compound_name = str(compound.HashedName)  
        reactant = compound.ReactantObject 
        product = compound.ProductObject
        reactant_termination, product_termination, storage_energy = run_geometry_optimization(reactant, product)


        if reactant_termination and product_termination:
            molecule_data = {
                'HashedName': compound_name,
                'BatchID': batch_id,
                'CalculationStage': 'storage_completed',
                'StorageEnergy': storage_energy
                }
            insert_data(molecule_data, "ORCAData")
        else:
            molecule_data = {
                'HashedName': compound_name,
                'BatchID': batch_id,
                'CalculationStage': 'storage'
                }
            insert_data(molecule_data, "ORCAData")
            
        molecule_dir = os.path.join('/groups/kemi/brq616/speciale/dft_calc', compound_name)
        os.makedirs(molecule_dir, exist_ok=True)
        reactant_dir = os.path.join(molecule_dir, reactant.label)
        os.makedirs(reactant_dir, exist_ok=True)
        product_dir = os.path.join(molecule_dir, product.label)
        os.makedirs(product_dir, exist_ok=True)
        
        source_directory = "./"  # Current directory

        for filename in glob.glob(f"{source_directory}*{reactant.label}*"):
            shutil.move(filename, reactant_dir)
        for filename in glob.glob(f"{source_directory}*{product.label}*"):
            shutil.move(filename, product_dir)
          
    end_time = time.time()
    elapsed_time = end_time - start_time
    formatted_elapsed_time = format_time(elapsed_time)
    print(f"""
       .
      ":"
    ___:____     |"\/"|
  ,'        `.    \  /
  |  O        \___/  |
~^~^~^~^~^~^~^~^~^~^~^~^~
{batch_id} finished:
Time: {formatted_elapsed_time}
          """)


if __name__ == '__main__':
    batch_id = sys.argv[1]
    main(batch_id)