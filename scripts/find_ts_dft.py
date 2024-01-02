import pandas as pd
import numpy as np
import sys
import os
import re
import math
import time
import shutil
import glob
from pathlib import Path

from rdkit import Chem
from rdkit.Chem import Descriptors
from rdkit.Chem import AllChem
from rdkit.Chem import rdmolfiles

# Add the directory of the script to the Python path
sys.path.append(os.path.dirname(os.path.realpath(__file__)))
from db_utils import retrieve_data, insert_data
from settings import QMC_PATH
sys.path.append(QMC_PATH)

from utils import format_time, run_orca, orca_terminated_normally, get_final_energy_orca

# from qmmol import QMMol
# from qmconf import QMConf
from QMC.qmmol import QMMol
from QMC.qmconf import QMConf


def write_orca_neb_input(reactant_xyz_file, product_xyz_file, input_file): 
   
    input_string = f"""! PAL2
! RKS cc-pVDZ M062X 
! NEB-TS

%neb
 neb_end_xyzfile "{product_xyz_file}"
 PreOpt_Ends true
 NImages 16
 Opt_Method VPO
end

%maxcore 4000
* xyzfile 0 1 {reactant_xyz_file} 
    """
    
    print(f"input orca: {input_string}")
    with open(input_file, 'w') as file:
            file.write(input_string) 
            

def run_neb_ts_search(reactant_xyz_file, product_xyz_file, compund_name):

    try:
        orca_input = f"{compund_name}_ts.inp" 
        write_orca_neb_input(reactant_xyz_file, product_xyz_file, orca_input)
        
        start_time_orca = time.time()
        orca_output, orca_error = run_orca(orca_input)
        normal_termination = orca_terminated_normally(orca_output)
        end_time_orca = time.time()
        elapsed_time_orca = end_time_orca - start_time_orca
        ts_energy = get_final_energy_orca(orca_output)
        with open(f"{compund_name}_ts.out", 'w') as file:
            file.write(orca_output.decode("utf-8"))
        print(f"Orca Time: {format_time(elapsed_time_orca)}")
    except Exception as e:
        print(f"Error running Geometry Optimization for reactant: {e}")
        
    return normal_termination, ts_energy



def main(batch_id):
    # Step 1: Load the data from the database
    filters = {'BatchID': batch_id, 'CalculationStage': 'storage_completed'}
    data = retrieve_data(filters, 'ORCAData')
    print(data)
    start_time = time.time()
    results = []
    # Step 2: Perform the calculations and store the new data
    for index, compound in data.iterrows():
        print(f"""
Index: {index+1}/{len(data)}, Compound: {compound.HashedName}
""")
        compound_name = str(compound.HashedName)  
        reactant = f"{compound_name}_r" 
        product = f"{compound_name}_p"
        ts_name = f"{compound_name}_ts"
        storage_energy = compound.StorageEnergy
        
        current_working_directory = os.getcwd()
        base_dir = '/groups/kemi/brq616/speciale/dft_calc'
        molecule_dir = os.path.join(base_dir, compound_name)

        ts_dir = os.path.join(molecule_dir, ts_name)
        # Check if the directory exists
        if os.path.exists(ts_dir):
            # Clear the contents of the directory
            shutil.rmtree(ts_dir)

        os.makedirs(ts_dir, exist_ok=True)
        reactant_dir = os.path.join(molecule_dir, reactant)
        product_dir = os.path.join(molecule_dir, product)

        reactant_xyz_files = glob.glob(os.path.join(reactant_dir, "*_opt.xyz"))
        product_xyz_files = glob.glob(os.path.join(product_dir, "*_opt.xyz"))
        product_txt_files = glob.glob(os.path.join(product_dir, "*.txt"))

        if not product_txt_files:
            print("No .txt files found in the specified directory.")
            return

        product_output = Path(product_txt_files[0]).read_text()
        product_energy = get_final_energy_orca(product_output, type="prop_file")

        if len(reactant_xyz_files) != 1:
            print("No .xyz file found or more than one reactant.xyz files found.")
            return

        if len(product_xyz_files) != 1:
            print("No .xyz file found or more than one product.xyz files found.")
            return

        reactant_xyz_path = os.path.abspath(reactant_xyz_files[0])
        product_xyz_path = os.path.abspath(product_xyz_files[0])

        reactant_xyz_relative = os.path.relpath(reactant_xyz_path, start=current_working_directory)
        product_xyz_relative = os.path.relpath(product_xyz_path, start=current_working_directory)

        termination, ts_energy = run_neb_ts_search(reactant_xyz_relative, product_xyz_relative, compound_name)
    
        if termination:
            tbr = ts_energy - product_energy
            molecule_data = {
                'HashedName': compound_name,
                'BatchID': batch_id,
                'CalculationStage': 'ts_completed',
                'StorageEnergy': storage_energy,
                'BackReactionBarrier': tbr
     
                }
            insert_data(molecule_data, "ORCAData")
        else:
            molecule_data = {
                'HashedName': compound_name,
                'BatchID': batch_id,
                'CalculationStage': 'storage_completed',
                'StorageEnergy': storage_energy,
                }
            insert_data(molecule_data, "ORCAData")
            
        source_directory = "./" 
        for filename in glob.glob(f"{source_directory}*{ts_name}*"):
            shutil.move(filename, ts_dir)
          
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
    results_df = pd.DataFrame(results)

if __name__ == '__main__':
    batch_id = sys.argv[1]
    main(batch_id)