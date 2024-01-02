from turtle import forward
import pandas as pd
import numpy as np
import sys
import os
import re
import math
import time

from rdkit import Chem
from rdkit.Chem import Descriptors
from rdkit.Chem import AllChem

# Add the directory of the script to the Python path
sys.path.append(os.path.dirname(os.path.realpath(__file__)))
from db_utils import retrieve_data, update_data
from settings import QMC_PATH, REACTION_PATH_TEMPLATE
sys.path.append(QMC_PATH)

from utils import execute_shell_command, get_total_energy_xtb, format_time, run_orca

# from qmmol import QMMol
# from qmconf import QMConf
from QMC.qmmol import QMMol
from QMC.qmconf import QMConf

def get_tbr_reaction_path_search(output):
    # The pattern for finding the back reaction barrier
    tbr_pattern = "backward barrier \(kcal\)  :\s*([\d\.]+)"
    tbr_matches = re.findall(tbr_pattern, output.decode('utf-8'))

    forward_barrier_pattern = "forward  barrier \(kcal\)  :\s*([\d\.]+)"
    forward_matches = re.findall(forward_barrier_pattern, output.decode('utf-8'))
    reaction_energy_pattern = "reaction energy \(kcal\)  :\s*([\d\.]+)"
    reaction_matches = re.findall(reaction_energy_pattern, output.decode('utf-8'))

    # The last match is the last backward barrier
    tbr = float(tbr_matches[-1]) if tbr_matches else None
    forward_barrier = float(forward_matches[-1]) if forward_matches else None
    reaction_energy = float(reaction_energy_pattern[-1]) if reaction_matches else None
    
    return tbr, forward_barrier, reaction_energy    

def frequency_calculation(ts_xyz_file):
    
    frequency_calc_command = f"xtb {ts_xyz_file} --bhess"
    start_time_freq = time.time()
    frequency_output = execute_shell_command(frequency_calc_command, shell=False)
    print(frequency_output.decode('utf-8'))
    end_time_freq = time.time()
    elapsed_time_freq = end_time_freq - start_time_freq
    print(f"The frequency calc executed in {elapsed_time_freq} seconds")
    

def find_transition_state(reactant, product, ts_name):
    """Find transition state by using reaction path method from Grimme lab"""
    path_file = 'path.inp'
    output = None
    try:
        # Existing code for setting up reactant and product
        reactant_xyz_file = f"{reactant.label}.xyz"
        reactant.write_xyz(to_file=True)
        product_xyz_file = f"{product.label}.xyz"
        product.write_xyz(to_file=True)
    except Exception as e:
        print(f"Error in setting up reactant and product: {e}")

    print(f"Finding transition state for {reactant.label}")
    command_input = f"xtb {reactant_xyz_file} --path {product_xyz_file} --input {path_file} --gfn 1"

    try:
        output = execute_shell_command(command_input, shell=False, timeout=2000)
    except TimeoutError:
        print("Reaction path timed out")
    except Exception as exc:
        print(f"Error: {exc}")
        

    if output:
        # print(output.decode('utf-8'))

        ts_xyz_file = "xtbpath_ts.xyz"
        try:
            ts_object = QMConf(ts_xyz_file, fmt='xyz', label=ts_name)
        except Exception as exc:
            print(f"Error Creating TS Object: {exc}")
            return None, None, None
        
        energy = get_total_energy_xtb(output)
        ts_object.results['energy'] = energy

        tbr, forward_barrier, reaction_energy = get_tbr_reaction_path_search(output)
        info_dictionary = {'forward_barrier': (forward_barrier or 0)*4.184,
                           'reaction_energy': (reaction_energy or 0)*4.184
        }

        return ts_object, tbr, info_dictionary

    else:
        return None, None, None



def main(batch_id, reaction_path_settings):
    # Step 1: Load the data from the database
    filters = {'BatchID': batch_id, 'CalculationStage': 'storage_completed'}
    data = retrieve_data(filters)
    print(data)
    start_time = time.time()
    # The path.inp file controls the setting for the reaction path search
    path_file = "path.inp"

    # Create the path.inp file
    with open(path_file, 'w') as file:
        file.write(reaction_path_settings)   
         
    # Step 2: Perform the calculations and store the new data
    results = []
    for index, compound in data.iterrows():
        print(f"""
Index: {index+1}/{len(data)}, Compound: {compound.HashedName}
""")
        compound_name = str(compound.HashedName)  
        reactant = compound.ReactantObject 
        product = compound.ProductObject
        ts_name = f"{compound_name}_ts"
        ts_object, tbr, info_dictionary = find_transition_state(reactant, product, ts_name)
        
        if tbr is None or ts_object is None:
            print("No TS found!")

            results.append({
                'HashedName': str(compound.HashedName),
                'CalculationStage': 'ts_failed',
                })
        else:
            tbr = tbr * 4.184
            print(f"{ts_name} back-reaction barrier: {tbr} kj/mol")  
        
            results.append({
                'HashedName': str(compound.HashedName),
                'TransitionStateObject': ts_object,
                'BackReactionBarrier': tbr,
                'CalculationStage': 'ts_completed',
                'TransistionStateStats': info_dictionary
            })
        
        
    end_time = time.time()
    elapsed_time = end_time - start_time
    formatted_elapsed_time = format_time(elapsed_time)
    # Gather results
    results_df = pd.DataFrame(results)
    print(f"""
──────▄▀▄─────▄▀▄
─────▄█░░▀▀▀▀▀░░█▄
─▄▄──█░░░░░░░░░░░█──▄▄
█▄▄█─█░░▀░░┬░░▀░░█─█▄▄█
{batch_id} finished:
{results_df}

Time: {formatted_elapsed_time}
          """)
    
    # Update the database with the results
    update_data(results_df)


if __name__ == '__main__':
    batch_id = sys.argv[1]
    main(batch_id, REACTION_PATH_TEMPLATE)