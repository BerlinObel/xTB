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

from utils import execute_shell_command, get_total_energy_xtb, format_time

# from qmmol import QMMol
# from qmconf import QMConf
from QMC.qmmol import QMMol
from QMC.qmconf import QMConf

def get_tbr_reaction_path_search(output):
    # The pattern for finding the back reaction barrier
    tbr_pattern = "backward barrier \(kcal\)  :\s*([\d\.]+)"
    tbr_matches = re.findall(tbr_pattern, output.decode('utf-8'))

    # The last match is the last backward barrier
    tbr = float(tbr_matches[-1]) if tbr_matches else None
    return tbr    

def check_transition_state(reactant_mol, product_mol, ts_mol):
    try:
        reactant_carbon_indices = find_carbon_connected_to_N_double(reactant_mol)
        reactant_carbon_distance = find_atom_distance(reactant_mol, reactant_carbon_indices)
        
        product_carbon_indices = find_carbon_connected_to_N_double(product_mol)
        product_carbon_distance = find_atom_distance(product_mol, product_carbon_indices)
        
        ts_carbon_indices = find_carbon_connected_to_N_double(ts_mol)
        ts_carbon_distance = find_atom_distance(ts_mol, ts_carbon_indices)
        
        print(f"Reactant distance: {reactant_carbon_distance:.3f}, TS distance: {ts_carbon_distance:.3f}, Product distance: {product_carbon_distance:.3f}")
        if reactant_carbon_distance > ts_carbon_distance > product_carbon_distance:
            return True
        else:
            return False
    except Exception as e:
        print(f"An error occurred in check_transistion_state: {e}")
        return False

def find_carbon_connected_to_N_double(mol):
    # Look for carbons connected to the nitrogens in the double bond
    nitrogen_indices = []
    carbon_indices = []
    try:
        # Iterate through atoms to find nitrogen atoms
        for atom in mol.GetAtoms():
            if atom.GetAtomicNum() == 7:  # Atomic number for Nitrogen
                for neighbor in atom.GetNeighbors():
                    if neighbor.GetAtomicNum() == 7:  # Atomic number for Nitrogen
                        # Check if the bond between them is double
                        bond = mol.GetBondBetweenAtoms(atom.GetIdx(), neighbor.GetIdx())
                        if bond.GetBondType() == Chem.rdchem.BondType.DOUBLE:
                            nitrogen_indices.append(atom.GetIdx())
                            break  # Stop checking other neighbors for this nitrogen
        # Find carbon atoms connected to the nitrogen atoms
        for idx in nitrogen_indices:
            atom = mol.GetAtomWithIdx(idx)
            for neighbor in atom.GetNeighbors():
                if neighbor.GetAtomicNum() == 6:  # Atomic number for Carbon
                    carbon_indices.append(neighbor.GetIdx()) 
        if not carbon_indices:
            print("No carbon indices found.")
        return carbon_indices
    except Exception as e:
        print(f"An error occurred in find_carbon_connected_to_N_double: {e}")
        return []

def find_atom_distance(rdKit_mol, atom_indices):
    # Get the conformer from the molecule
    conf = rdKit_mol.GetConformer()
    try:
        # Get the coordinates of the two atoms
        coord1 = conf.GetAtomPosition(atom_indices[0])
        coord2 = conf.GetAtomPosition(atom_indices[1])

        # Calculate the distance between the two atoms
        dist = math.sqrt((coord1.x - coord2.x)**2 + (coord1.y - coord2.y)**2 + (coord1.z - coord2.z)**2)
        if dist is None:
            print("Distance could not be calculated.")
        return dist
    except Exception as e:
        print(f"An error occurred in find_atom_distance: {e}")
        return None
    
def frequency_calculation(ts_xyz_file):
    
    frequency_calc_command = f"xtb {ts_xyz_file} --bhess"
    start_time_freq = time.time()
    frequency_output = execute_shell_command(frequency_calc_command, shell=False)
    print(frequency_output.decode('utf-8'))
    end_time_freq = time.time()
    elapsed_time_freq = end_time_freq - start_time_freq
    print(f"The frequency calc executed in {elapsed_time_freq} seconds")
    
def find_transistion_state(reactant, product, ts_name):
    """Find transistion state by using reaction path method from Grimme lab"""
    path_file = 'path.inp'

    # Creating the xyz files needed fot the xtb calculations
    reactant_xyz_file = f"{reactant.label}.xyz"
    reactant.write_xyz(to_file=True)
    product_xyz_file = f"{product.label}.xyz"
    product.write_xyz(to_file=True)
    
    print(f"Finding transistion state for {reactant.label}")
    # The command takes a start xyz file (the reactant) and an end xyz file (the product)
    command_input = f"xtb {reactant_xyz_file} --path {product_xyz_file} --input {path_file} --gfn 1"

    # Start the reaction path calculation
    output = execute_shell_command(command_input, shell=False)
    print(output.decode('utf-8'))
   
    ts_xyz_file = "xtbpath_ts.xyz"
    ts_object = QMConf(ts_xyz_file, fmt='xyz', label=ts_name)
    
    # Save the energy under results in the QMConf object
    energy = get_total_energy_xtb(output)
    ts_object.results['energy'] = energy

    # Back reaction barrier
    tbr = get_tbr_reaction_path_search(output)
    
    # frequency_calculation(ts_xyz_file)

    return ts_object, tbr


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
        ts_object, tbr = find_transistion_state(reactant, product, ts_name)
        
        if tbr == None or ts_object == None:
            print("No TS found!")

            results.append({
                'HashedName': str(compound.HashedName)
                })
        else:
            tbr = tbr * 4.184
            print(f"{ts_name} back-reaction barrier: {tbr}")  
        
            results.append({
                'HashedName': str(compound.HashedName),
                'BatchID': str(compound.BatchID),
                'TransitionStateObject': ts_object,
                'BackReactionBarrier': tbr,
                'CalculationStage': 'ts_completed'
            })
        
            ts_valid = check_transition_state(reactant.get_rdkit_mol(), product.get_rdkit_mol(), ts_object.get_rdkit_mol())
            if ts_valid:
                print("Carbon distance is inbetween product and reactant")
            else:
                print("Carbon distance is the same")
        
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