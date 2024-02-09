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
from db_utils import retrieve_data, update_data
from settings import QMC_PATH
sys.path.append(QMC_PATH)
import calc_sce

from utils import format_time, run_orca, orca_terminated_normally

# from qmmol import QMMol
# from qmconf import QMConf
from QMC.qmmol import QMMol
from QMC.qmconf import QMConf


def write_orca_td_dft_input(xyz_file, input_file): 
   
    input_string = f"""! PAL2
! RKS cc-pVDZ M062X 
! NEB-TS

%TDDFT
   NROOTS   30
END
* xyzfile 0 1 {xyz_file} 
    """
    
    print(f"input orca: {input_string}")
    with open(input_file, 'w') as file:
            file.write(input_string) 
            
def find_max_absorption_wavelength(output):
    # Shows the calculation outputs
    output_string = output.decode('utf-8')

    # Regex patters that picks up the excitation energies, oscillator strengths and so on..
    data_pattern = re.compile(r"\s+(\d+)\s+([-?\d.]+)\s+([-?\d.]+)\s+([-?\d.]+)\s+([-?\d.]+)\s+([-?\d.]+)\s+([-?\d.]+)\s+([-?\d.]+)")

    wavelengths = []
    osc_strengths = []

    start_reading = False
    for line in output_string.split('\n'):
        # Ignore everything before "ABSORPTION SPECTRUM VIA TRANSITION ELECTRIC DIPOLE MOMENTS" line
        if 'ABSORPTION SPECTRUM VIA TRANSITION ELECTRIC DIPOLE MOMENTS' in line:
            start_reading = True
            continue
        if not start_reading:
            continue
        # If there is an empty line, stop reading
        if line.strip() == '':
            break
        # After the starting line, apply regex
        match = data_pattern.match(line)
        if match:
            # match[0] is state
            # match[1] is the energy in cm-1
            # match[2] is the wavelength in nm
            # match[3] is the oscillator strength (fosc)
            # rest is transition dipole moments
            wavelengths.append(float(match.group(2)))
            osc_strengths.append(float(match.group(3)))
            
    max_index = osc_strengths.index(max(osc_strengths[1:]))
    max_osc, max_wavelength = osc_strengths[max_index], wavelengths[max_index]

    return max_osc, max_wavelength    

def run_single_td_dft_calculation(xyz_file, compound_name, label):
    orca_input = f"{compound_name}_{label}_td_dft.inp"
    write_orca_td_dft_input(xyz_file, orca_input)
    
    start_time_orca = time.time()
    orca_output, orca_error = run_orca(orca_input)
    normal_termination = orca_terminated_normally(orca_output)
    end_time_orca = time.time()
    elapsed_time_orca = end_time_orca - start_time_orca
    max_osc, max_wavelength = find_max_absorption_wavelength(orca_output)
    print(f"Orca Time: {format_time(elapsed_time_orca)}")
    
    return normal_termination, max_osc, max_wavelength


def run_td_dft_calculations(reactant_xyz_file, product_xyz_file, compound_name):
    
    try:
        normal_termination_reactant, max_osc_reactant, max_wavelength_reactant = run_single_td_dft_calculation(
            reactant_xyz_file, compound_name, "r"
        )
    except Exception as e:
        print(f"Error running Geometry Optimization for reactant: {e}")
        normal_termination_reactant, max_osc_reactant, max_wavelength_reactant = None, None, None
    
    try:
        normal_termination_product, max_osc_product, max_wavelength_product = run_single_td_dft_calculation(
            product_xyz_file, compound_name, "p"
        )
    except Exception as e:
        print(f"Error running Geometry Optimization for product: {e}")
        normal_termination_product, max_osc_product, max_wavelength_product = None, None, None

    return normal_termination_reactant, max_osc_reactant, max_wavelength_reactant, normal_termination_product, max_osc_product, max_wavelength_product

def main(batch_id):
    # Step 1: Load the data from the database
    filters = {'BatchID': batch_id, 'CalculationStage': 'storage_completed'}
    data = retrieve_data(filters, "ORCAData")
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
        reactant_name = f"{reactant}_td_dft"
        product_name = f"{product}_td_dft"
        
        molecule_dir = os.path.join('/groups/kemi/brq616/speciale/dft_calc', compound_name)
        reactant_dir = os.path.join(molecule_dir, reactant)
        td_dft_reactant_dir = os.path.join(reactant_dir, "td_dft")
        os.makedirs(td_dft_reactant_dir, exist_ok=True)
        
        product_dir = os.path.join(molecule_dir, product)
        td_dft_product_dir = os.path.join(product_dir, "td_dft")
        os.makedirs(td_dft_product_dir, exist_ok=True)

        reactant_xyz_file = glob.glob(os.path.join(reactant_dir, "*.xyz"))

        # Get the absolute path of the single .xyz file
        if len(reactant_xyz_file) == 1:
            reactant_xyz_path = os.path.abspath(reactant_xyz_file[0])
        else:
            print("No .xyz file found or more than one reactant.xyz files found.")
            
        product_xyz_file = glob.glob(os.path.join(product_dir, "*.xyz"))

        # Get the absolute path of the single .xyz file
        if len(product_xyz_file) == 1:
            product_xyz_path = os.path.abspath(product_xyz_file[0])
        else:
            print("No .xyz file found or more than one product.xyz files found.")
   
        normal_termination_reactant, max_osc_reactant, max_wavelength_reactant, normal_termination_product, max_osc_product, max_wavelength_product = run_td_dft_calculations(reactant_xyz_path, product_xyz_path)

        solarCE = calc_sce.calculate_SCE(compound.StorageEnergy, compound.BackReactionBarrier,max_wavelength_reactant,max_wavelength_product, max_osc_reactant,max_osc_product)
        if normal_termination_reactant and normal_termination_product:
            results.append({
                'HashedName': compound_name,
                'BatchID': batch_id,
                'CalculationStage': 'abs_completed',
                'MaxAbsorptionProduct': max_wavelength_product,
                'MaxOscillatorStrengthProduct': max_osc_product,
                'MaxAbsorptionReactant': max_wavelength_reactant,
                'MaxOscillatorStrengthReactant': max_osc_reactant,
                'SolarConversionEfficiency': solarCE
                })
            
        else:
            results.append({
                'HashedName': compound_name,
                'BatchID': batch_id
                })

        source_directory = "./"  # Current directory

        for filename in glob.glob(f"{source_directory}*{reactant_name}*"):
            shutil.move(filename, td_dft_reactant_dir)
        for filename in glob.glob(f"{source_directory}*{product_name}*"):
            shutil.move(filename, td_dft_product_dir)
          
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
    update_data(results_df, "ORCAData")

if __name__ == '__main__':
    batch_id = sys.argv[1]
    main(batch_id)