import pandas as pd
import numpy as np
import sys
import re
from rdkit import Chem
from rdkit.Chem import Descriptors
from rdkit.Chem import AllChem

sys.path.append("/groups/kemi/brq616/speciale/opt/xTB/tQMC/QMC")
sys.path.append("/groups/kemi/brq616/speciale/opt/xTB/scripts")

from useful_functions import execute_shell_command
from qmmol import QMMol
from qmconf import QMConf

def find_transistion_state(molecule):
    """Find transistion state by using reaction path method from Grimme lab"""
    path_file = 'path.inp'

    # Creating the xyz files needed fot the xtb calculations
    reactant_xyz_file = f"{molecule.reac.label}.xyz"
    molecule.reac.write_xyz(to_file=True)
    product_xyz_file = f"{molecule.prod.label}.xyz"
    molecule.prod.write_xyz(to_file=True)
    
    print(f"Finding transistion state for {molecule.reac.label}")
    # The command takes a start xyz file (the reactant) and an end xyz file (the product)
    command_input = f"xtb {reactant_xyz_file} --path {product_xyz_file} --input {path_file}"

    # Start the reaction path calculation
    output = execute_shell_command(command_input, shell=False)
    print(output.decode('utf-8'))
   
    ts_xyz_file = "xtbpath_ts.xyz"
    ts_name  = f"{str(molecule.reac.label.split('_r')[0])}_ts" 
    ts_object = QMConf(ts_xyz_file, fmt='xyz', label=ts_name)
    
    # The pattern for finding the back reaction barrier
    pattern = "backward barrier \(kcal\)  :\s*([\d\.]+)"

    # Find all matches in the output string
    matches = re.findall(pattern, output.decode('utf-8'))

    # The last match is the last backward barrier
    tbr = float(matches[-1]) if matches else None

    return ts_object, tbr


def main(input_filename):
    # Load the data
    data = pd.read_pickle(input_filename)

    # The path.inp file controls the setting for the reaction path search
    path_file = "path.inp"
    reaction_path_settings = '''$path
    nrun=1
    npoint=25
    anopt=10
    kpush=0.003
    kpull=-0.015
    ppull=0.05
    alp=0.9
    $end'''

    # Create the path.inp file
    with open(path_file, 'w') as file:
        file.write(reaction_path_settings)    

    # Find transistions states for each compound     
    compound_list = []
    for compound in data.itertuples():

        compound_name = str(compound.reac.label.split('_r')[0])  
              
        ts_name = f"{compound_name}_ts"
        ts_object, tbr = find_transistion_state(compound)
        if tbr == None:
            print("No TS found!")
            tbr = 0
        
        tbr = tbr * 4.184
        print(f"{ts_name} back-reaction barrier: {tbr}")  
        
        compound_list.append({
            'rep': compound.rep,
            'reac': compound.reac,
            'prod': compound.prod,
            'ts': ts_object,
            'storage': compound.storage * 2625.5,
            'tbr': tbr
        })

    # Save the results
    results_df = pd.DataFrame(compound_list)
    pickle_name = f"{input_filename.split('/')[0].split('.')[0]}_new.pkl"
    print(f"pickle name: {pickle_name}")
    print(results_df)
    results_df.to_pickle(pickle_name)


if __name__ == "__main__":
    main(sys.argv[1])