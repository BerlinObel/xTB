import pandas as pd
import sys
import os
import re
import time


import calc_sce
# Add the directory of the script to the Python path
sys.path.append(os.path.dirname(os.path.realpath(__file__)))
from utils import execute_shell_command, get_statistics, format_time
from settings import XTB4STDA_PATH, STDA_PATH, QMC_PATH
from db_utils import retrieve_data, update_data
sys.path.append(QMC_PATH)

def repeat_stda_calc(compound, n_runs):
    max_abs_reactant_lst = []
    osc_str_reactant_lst = []
    max_abs_product_lst = []
    osc_str_product_lst = []
    while len(max_abs_product_lst) < n_runs:
        max_abs_reactant_value, osc_str_reactant_value = calculate_absorbtion(compound.ReactantObject, print_results=False)
        max_abs_product_value, osc_str_product_value = calculate_absorbtion(compound.ProductObject, print_results=False)
        max_abs_reactant_lst.append(max_abs_reactant_value)
        osc_str_reactant_lst.append(osc_str_reactant_value)
        max_abs_product_lst.append(max_abs_product_value)
        osc_str_product_lst.append(osc_str_product_value)
    return max_abs_reactant_lst, osc_str_reactant_lst, max_abs_product_lst, osc_str_product_lst

def find_excited_states(molecule_file_path, charge, spin, print_command=True):
    """
    Calculates excited states for a given molecule using stda return the wavelenghts(nm) and oscillator strengths.
    """
    # Runs a ground state xtb calculation that writes a wfn.xtb with coordinates that is used in the stda step.
    shell_command = f"{XTB4STDA_PATH} {molecule_file_path} -chrg {charge} -uhf {spin}"
    if print_command:
        print(f"Shell Command: {shell_command}")
    execute_shell_command(shell_command, shell=False)

    # Executes the stda calculation for the excitation energies
    # -e 10 means that all excited stated up to 10 eV are calculated
    output = execute_shell_command(f'{STDA_PATH} -xtb -e 10', shell=False)

    # Shows the calculation outputs
    output_string = output.decode('utf-8')
    # print(output_string)

    # Regex patters that picks up the excitation energies, oscillator strengths and so on..
    data_pattern = re.compile(r"\s+(\d+)\s+([\d.]+)\s+([\d.]+)\s+([\d.]+)")

    wavelengths = []
    osc_strengths = []

    start_reading = False
    for line in output_string.split('\n'):
        # Ignore everything before "state    eV      nm       fL        Rv(corr)" line
        if 'state    eV      nm       fL        Rv(corr)' in line:
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
            # match[2] is the energy in eV
            # match[3] is the wavelength in nm
            # match[4] is the oscillator strength (fL)
            wavelengths.append(float(match.group(3)))
            osc_strengths.append(float(match.group(4)))

    return wavelengths, osc_strengths


def calculate_absorbtion(compound, print_results=True):
    charge = '0'
    spin = '0'

    # Create xyz file to use as input in the xtb/stda calculations
    xyz_file_path = str(compound.label)+'.xyz'
    compound.write_xyz(to_file=True)
    if print_results:
        print_command_line = True
    else:
        print_command_line = False
        
    # The results of the calculation
    wavelengths, osc_strengths = find_excited_states(
        xyz_file_path, charge, spin, print_command_line)

    # Find index of oscillator strength at the first peak and corresponding wavelength
    max_index = 1
    while True:
        if osc_strengths[max_index] > osc_strengths[max_index+1]:
            break
        else:
            max_index = max_index +1
    max_osc, max_wavelength = osc_strengths[max_index], wavelengths[max_index]
    if print_results:   
        print(f"Max oscillator strength: {max_osc}")
        print(f"Corresponding wavelength: {max_wavelength}")

    # os.remove(xyz_file_path)

    return max_wavelength, max_osc

def main(batch_id):
    # Step 1: Load the data from the database
    filters = {'BatchID': batch_id, 'CalculationStage': '%ts_complete%'}
    data = retrieve_data(filters)
    print(data)
    start_time = time.time()
    # Step 2: Perform the calculations and store the new data
    results = []
    for index, compound in data.iterrows():
        print(f"""
Index: {index+1}/{len(data)}, Compound: {compound.HashedName}
""")
        stat_dictionary = {}
        print(f"""
    Reactant
""")
        max_abs_reactant, osc_str_reactant = calculate_absorbtion(
            compound.ReactantObject)
        print(f"""
    Product
""")
        max_abs_product, osc_str_product = calculate_absorbtion(compound.ProductObject)

        max_abs_reactant_lst, osc_str_reactant_lst, max_abs_product_lst, osc_str_product_lst = repeat_stda_calc(compound, 10)
        stat_dictionary["Reactant_Max_Abs"] = get_statistics(max_abs_reactant_lst)
        stat_dictionary["Reactant_Osc"] = get_statistics(osc_str_reactant_lst)
        stat_dictionary["Product_Max_Abs"] = get_statistics(max_abs_product_lst)
        stat_dictionary["Product_Osc"] = get_statistics(osc_str_product_lst)  
    
        #Calculate Solar conversion efficiency
        solarCE = calc_sce.calculate_SCE(compound.StorageEnergy, compound.BackReactionBarrier, max_abs_reactant, max_abs_product, osc_str_reactant,osc_str_product)
        print(f"SCE: {solarCE}")
        
        results.append({
            'HashedName': str(compound.HashedName),
            'BatchID': str(compound.BatchID),
            'MaxAbsorptionProduct': max_abs_product,
            'MaxOscillatorStrengthProduct': osc_str_product,
            'MaxAbsorptionReactant': max_abs_reactant,
            'MaxOscillatorStrengthReactant': osc_str_reactant,
            'CalculationStage': 'abs_completed',
            'SolarConversionEfficiency': solarCE,
            'ExcitationStats': stat_dictionary
        })

    end_time = time.time()
    elapsed_time = end_time - start_time
    formatted_elapsed_time = format_time(elapsed_time)
    # Gather results
    results_df = pd.DataFrame(results)
    print(f"""
───▄▀▀▀▄▄▄▄▄▄▄▀▀▀▄───
───█▒▒░░░░░░░░░▒▒█───
────█░░█░░░░░█░░█────
─▄▄──█░░░▀█▀░░░█──▄▄─
█░░█─▀▄░░░░░░░▄▀─█░░█
{batch_id} finished:
{results_df}

Time: {formatted_elapsed_time}""")

    # Update the database with the results
    update_data(results_df)


if __name__ == '__main__':
    batch_id = sys.argv[1]
    main(batch_id)