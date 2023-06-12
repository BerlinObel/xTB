import pandas as pd
import os
import sys
import re

sys.path.append("/groups/kemi/brq616/speciale/opt/xTB/tQMC/QMC")
sys.path.append("/groups/kemi/brq616/speciale/opt/xTB/scripts")

from useful_functions import execute_shell_command
from useful_functions import get_max_absorption


def find_excited_states(molecule_file_path, charge, spin, n):
    """
    Calculates excited states for a given molecule using stda return the wavelenghts(nm) and oscillator strengths.
    """
    
    print(f"{molecule_file_path}, {charge}, {spin}, {n}")
    print("/groups/kemi/brq616/speciale/opt/stda/xtb4stda "+str(molecule_file_path)+" -chrg "+str(charge)+" -uhf "+str(spin))
    execute_shell_command(f"/groups/kemi/brq616/speciale/opt/stda/xtb4stda {molecule_file_path} -chrg {charge} -uhf {spin}",shell=False)
    output = execute_shell_command('/groups/kemi/brq616/speciale/opt/stda/stda_v1.6.1 -xtb -e 10',shell=False)

    output_string = output.decode('utf-8')
    print(output_string)

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
        # After the starting line, apply regex
        match = data_pattern.match(line)
        if match:
            # match[2] is the energy in eV
            # match[3] is the wavelength in nm
            # match[4] is the oscillator strength (fL)
            wavelengths.append(float(match.group(3)))
            osc_strengths.append(float(match.group(4)))
    
    return wavelengths, osc_strengths


def calculate_absorbtion(compound):
    # how to implement to control number of excited states?
    n = 75 #select the number of excited states
    charge = '0'
    spin = '0'

    xyz_file_path = str(compound.label)+'.xyz'
    compound.write_xyz(to_file=True)
    wavelengths, osc_strengths = find_excited_states(xyz_file_path, charge, spin, n)

    # Get max absorption
    max_wavelength, max_absorption, corresponding_max_os_strength = get_max_absorption(wavelengths, osc_strengths)

    print(f"Max absorption wavelength: {max_wavelength}")
    print(f"Corresponding absorption value: {max_absorption}")
    print(f"Corresponding oscillator strength: {corresponding_max_os_strength}")

    os.remove(xyz_file_path)
    
    return max_wavelength, corresponding_max_os_strength

if __name__ == '__main__':
    data_df = pd.read_pickle(sys.argv[1])
    compound_list = list()

    for compound in data_df.itertuples():
        reac_qmconf = compound.reac
        prod_qmconf = compound.prod
        ts_qmconf = compound.ts
        storage = compound.storage
        tbr = compound.tbr
        max_abs_reactant, osc_str_reactant = calculate_absorbtion(compound.reac)
        max_abs_product, osc_str_product = calculate_absorbtion(compound.prod)

        compound_list.append({'rep': compound.rep,
                              'reac': reac_qmconf,
                              'prod': prod_qmconf,
                              'ts': ts_qmconf,
                              'storage': storage,
                              'tbr': tbr,
                              'max_abs_reac': max_abs_reactant,
                              'osc_str_reac': osc_str_reactant,
                              'max_abs_prod': max_abs_product,
                              'osc_str_prod': osc_str_product})
        

    results_df = pd.DataFrame(compound_list)
    pickle_name = sys.argv[1].split('/')[0].split('.')[0] + '_abs.pkl'
    print(f"pickle name: {pickle_name}")
    print(results_df)
    results_df.to_pickle(pickle_name)