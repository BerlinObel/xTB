import os
import sys
import re
import subprocess
import numpy as np
import matplotlib.pyplot as plt
from scipy.signal import gaussian
import pandas as pd
from rdkit import Chem
import copy

# Add the directory of the script to the Python path
sys.path.append(os.path.dirname(os.path.realpath(__file__)))
from settings import WAVELENGTH_RANGE, FACTOR, SIGMA_CM

from datetime import timedelta

def find_atom_mapping_reactant_to_product(reactant, product):
    """Find and return new atom order in product based on reactant"""
    # Copy molecules to prevent mutation of original molecules
    reactant_copy = copy.deepcopy(reactant)
    product_copy = copy.deepcopy(product)

    # Kekulize molecules to make aromaticity explicit
    Chem.Kekulize(reactant_copy, clearAromaticFlags=True)
    Chem.Kekulize(product_copy, clearAromaticFlags=True)

    # Normalize molecules to only compare connectivity
    reactant_norm = normalize_molecule(reactant_copy)
    product_norm = normalize_molecule(product_copy)

    # Attempt to match structure of product to reactant by breaking a bond
    smarts_bond = Chem.MolFromSmarts('[CX4;H0;R]-[CX4;H1;R]')
    atom_indices = list(reactant_norm.GetSubstructMatch(smarts_bond))

    if len(atom_indices) != 0:
        bond = reactant_norm.GetBondBetweenAtoms(
            atom_indices[0], atom_indices[1])
        broken_bond_reactant = Chem.FragmentOnBonds(
            reactant_norm, [bond.GetIdx()], addDummies=False)

        # find new atom order for product
        product_order = product_norm.GetSubstructMatch(broken_bond_reactant)
    else:
        product_order = product_norm.GetSubstructMatch(reactant_norm)

    return product_order


def normalize_molecule(mol):
    """Normalize molecule to compare only connectivity"""
    # Change all bond types to single
    for bond in mol.GetBonds():
        if bond.GetBondType() != Chem.BondType.SINGLE:
            bond.SetBondType(Chem.BondType.SINGLE)
    # Remove formal charges
    for atom in mol.GetAtoms():
        if atom.GetFormalCharge() != 0:
            atom.SetFormalCharge(0)
    return mol


def reorder_product_to_match_reactant(reactant, product):
    """Change atom order of product to match reactant"""
    new_product_order = find_atom_mapping_reactant_to_product(
        reactant, product)
    reordered_product = Chem.RenumberAtoms(product, new_product_order)
    return reordered_product

def format_time(seconds):

    delta = timedelta(seconds=seconds)
    days = delta.days
    hours, remainder = divmod(delta.seconds, 3600)
    minutes, seconds = divmod(remainder, 60)
    time_str = f"{days} days, {hours} hours, {minutes} minutes, {seconds} seconds"

    return time_str

def get_statistics(value_list):
    values = np.array(value_list)

    statistics_data = {
    "value_mean": np.mean(values),
    "value_median": np.median(values),
    "value_std": np.std(values),
    "value_variance": np.var(values),
    "value_range": np.ptp(values),
    "value_len": len(values)
    }

    statistics_dataframe = pd.DataFrame([statistics_data]) 
    return statistics_dataframe
    

def gaussian_distribution(wavelength_range, center_wavelength, oscillator_strength, sigma_cm):
    """
    This function creates a Gaussian distribution representing a single absorption peak. 
    The center of the Gaussian is the transition wavelength, the width is determined by sigma_cm 
    (the standard deviation), and the height is set by the oscillator strength.
    """
    distribution = oscillator_strength * np.exp(-4 * np.log(2) * ((1 / wavelength_range - 1 / center_wavelength) / (1e-7 * sigma_cm))**2)
    return distribution


def execute_shell_command(cmd, shell=False, timeout=None):
    """
    Executes a shell command and waits for it to finish execution.
    """
    if shell:
        p = subprocess.Popen(cmd, shell=True, stdin=subprocess.PIPE, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    else:
        cmd = cmd.split()
        p = subprocess.Popen(cmd, stdin=subprocess.PIPE, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        
    try:
        output, error = p.communicate(timeout=timeout)
        # if error:
            # print(f"Standard Error: {error.decode('utf-8')}")
        # if output:
            # print(f"{output.decode('utf-8')}")
    except subprocess.TimeoutExpired:
        # print(f"The command '{cmd}' timed out")
        p.kill()
        output = None
    except Exception as e:
        print(f"An error occurred while executing the command '{cmd}': {e}")
        output = 999

    return output

def run_orca(input_file, timeout=None):
    ORCA="/software/kemi/Orca/orca_5_0_1_linux_x86-64_openmpi411/orca"

    cmd = [ORCA, input_file]
    p = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE)

    output, error = p.communicate(timeout=timeout)
    if error:
        print(f"Standard Error: {error.decode('utf-8')}")
    if p.returncode != 0:
        print(f"ORCA exited with status {p.returncode}")
    # if output:
        # print(f"{output.decode('utf-8')}")
    return output, error

def orca_terminated_normally(output):
    normal_termination_pattern = r"\s*\*\*\*\*ORCA TERMINATED NORMALLY\*\*\*\*\s*"
    energy_match = re.search(normal_termination_pattern, output.decode('utf-8'))
    if energy_match:
        return True
    else:
        return False
    
def get_total_energy_xtb(output):

        if output == None:
            # print("Output is None")
            energy = None
        else:
            # The pattern to find the energy
            energy_pattern = r"\|\s*TOTAL ENERGY\s+(-?\d+\.\d+)\s+Eh\s*\|"
            energy_match = re.search(energy_pattern, output.decode('utf-8'))

            # Save the energy under results in the QMConf object
            energy = float(energy_match.group(1)) if energy_match else 0
        return energy
    
def get_final_energy_orca(output, type="out_file"):
    if output is None:
        return None
    
    # Decode the output if it's in bytes
    decoded_output = output.decode('utf-8') if isinstance(output, bytes) else output
    if type == "prop_file":
        orca_energy_pattern = r"\s*SCF Energy:\s+(-?\d+\.\d+)\s*"
    else:
        orca_energy_pattern = r"FINAL SINGLE POINT ENERGY\s+(-?\d+\.\d+)\s*"
    energy_matches = re.findall(orca_energy_pattern, decoded_output)
    
    # Get the last match
    return float(energy_matches[-1]) if energy_matches else None
    

def write_xyz(compound, name):
    filename = f"{name}_{compound.label}.xyz"
    compound.write_xyz(to_file=True)
    print(f"Written file: {filename}")