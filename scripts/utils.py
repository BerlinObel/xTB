import os
import sys
import re
import subprocess
import numpy as np
import matplotlib.pyplot as plt
from scipy.signal import gaussian
import pandas as pd

# Add the directory of the script to the Python path
sys.path.append(os.path.dirname(os.path.realpath(__file__)))
from settings import WAVELENGTH_RANGE, FACTOR, SIGMA_CM

from datetime import timedelta

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

def calculate_absorption_and_max_os_strength(wavelengths, oscillator_strengths):
    """
    This function calculates the overall absorption spectrum by summing the contributions from 
    Gaussians centered at each transition wavelength.
    """
    absorption_spectrum = np.zeros_like(WAVELENGTH_RANGE)
    max_os_strength = np.zeros_like(WAVELENGTH_RANGE)
    
    for i, wavelength in enumerate(wavelengths):
        contribution = gaussian_distribution(WAVELENGTH_RANGE, wavelength, oscillator_strengths[i], SIGMA_CM)
        absorption_spectrum += contribution
        max_os_strength[np.argmax(contribution)] = oscillator_strengths[i]

    absorption_spectrum *= FACTOR / SIGMA_CM

    # plt.plot(WAVELENGTH_RANGE, absorption_spectrum)

    return absorption_spectrum, max_os_strength

def get_max_absorption(wavelengths, oscillator_strengths):
    """
    This function finds the maximum absorption wavelength, its value, and the corresponding oscillator strength. 
    It specifically searches within the range of 185th to 500th wavelength assuming that the range covers the visible light spectrum.
    """
    absorption_spectrum, max_os_strength = calculate_absorption_and_max_os_strength(wavelengths, oscillator_strengths)
    max_index = np.argmax(absorption_spectrum[185:500]) + 185
    max_wavelength = WAVELENGTH_RANGE[max_index]
    max_absorption = absorption_spectrum[max_index]
    corresponding_max_os_strength = max_os_strength[max_index]

    return max_wavelength, max_absorption, corresponding_max_os_strength

# Execute programs in the shell
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
        p.wait(timeout=timeout)
        output, _ = p.communicate()
    except subprocess.TimeoutExpired:
        # print(f"The command '{cmd}' timed out")
        p.kill()
        output = None
    except Exception as e:
        print(f"An error occurred while executing the command '{cmd}': {e}")
        output = 999

    return output

def get_total_energy_xtb(output):

        if output == None:
            # print("Output is None")
            energy = None
        else:
            # The pattern to find the energy
            energy_pattern = r"\|\s*TOTAL ENERGY\s+(-?\d+\.\d+)\s+Eh\s*\|"
            energy_match = re.search(energy_pattern, output.decode('utf-8'))

            # Save the energy under results in the QMConf object
            energy = float(energy_match.group(1)) if energy_match else None
        return energy

def write_xyz(compound, name):
    filename = f"{name}_{compound.label}.xyz"
    compound.write_xyz(to_file=True)
    print(f"Written file: {filename}")