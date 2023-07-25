import os
import sys
import subprocess
import numpy as np
from scipy.signal import gaussian

script_dir = os.path.dirname(os.path.abspath(__file__))
sys.path.append(script_dir)
from settings import WAVELENGTH_RANGE, FACTOR, SIGMA_CM


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
def execute_shell_command(cmd, shell=False):
    """
    Executes a shell command and waits for it to finish execution.
    """
    if shell:
        p = subprocess.Popen(cmd, shell=True, stdin=subprocess.PIPE, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    else:
        cmd = cmd.split()
        p = subprocess.Popen(cmd, stdin=subprocess.PIPE, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    p.wait()
    output, _ = p.communicate()
    return output