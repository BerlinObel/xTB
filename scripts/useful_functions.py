import subprocess
import numpy as np
from scipy.signal import gaussian

# Constants
AVOGADROS_NUMBER = 6.02214199e23
SPEED_OF_LIGHT = 299792458
ELECTRON_CHARGE = 1.60217662e-19
ELECTRON_MASS = 9.10938e-31
PERMITTIVITY_OF_VACUUM = 8.8541878176e-12
SIGMA_ELECTRON_VOLT = 0.4
SIGMA_CM = SIGMA_ELECTRON_VOLT * 8065.544
FACTOR = (AVOGADROS_NUMBER * ELECTRON_CHARGE**2) / (np.log(10) * 2 * ELECTRON_MASS * SPEED_OF_LIGHT**2 * PERMITTIVITY_OF_VACUUM) * np.sqrt(np.log(2)/np.pi) * 1e-1

# Define wavelength range
WAVELENGTH_RANGE = np.linspace(125, 525, 593)

def gaussian_distribution(wavelength_range, center_wavelength, oscillator_strength, sigma_cm):
    """Create a gaussian distribution."""
    distribution = oscillator_strength * np.exp(-4 * np.log(2) * ((1 / wavelength_range - 1 / center_wavelength) / (1e-7 * sigma_cm))**2)
    return distribution

def calculate_absorption_and_max_os_strength(wavelengths, oscillator_strengths):
    """Calculate the absorption spectrum and max oscillator strength."""
    absorption_spectrum = np.zeros_like(WAVELENGTH_RANGE)
    max_os_strength = np.zeros_like(WAVELENGTH_RANGE)

    for i, wavelength in enumerate(wavelengths):
        contribution = gaussian_distribution(WAVELENGTH_RANGE, wavelength, oscillator_strengths[i], SIGMA_CM)
        absorption_spectrum += contribution
        max_os_strength[np.argmax(contribution)] = oscillator_strengths[i]

    absorption_spectrum *= FACTOR / SIGMA_CM

    return absorption_spectrum, max_os_strength

def get_max_absorption(wavelengths, oscillator_strengths):
    """Find the max absorption wavelength, value and corresponding oscillator strength."""
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