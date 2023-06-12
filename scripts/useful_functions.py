import subprocess
import numpy as np

# Constants
AVOGADROS_NUMBER = 6.02214199 * 10**23
SPEED_OF_LIGHT = 299792458
ELECTRON_CHARGE = 1.60217662 * 10**(-19)
ELECTRON_MASS = 9.10938 * 10**(-31)
PERMITTIVITY_OF_VACUUM = 8.8541878176 * 10**(-12)
SIGMA_ELECTRON_VOLT = 0.4
SIGMA_CM = SIGMA_ELECTRON_VOLT * 8065.544
FACTOR = (AVOGADROS_NUMBER * ELECTRON_CHARGE**2) / (np.log(10) * 2 * ELECTRON_MASS * SPEED_OF_LIGHT**2 * PERMITTIVITY_OF_VACUUM) * np.sqrt(np.log(2)/np.pi) * 10**(-1)

# Define wavelength range
wavelength_range = np.linspace(125, 525, 593)

# Function to calculate absorption spectrum
def calculate_absorption_and_max_os_strength(wavelengths, wavelength_range, oscillator_strengths):
    absorption_spectrum = np.zeros(len(wavelength_range))
    max_os_strength = np.zeros(len(wavelength_range))

    for x in range(1, len(wavelength_range)):
        absorption_contributions = np.zeros(len(wavelengths))
        for i in range(len(wavelengths)):
            absorption_contributions[i] = (FACTOR / SIGMA_CM) * oscillator_strengths[i] * np.exp(-4 * np.log(2) * ((1 / wavelength_range[x] - 1 / wavelengths[i]) / (10**(-7) * SIGMA_CM))**2)
        absorption_spectrum[x] = absorption_contributions.sum()
        max_os_strength[x] = oscillator_strengths[absorption_contributions.argmax()]
    
    return absorption_spectrum, max_os_strength

def get_max_absorption(wavelengths, oscillator_strengths):
    absorption_spectrum, max_os_strength = calculate_absorption_and_max_os_strength(wavelengths, wavelength_range, oscillator_strengths)
    max_index = np.argmax(absorption_spectrum[185:500]) + 185
    max_wavelength = wavelength_range[max_index]
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