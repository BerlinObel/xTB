
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import math




try:
    am15_sheet = pd.read_pickle('am15_SMARTS2.pkl')
except:
    am15_path = "/lustre/hpc/kemi/elholm/bod_calc/alle/eff/am_15g.xls"
    am15_sheet = pd.read_excel(am15_path,sheet_name = "SMARTS2",header=None)
    pd.to_pickle(am15_sheet,'am15_SMARTS2.pkl')


# Read solar spectrum - chosen spectrum normalized to 1000 (1000 W/m2)
solar_spectrum = am15_sheet.to_numpy()
wavelengths = solar_spectrum[2:, 0]
solar_irradiance = solar_spectrum[2:, 2]


# Natural constants
hbar = 1.05457182 * 10 ** (-34)
planck_constant = 6.6261 * 10 ** (-34)
speed_of_light = 299792458
avogadro_number = 6.0221408 * 10 ** (23)
electronvolt_to_joule = 1.602176634 * 10 ** (-19)

def calculate_total_solar_power(irradiance_spectrum, wavelength_spectrum):
    """
    Calculate the total solar power based on the solar irradiance spectrum and wavelength spectrum.
    :param irradiance_spectrum: Array of solar irradiance values (in W/m2).
    :param wavelength_spectrum: Array of solar wavelengths (in nm).
    :return: Total solar power (in W).
    """
    power_intervals = [
        irradiance * (next_wavelength - current_wavelength)
        for irradiance, current_wavelength, next_wavelength 
        in zip(irradiance_spectrum[1:], wavelength_spectrum[:-1], wavelength_spectrum[1:])
    ]
    total_power = irradiance_spectrum[0] * 0.5 + sum(power_intervals)
    return total_power

def calculate_photon_absorption_fraction(epsilon_reactant, epsilon_product, concentration_reactant=0.5, concentration_product=0.5):
    """
    Calculate the number of photons absorbed by the reactant.
    :param epsilon_reactant: Molar extinction coefficient of the reactant (in M-1 cm-1).
    :param epsilon_product: Molar extinction coefficient of the product (in M-1 cm-1).
    :param concentration_reactant: Concentration of the reactant (in M).
    :param concentration_product: Concentration of the product (in M).
    :return: Fraction of photons absorbed by the reactant.
    """
    return (epsilon_reactant * concentration_reactant) / ((epsilon_reactant * concentration_reactant) + (epsilon_product * concentration_product))

def apply_gaussian_broadening(excitation_energy, transition_strength, solar_wavelengths, broadening_factor):
    """
    Applies Gaussian broadening to the solar spectrum for a specific molecular transition.
    :param excitation_energy: Energy at which the molecule is excited (in eV).
    :param transition_strength: Strength of the molecular transition.
    :param solar_wavelengths: Array of solar wavelengths (in nm).
    :param broadening_factor: Standard deviation of the Gaussian broadening (in eV).
    :return: Tuple of energy values (in eV) and broadened transition intensities.
    """
    energy_values = planck_constant * speed_of_light / (solar_wavelengths * 1e-9 * electronvolt_to_joule)
    broadened_intensities = np.zeros_like(solar_wavelengths)
    normalization_factor = 1.3062974 * 10 ** 8.0
    for i, energy in enumerate(energy_values):
        gaussian_term = np.exp(-((energy - excitation_energy) ** 2) / (broadening_factor ** 2))
        normalization = normalization_factor * (transition_strength / (8065.73 * broadening_factor))
        broadened_intensities[i] += normalization * gaussian_term

    return energy_values, broadened_intensities

def calculate_SCE(storage_energy_Hartree, reaction_barrier_kJ, excitation_wavelength_nm, product_excitation_wavelength_nm, transition_strength, product_transition_strength, irradiance_spectrum=solar_irradiance, wavelength_spectrum=wavelengths):
    """
    Calculates the Solar Conversion Efficiency for a molecule based on its energy storage and reaction properties.
    :param storage_energy_Hartree: Energy stored in the molecule (in Hartree).
    :param reaction_barrier_kJ: Back reaction energy barrier (in kJ/mol).
    :param excitation_wavelength_nm: Wavelength of maximum excitation (in nm).
    :param product_excitation_wavelength_nm: Wavelength of maximum excitation of product (in nm).
    :param transition_strength: Strength of the excitation transition.
    :param product_transition_strength: Strength of the excitation transition of product.
    :param irradiance_spectrum: Array of solar irradiance values (in W/m2).
    :param wavelength_spectrum: Array of solar wavelengths (in nm).
    :param total_solar_power: Total solar power hitting the Earth's surface (in W).
    :return: Solar Conversion Efficiency (SCE) as a percentage.
    """
    # chech if reaction barrier is given in Hartree or kJ/mol
    if reaction_barrier_kJ < 5:
        reaction_barrier_kJ = reaction_barrier_kJ * 2625.5

    # Convert storage energy from Hartree to kJ/mol
    storage_energy_kJ = storage_energy_Hartree * 2625.5
    
    # Calculate the energy cutoff for the solar spectrum
    energy_cutoff = (planck_constant * speed_of_light) / (1e-6 * (storage_energy_kJ + reaction_barrier_kJ) / avogadro_number)

    total_solar_power = calculate_total_solar_power(irradiance_spectrum, wavelength_spectrum)

    # Apply Gaussian broadening to simulate the absorption profile
    excitation_energy_eV = 1239.8 / excitation_wavelength_nm  # Convert wavelength to energy
    energy_values, broadened_intensities = apply_gaussian_broadening(excitation_energy_eV, transition_strength, wavelength_spectrum, 0.25)
    # Apply Gaussian broadening to simulate the absorption of product
    product_energy_eV = 1239.8 / product_excitation_wavelength_nm  # Convert wavelength to energy
    product_energy_values, product_broadened_intensities = apply_gaussian_broadening(product_energy_eV, product_transition_strength, wavelength_spectrum, 0.25)
    
    # Calculate the attenuation factor
    attenuation_factor = 1.0 - 10 ** (-(broadened_intensities + product_broadened_intensities))

    # CAUTION: This is a very crude approximation of the fraction of photons absorbed by the reactant
    #          It is only valid for a single transition and should be replaced by a more general expression
    #          for the fraction of photons absorbed by the reactant
    photon_absorption_fraction = calculate_photon_absorption_fraction(broadened_intensities, product_broadened_intensities)  

    # Calculate absorbed photon count
    absorbed_photon_count = 0.0
    for i, wavelength in enumerate(wavelength_spectrum):
        if wavelength >= energy_cutoff:
            break
        photon_energy = planck_constant * speed_of_light / (wavelength * 1e-9)
        if i == 0:
            delta_wavelength = 0.5
        else:
            delta_wavelength = wavelength - wavelength_spectrum[i - 1]
        absorbed_photon_count += (irradiance_spectrum[i] / photon_energy) * delta_wavelength * attenuation_factor[i] * photon_absorption_fraction[i]

    # Calculate solar conversion efficiency
    sce_percentage = max(0.0, (storage_energy_kJ * 1e3 / avogadro_number) * absorbed_photon_count / total_solar_power) * 100

    return sce_percentage

