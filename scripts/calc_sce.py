
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import math
from rdkit import Chem



try:
    am15_sheet = pd.read_pickle('am15_SMARTS2.pkl')
except:
    am15_path = "/lustre/hpc/kemi/elholm/bod_calc/alle/eff/am_15g.xls"
    am15_sheet = pd.read_excel(am15_path,sheet_name = "SMARTS2",header=None)
    pd.to_pickle(am15_sheet,'am15_SMARTS2.pkl')


# Read solar spectrum - choosen spectra normalized to 1000 (1000 W/m2)
solspec = am15_sheet.to_numpy()
wavelength = solspec[2:,0]
Esun = solspec[2:,2]

# Calculate Edot
def calc_Edot():
    return Esun[0] * 0.5 + sum(E * (w2 - w1) for E, w1, w2 in zip(Esun[1:], wavelength[:-1], wavelength[1:]))

Edot = calc_Edot() #should be 999.7102563234429

# Natural constants
hbar = 1.05457182*10**(-34)
h = 6.6261*10**(-34)
c = 299792458
Na = 6.0221408*10**(23)
evtoJ = 1.602176634*10**(-19)


def gaussian(x, y, wavelength, gamma):
    '''
    Gaussian broadening function
     
    Call: xi,yi = gaussian(energies, intensities, start energy, end stop energy, energy step, gamma)
    '''
 
    factor = 1.3062974*10**8.0
    xi = h*c/(wavelength*10**(-9)*evtoJ)
    yi=np.zeros_like(wavelength)
    for i in range(len(xi)):
        for k in range(len(x)): yi[i] = yi[i] + factor * (y[k]/(8065.73*gamma)) * np.e**(-((xi[i]-x[k])**2)/(gamma**2))
     
    return xi,yi

def SCE(storage,barrier,exci,osc,Esun=Esun,wavelength=wavelength):
    print(f'###----------------###\nStorage Energy: {storage} (Hartree)\n{storage*2625.5} (kJ/mol)\n{(storage*2625.5)*10**(3)/Na}(eV ?)\nTBR Barrier: {barrier}\nAbsorptions wavelength: {exci}\nOscilator strength: {osc}\nSpectra:\n{wavelength}\n{Esun}\n')
    
    storage = storage*2625.5 #Convert Hartree to kJ/mol
    Ndot = 0.0
    sce = 0.0
    Ecut = (h*c)/(10**(-6)*(storage+barrier)/Na)


    xi,yi = gaussian(np.array([1239.8/exci]),np.array([osc]),wavelength,0.25)

    attenuation = (1.0 - 10**(-yi))

    # Calculate Ndot times attenuation
    for i in range(len(wavelength[wavelength<Ecut])):
        if i==0:
            Ndot += (Esun[i]/(h*c/(wavelength[i]*10**(-9))))* 0.5 * attenuation[i]
        else:
            Ndot += (Esun[i]/(h*c/(wavelength[i]*10**(-9))))* (wavelength[i]-wavelength[i-1]) * attenuation[i]

    sce = ((storage*10**(3)/Na)) * Ndot/Edot


    if sce < 0.0:
        sce = 0.0

    return sce*100
