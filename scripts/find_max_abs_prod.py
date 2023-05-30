#!/groups/kemi/koerstz/anaconda3/envs/quantum_ml/bin/python

import pickle
import pandas as pd
import numpy as np
import subprocess
import sys
import os

from rdkit import Chem
from rdkit.Chem import Descriptors
from rdkit.Chem import AllChem
sys.path.append('/groups/kemi/elholm/opt/xTB_sTDA')
import vis_osc_pkl
from vis_osc_pkl import get_max_abs
sys.path.append("./tQMC/QMC")
import qmconf
from qmmol import QMMol
from qmconf import QMConf


def shell(cmd, shell=False):

    if shell:
        p = subprocess.Popen(cmd, shell=True, stdin=subprocess.PIPE, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        p.wait()
    else:
        cmd = cmd.split()
        p = subprocess.Popen(cmd, stdin=subprocess.PIPE, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        p.wait()

    output, err = p.communicate()
    return output



def compute_absorbance_2(mol, charge, spin, n):
    print(f"{mol}, {charge}, {spin}, {n}")
    #shell("/groups/kemi/koerstz/opt/xtb/5.8/xtb "+str(mol)+" -chrg "+str(charge)+" -uhf "+str(spin),shell=False)
    print("/groups/kemi/ree/opt/stda/xtb4stda "+str(mol)+" -chrg "+str(charge)+" -uhf "+str(spin))
    shell("/groups/kemi/ree/opt/stda/xtb4stda "+str(mol)+" -chrg "+str(charge)+" -uhf "+str(spin),shell=False)
    #out = shell('/groups/kemi/ree/opt/xTB_sTDA/stda_1.6 -xtb -e 10',shell=False)
    out = shell('/groups/kemi/ree/opt/stda/stda_v1.6.1 -xtb -e 10',shell=False)
    print(out)
    data = str(out).split('Rv(corr)\\n')[1].split('(')[0]
    wavelength, osc_strength = float(data.split()[2]), float(data.split()[3])

    data2 = str(out).split('Rv(corr)\\n')[1].split('\\n')
    if len(data2) > 100:
        for i in range(0,n): #number of excited states
            print(float(data2[i].split()[2]), float(data2[i].split()[3]), file=open("plot_info.txt", "a"))
    else:
        print(mol, "do not have", n, "excitations, but", int(len(data2)-25), "excitations")

        for i in range(0,len(data2)-25):
            print(float(data2[i].split()[2]), float(data2[i].split()[3]), file=open("plot_info.txt", "a"))
      
    return wavelength, osc_strength


def run_abs_calc(df_x):
    #Define threshold here
    threshold = 500

    n = 75 #select the number of excited states
    charge = '0' #charge of the system for the xTB calculation
    spin = '0' #spin of the system for the xTB calculation

  
    xyz = str(df_x.prod.label)+'.xyz'
    df_x.prod.write_xyz(to_file=True)
    wavelength, osc_strength = compute_absorbance_2(xyz, charge, spin, n)
    shell("sed -i '1i \ UV' tda.dat",shell=True)
    shell("/groups/kemi/elholm/opt/xTB_sTDA/g_spec < tda.dat",shell=True)
    
    x1,x2 = get_max_abs('plot_info.txt')
    for x in x2:
      if x < threshold:
        ind = np.where(x2==x)
        np.delete(x1,ind)
        np.delete(x2,ind)
    max_ind = np.where(x1==x1.max())
  
    max_abs = float(x1[max_ind])
    osc_str = float(x2[max_ind])
  
    os.remove('plot_info.txt')
    os.remove('charges')
    os.remove('tda.dat')
    os.remove('rots.dat')
    os.remove('spec.dat')
    os.remove('wbo')
    os.remove('wfn.xtb')
    os.remove(xyz)

    return max_abs, osc_str



def run_abs_calc_first_exc(df_x):
    #Define threshold here
    threshold = 700

    n = 75 #select the number of excited states
    charge = '0' #charge of the system for the xTB calculation
    spin = '0' #spin of the system for the xTB calculation

    xyz = str(df_x.prod.label)+'.xyz'
    df_x.prod.write_xyz(to_file=True)
    wavelength, osc_strength = compute_absorbance_2(xyz, charge, spin, n)
    shell("sed -i '1i \ UV' tda.dat",shell=True)
    shell("/groups/kemi/elholm/opt/xTB_sTDA/g_spec < tda.dat",shell=True)
    
    data1 = np.loadtxt('plot_info.txt')
    l1,f1 = data1[:,0],data1[:,1]
    l1,f1 = l1[::-1],f1[::-1]
    for i in range(len(f1)):
      if f1[i]>= 0.02 and l1[i] <= 700:
        first_abs = l1[i]
        first_osc = f1[i]
    
    first_abs = float(first_abs)
    osc_str = float(first_osc)
  
    os.remove('plot_info.txt')
    os.remove('charges')
    os.remove('tda.dat')
    os.remove('rots.dat')
    os.remove('spec.dat')
    os.remove('wbo')
    os.remove('wfn.xtb')
    os.remove(xyz)

    return first_abs, osc_str



if __name__ == '__main__':

    df = pd.read_pickle(sys.argv[1])
    compound_list = list()

    for x in df.itertuples():
        reac_qmconf = x.reac
        prod_qmconf = x.prod
        ts_qmconf = x.scan_ts
        storage = x.storage
        tbr = x.tbr
        max_abs, osc_str = run_abs_calc_first_exc(x)

        compound_list.append({'rep': x.rep,
                              'reac': reac_qmconf,
                              'prod': prod_qmconf,
                              'scan_ts': ts_qmconf,
                              'storage': storage,
                              'tbr': tbr,
                              'max_abs': max_abs,
                              'osc_str': osc_str})
        

    data = pd.DataFrame(compound_list)
    data.to_pickle(sys.argv[1].split('/')[0].split('.')[0] + '_abs_prod.pkl')
    #data.to_pickle(sys.argv[1].split('.')[0] + '_abs.pkl') #IF running on frontend this should be used!
