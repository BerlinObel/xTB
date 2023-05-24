#!/groups/kemi/koerstz/anaconda3/envs/quantum_ml/bin/python

import pickle
import pandas as pd
import numpy as np
import sys
import os
import xtb
from rdkit import Chem
from rdkit.Chem import AllChem
sys.path.append("/groups/kemi/obel/opt/tQMC/QMC")
from qmmol import QMMol
from qmconf import QMConf
import myio.io_ts_scan_xtb as ts_scan_xtb
from calculator.xtb import xTB

tstype=['a','b','c']

def create_xcontrol(mol,xcontrol_name,scan):
    
    #Find the pairs
    smart_mol = Chem.MolFromSmarts('[#6]/[#7]=[#7]\\[#6]')
    mol = Chem.AddHs(mol)
    
    sub = np.array(mol.GetSubstructMatches(smart_mol))
    mol.Get
    print(sub)
    
    #Write the xcontrol file
    fo = open(xcontrol_name, "w")
    if scan == 'a':
        fo.write("$constrain\n")
        fo.write(" force constant=0.5\n")
        fo.write(" dihedral: "+str(sub[0][0])+", "+str(sub[0][1])+", "+str(sub[0][2])+", "+str(sub[0][3])+", auto\n")
        fo.write("$scan\n")
        fo.write(" 1: 0, 180, 100\n")
        fo.write("$end\n")
    elif scan == 'b':
        fo.write("$constrain\n")
        fo.write(" force constant=0.5\n")
        fo.write(" angle: "+str(sub[0][0])+", "+str(sub[0][1])+", "+str(sub[0][2])+", auto\n")
        fo.write("$scan\n")
        fo.write(" 1: 120, 240, 100\n")
        fo.write("$end\n")
    elif scan == 'c':
        fo.write("$constrain\n")
        fo.write(" force constant=0.5\n")
        fo.write(" angle: "+str(sub[0][1])+", "+str(sub[0][2])+", "+str(sub[0][3])+", auto\n")
        fo.write("$scan\n")
        fo.write(" 1: 120, 240, 100\n")
        fo.write("$end\n")
    fo.close()
    return

def ts_scan(prod, name, chrg, mult, xcontrol_name):
    """ scan for transition state """
    
    charged = True # hard coded for mogens
    
    # create conformers
    ts_qmmol = QMMol()
    ts_qmmol.add_conformer(prod.write_xyz(to_file=False), fmt='xyz', label=name, charged_fragments=charged, set_initial=True)
    
    xtb_params = {'method': 'gfn1',
                  'opt': 'opt',
                  'cpus': 1,
                  'input': '../' + str(xcontrol_name)}

    ts_qmmol.calc = xTB(parameters=xtb_params)
    ts_conf = ts_qmmol.conformers[0]

    #ts_conf.conf_calculate(quantities=['energy', 'structure'], keep_files=True)
    #ts_conf.conf_calculate(quantities=['energy'], keep_files=True)
    ts_conf.conf_calculate(quantities=['energy', 'ts_guess'], keep_files=True)

def ts_path(prod, reac):
    """ scan for transition state """
    
    charged = True # hard coded for mogens
    
    # create conformers
    ts_qmmol = QMMol()
    ts_qmmol.add_conformer(prod.write_xyz(to_file=False), fmt='xyz', label=name, charged_fragments=charged, set_initial=True)
    
    xtb_params = {'method': 'gfn2',
                  'path': 'path.inp',
                  'cpus': 1,
                  'input': '../' + str(xcontrol_name)}

    ts_qmmol.calc = xTB(parameters=xtb_params)
    ts_conf = ts_qmmol.conformers[0]

    #ts_conf.conf_calculate(quantities=['energy', 'structure'], keep_files=True)
    #ts_conf.conf_calculate(quantities=['energy'], keep_files=True)
    ts_conf.conf_calculate(quantities=['energy', 'ts_guess'], keep_files=True)
    
    return



if __name__ == '__main__':

    df = pd.read_pickle(sys.argv[1])
    charge = 0
    mult = 1
    compound_list = list()
    

    for x in df.itertuples():
        scan_list = list()
        ts_qmconf_a = QMConf()
        ts_qmconf_b = QMConf()
        ts_qmconf_c = QMConf()
        tbr_a = 0
        tbr_b = 0
        tbr_c = 0
        for q in tstype:
    
            name = str(x.reac.label.split('_r')[0])
            mol = x.prod.get_rdkit_mol()
            print(name,q)
            xcontrol_name = name+"_xcontrol"
            create_xcontrol(mol, xcontrol_name, q)

            n = name + "_ts-" + q
            
            ts_scan(x.prod, n, charge, mult, xcontrol_name)
            
            with open(n+"_xtbscan.log", 'r') as out:
                output = out.read() 
            
            print(output)
            ts_qmmol = QMMol()
            
            xyz_file = ts_scan_xtb.read_ts_guess_structure2(output)
            ts_qmmol.add_conformer(input_mol=xyz_file, fmt='xyz', label=n, charge=charge, multiplicity=mult, read_freq=False, charged_fragments=True)
            ts_conf = ts_qmmol.conformers[0]

            ts_conf.results['energy'] = ts_scan_xtb.read_ts_guess_energy(output)[0]
            
            #print(ts_conf.results['energy'])
            reac_qmconf = x.reac
            prod_qmconf = x.prod
            ts_qmconf = ts_conf
            storage = x.storage*2625.5 #convert from Hartree to kJ/mol
            tbr = (ts_conf.results['energy']-x.prod.results['energy'])*2625.5 #convert from Hartree to kJ/mol
            
            print(tbr, ts_conf.results['energy'], x.prod.results['energy'])
            
            if q == 'a':
                ts_qmconf_a = ts_qmconf
                tbr_a = tbr
            elif q == 'b':
                ts_qmconf_b = ts_qmconf
                tbr_b = tbr
            elif q == 'c':
                ts_qmconf_c = ts_qmconf
                tbr_c = tbr
            #Clean up
            os.remove(str(n) + ".out")
            os.remove(str(n) + "_xtbscan.log")
            os.remove(str(n) + ".xyz")

        

        compound_list.append({'rep': x.rep,
                                'reac': reac_qmconf,
                                'prod': prod_qmconf,
                                'scan_TSa': ts_qmconf_a,
                                'scan_TSb': ts_qmconf_b,
                                'scan_TSc': ts_qmconf_c,
                                'storage': storage,
                                'tbr_a': tbr_a,
                                'tbr_b': tbr_b,
                                'tbr_c': tbr_c})


    structures = compound_list
    data = pd.DataFrame(structures)
    data.to_pickle(sys.argv[1].split('/')[0].split('.')[0] + '_new.pkl')
    #data.to_pickle(sys.argv[1].split('.')[0] + '_new.pkl') #IF running on frontend this should be used!
