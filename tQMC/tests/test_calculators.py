import sys
import os

sys.path.append('/groups/kemi/ree/opt/tQMC/QMC')
from qmconf import QMConf

import numpy as np

def get_qmconf():
    """Return QMConf object."""
    
    test_dir = os.path.dirname(os.path.realpath(__file__))
    conf = QMConf(test_dir + "/data/methane.xyz", fmt='xyz',
                  label='test')
    return conf 


def compare(a, b):
    
    a = np.asarray(a)
    b = np.asarray(b)
    
    for pair in zip(a,b):
        if np.all(abs(pair[0] - pair[1]) >= 1e-6):
            return False
    
    return True

def test_gaussian():
    from calculator.gaussian import Gaussian
    
    # check geom optimization
    ref_opt_energy = -0.0207553882901 
    ref_structure = [[-0.000000,  0.000000,  0.000000],
                     [-0.190177, -0.869811,  0.623603],
                     [-0.923345,  0.557314, -0.135846],
                     [ 0.374014, -0.323072, -0.968172],
                     [ 0.739510,  0.635567,  0.480415]]

    param = {'method': 'pm3',
             'basis': '',
             'opt': 'opt',
             'nproc': '2',
             'mem': '4GB'}

    qmconf = get_qmconf()
    qmconf.set_calculator(Gaussian(parameters=param))
    qmconf.conf_calculate(quantities=['energy', 'structure'])
    
    assert compare(ref_structure, qmconf.structure), "Gaussian opt failed"
    assert compare([ref_opt_energy], [qmconf.results['energy']]), "Gaussian opt failed"    

    print("Gaussian: passed.")


def test_orca():
    from calculator.orca import ORCA
    
    # check geometry optimization
    ref_opt_energy = -6.63457504538, 
    ref_structure = [[-2.073347e+00, 5.140800e-02,-9.000000e-06],
                     [-9.863040e-01, 5.139900e-02,-2.000000e-06],
                     [-2.435663e+00, 3.932880e-01, 9.661570e-01],
                     [-2.435660e+00,-9.562700e-01,-1.869700e-01],
                     [-2.435666e+00, 7.171550e-01,-7.791860e-01]]

    param = {'method': 'pm3',
             'basis': '',
             'opt': 'opt'}

    qmconf = get_qmconf()
    qmconf.set_calculator(ORCA(parameters=param))
    qmconf.conf_calculate(quantities=['energy', 'structure'])
    
    assert compare(ref_structure, qmconf.structure), "ORCA opt failed"
    assert compare([ref_opt_energy], [qmconf.results['energy']]), "ORCA opt failed"

    print("ORCA: passed.")


def test_dftb():
    from calculator.dftb import dftb
    
    # check geom optimization
    ref_opt_energy = -3.2311624309 
    ref_structure = [[-2.07333610e+00, 5.14116060e-02, 4.81180000e-06],
                     [-9.92940722e-01, 5.13919430e-02,-9.62400000e-07],
                     [-2.43345247e+00, 3.91206982e-01, 9.60263464e-01],
                     [-2.43344624e+00,-9.50105951e-01,-1.85844872e-01],
                     [-2.43346447e+00, 7.13075420e-01,-7.74432441e-01]]

    qmconf = get_qmconf()
    qmconf.set_calculator(dftb())
    qmconf.conf_calculate(quantities=['energy', 'structure'])
    
    assert compare(ref_structure, qmconf.structure), "GAMESS DFTB opt failed"
    assert compare([ref_opt_energy], [qmconf.results['energy']]), "GAMESS DFTB  opt failed"
    
    print("GAMESS DFTB: passed.")


def test_xtb():
    from calculator.xtb import xTB
    
    # check gecom optimization
    ref_opt_energy = -4.1752199 
    ref_structure = [[-2.07332805e+00, 5.13961598e-02,-1.90470483e-06],
                     [-9.91160210e-01, 5.14007859e-02,-3.47724511e-06],
                     [-2.43405082e+00, 3.91769333e-01, 9.61825533e-01],
                     [-2.43404631e+00,-9.51759736e-01,-1.86141695e-01],
                     [-2.43405461e+00, 7.14173457e-01,-7.75688456e-01]]

    qmconf = get_qmconf()
    qmconf.set_calculator(xTB())
    qmconf.conf_calculate(quantities=['energy', 'structure'])
    
    assert compare(ref_structure, qmconf.structure), "xTB opt failed"
    assert compare([ref_opt_energy], [qmconf.results['energy']]), "xTB  opt failed"

    print("xTB: passed.")


if __name__ == '__main__':
    test_gaussian()
    #test_orca()
    #test_dftb()
    test_xtb()
    
