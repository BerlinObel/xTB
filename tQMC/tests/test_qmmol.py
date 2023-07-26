import sys
import os

import numpy as np

sys.path.append('/groups/kemi/ree/opt/tQMC/QMC')
from qmmol import QMMol




def create_qmmol():
    """ Create qmmol, if qmconf works this is ok """
    test_dir = os.path.dirname(os.path.realpath(__file__)) + '/data/formats/'

    qmmol = QMMol()
    qmmol.add_conformer(input_mol=test_dir+'propanol.xyz', fmt='xyz',
                        label='propanol', charge=0, multiplicity=1,
                        charged_fragments=True, set_initial=True)

    return qmmol


def test_systematic_conf_generation():
    """ """
    qmmol = create_qmmol()
    qmmol.create_systematic_conformers(theta=120, threads=1, ff_variant='uff',
                                       max_iters=200, check_stero=True)
    
    return qmmol

def test_random_conf_generations():
    """ """
    qmmol = create_qmmol()
    qmmol.create_random_conformers(num_confs=20, threads=1)

    return qmmol


def test_parallel_execution_of_qm_calc():
    """ """
    from calculator.gaussian import Gaussian

    qmmol = test_systematic_conf_generation()
    
    qmmol.set_calculator(Gaussian())
    qmmol.calculate(num_procs=4, keep_files=True)

    for conf in qmmol.conformers:
        print(conf.results)
        

def test_clustering():
    """ """
    pass


if __name__ == '__main__':
    test_systematic_conf_generation()
    test_parallel_execution_of_qm_calc()


