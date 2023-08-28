#!/groups/kemi/koerstz/anaconda3/envs/quantum_ml/bin/python

import pickle
import pandas as pd
import numpy as np
import sys
import os
from rdkit import Chem
from rdkit.Chem import Descriptors
# Add the directory of the script to the Python path
sys.path.append(os.path.dirname(os.path.realpath(__file__)))

from settings import QMC_PATH
sys.path.append(QMC_PATH)

from QMC.qmconf import QMConf


pd.set_option('display.max_rows', None)
pd.set_option('display.max_columns', None)

df = pd.read_pickle(sys.argv[1])
print(df)
