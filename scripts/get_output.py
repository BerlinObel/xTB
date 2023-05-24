#!/groups/kemi/koerstz/anaconda3/envs/quantum_ml/bin/python

import pickle
import pandas as pd
import numpy as np
import sys
from rdkit import Chem
from rdkit.Chem import Descriptors
sys.path.append("/groups/kemi/obel/opt/tQMC/QMC")
import qmconf
pd.set_option('display.max_rows', None)
pd.set_option('display.max_columns', None)

df = pd.read_pickle(sys.argv[1])
print(df)
