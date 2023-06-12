#!/groups/kemi/koerstz/anaconda3/envs/quantum_ml/bin/python

import pickle
import pandas as pd
import numpy as np
import sys
from rdkit import Chem
from rdkit.Chem import Descriptors
sys.path.append("/groups/kemi/brq616/speciale/opt/xTB/tQMC/QMC")
import qmconf
import os


df = pd.read_pickle(sys.argv[1])
if sys.argv[2] == 'all':
    print(df.head())
    for x in df.itertuples():
        print(x)
        x.prod.write_xyz(to_file=True)
        x.reac.write_xyz(to_file=True)
        print(x.reac.label)
        x.ts.write_xyz(to_file=True)

elif sys.argv[2] == 'ts':
    print(df.head())
    for x in df.itertuples():
        print(x)
        x.ts.write_xyz(to_file=True)
else:
    print(df.head())
    for x in df.itertuples():
        print(x)
        x.prod.write_xyz(to_file=True)
        x.reac.write_xyz(to_file=True)

print("All done!")

