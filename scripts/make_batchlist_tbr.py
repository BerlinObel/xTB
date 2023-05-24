import pandas as pd
from glob import glob

files = glob("smiles_batch*.pkl")
names = [x.split(".")[0] for x in files]
df = pd.DataFrame(names)
df.to_csv("batch_list.csv",index=False,header=False)
