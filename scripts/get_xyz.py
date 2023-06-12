import pandas as pd
import sys
from rdkit import Chem
from rdkit.Chem import Descriptors
sys.path.append("/groups/kemi/brq616/speciale/opt/xTB/tQMC/QMC")
import qmconf
import os

def write_xyz(compound, name):
    filename = f"{name}_{compound.label}.xyz"
    compound.write_xyz(to_file=True)
    print(f"Written file: {filename}")

def process_compound(compound, mode):
    print(compound)

    if mode in ['all', 'ts']:
        write_xyz(compound.ts, "ts")

    if mode in ['all']:
        write_xyz(compound.prod, "prod")
        write_xyz(compound.reac, "reac")

def main():
    df = pd.read_pickle(sys.argv[1])
    print(df.head())

    mode = sys.argv[2]

    if mode not in ['all', 'ts']:
        raise ValueError("Invalid mode. Use 'all' or 'ts'.")

    for compound in df.itertuples():
        process_compound(compound, mode)

    print("All done!")

if __name__ == "__main__":
    main()
