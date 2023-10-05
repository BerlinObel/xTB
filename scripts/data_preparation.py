import sys
import os
import pandas as pd
from glob import glob
from tqdm import tqdm


# Add the directory of the script to the Python path
sys.path.append(os.path.dirname(os.path.realpath(__file__)))

from db_utils import create_table, insert_data, set_calculation_stage
from settings import QMC_PATH
sys.path.append(QMC_PATH)

from qmmol import QMMol
from qmconf import QMConf

def create_batches(input_file):
    chunk_size = int(input("Batch size: "))
    data = pd.read_csv(input_file)
    chunks = [data[i:i+chunk_size] for i in range(0, data.shape[0], chunk_size)]

    for idx, chunk in enumerate(chunks):
        chunk_name = f"batch_{idx}"
        for _, row in chunk.iterrows():
            molecule_data = {
                'HashedName': row['comp_name'],
                'Smiles': row['smiles'],
                'Multiplicity': row['multiplicity'],
                'Charge': row['charge'],
                'BatchID': chunk_name,
                'CalculationStage': 'storage'
                }
            insert_data(molecule_data)

if __name__ == "__main__":
    job_type = sys.argv[1]
    
    if job_type == "fill":
        input_file = sys.argv[2]
    else:
        input_file = None
        
    if job_type == "create":
        create_table()
    elif job_type == "fill":
        if input_file is None:
            print("Data file must be provided for storage job type")
            sys.exit(1)
        create_batches(input_file)
        
    elif job_type == "change":
        wanted_stage = input("CalculationStage: ")
        set_calculation_stage(wanted_stage)
    else:
        print(f"Unknown job type: {job_type}")


    

