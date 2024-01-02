import sys
import os
import pandas as pd
from glob import glob
from tqdm import tqdm
import datetime

# Add the directory of the script to the Python path
sys.path.append(os.path.dirname(os.path.realpath(__file__)))

from db_utils import create_table, delete_table, insert_data, set_calculation_stage, delete_table, create_orca_table, rename_table, backup_table, add_batch
from settings import QMC_PATH
sys.path.append(QMC_PATH)

# from qmmol import QMMol
# from qmconf import QMConf
from QMC.qmmol import QMMol
from QMC.qmconf import QMConf

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

    elif job_type == "orca":
        create_orca_table()
        
    elif job_type == "fill":
        if input_file is None:
            print("Data file must be provided for storage job type")
            sys.exit(1)
        create_batches(input_file)
        
    elif job_type == "change":
        wanted_stage = input("CalculationStage: ")
        set_calculation_stage(wanted_stage)

    elif job_type == "delete":
        delete_table()
        
    elif job_type == "rename":
        rename_table()
        
    elif job_type == "backup":
        original_table_name = input("Which table: ")
        # Get the current date and time
        current_time = datetime.datetime.now()
        # Format the date and time into a string (YYYYMMDD_HHMMSS)
        formatted_time = current_time.strftime("%Y%m%d_%H%M%S")
        backup_table_name = f"{original_table_name}_{formatted_time}"
        backup_table(original_table_name, backup_table_name)

    elif job_type == "add_batch":
        if len(sys.argv) < 3:
            print("Data file and Batch ID must be provided for add_batch job type")
            sys.exit(1)
        input_file = sys.argv[2]
        batch_id = sys.argv[3]
        add_batch(input_file, batch_id)

    else:
        print(f"Unknown job type: {job_type}")


    

