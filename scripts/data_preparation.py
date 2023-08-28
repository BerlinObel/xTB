import sys
import os
import pandas as pd
from glob import glob
from tqdm import tqdm


# Add the directory of the script to the Python path
sys.path.append(os.path.dirname(os.path.realpath(__file__)))

from settings import QMC_PATH
sys.path.append(QMC_PATH)

from QMC.qmmol import QMMol
from QMC.qmconf import QMConf

def create_batches(pattern, file_suffix):
    file_paths = pattern
    files = glob(file_paths)
    batch_names = [file_name.split(".")[0] for file_name in files]
    batch_name_df = pd.DataFrame(batch_names)
    batch_name_df.to_csv(f"batch_list{file_suffix}.csv",index=False,header=False)

def gather_results(pattern):
    file_paths = pattern
    files = glob(file_paths)
    all_results = []
    batch_names = [file_name.split(".")[0] for file_name in files]
    
    for batch_name in tqdm(batch_names):
        file_name = batch_name + ".pkl"
        results = pd.read_pickle(file_name) 
        all_results.append(results)
        
    all_results_df = pd.concat(all_results)
    all_results_df.to_pickle('final_result.pkl')
    all_results_df.to_csv("final_result.csv")

if __name__ == "__main__":
    job_type = sys.argv[1]

    if job_type == "tbr":
        create_batches("smiles_batch*.pkl", "")
    elif job_type == "abs":
        create_batches("smiles_batch*_new.pkl", "_abs")
    elif job_type == "final":
        gather_results("smiles_batch*_abs.pkl")
    else:
        print(f"Unknown job type: {job_type}")


    

