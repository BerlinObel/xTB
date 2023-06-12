import pandas as pd
from glob import glob
import sys
from tqdm import tqdm

sys.path.append("/groups/kemi/brq616/speciale/opt/xTB/tQMC/QMC")

from qmmol import QMMol
from qmconf import QMConf

if __name__ == "__main__":
    

    if sys.argv[1] == "tbr":
        file_paths = "smiles_batch*.pkl"
        files = glob(file_paths)
        batch_list_name = "batch_list.csv"
        batch_names = [file_name.split(".")[0] for file_name in files]
        batch_name_df = pd.DataFrame(batch_names)
        batch_name_df.to_csv(batch_list_name,index=False,header=False)

    if sys.argv[1] == "abs": 
        file_paths = "smiles_batch*_new.pkl"
        files = glob(file_paths)
        batch_list_name = "batch_list_abs.csv"
        batch_names = [file_name.split(".")[0] for file_name in files]
        batch_name_df = pd.DataFrame(batch_names)
        batch_name_df.to_csv(batch_list_name,index=False,header=False)       

    if sys.argv[1] == "final":
        file_paths = "smiles_batch*_abs.pkl"
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


    

