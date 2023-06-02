#!/groups/kemi/koerstz/anaconda3/envs/quantum_ml/bin/python

import pickle
import pandas as pd
import numpy as np
import sys
sys.path.append("/groups/kemi/brq616/speciale/opt/xTB/tQMC/QMC")
import qmconf


file_names = np.loadtxt("batch_list_final.csv", dtype=str) #a csv file with all the batch names without .pkl
file_names = file_names.tolist()

print(file_names)
frames = []

for i, file_name in enumerate(file_names):
    name_file = file_name + ".pkl"
    frame = pd.read_pickle(name_file) 
    frames.append(frame)
    if i%100 == 0:
        print(f"reac {i}")

frames_prod = []

for i, file_name in enumerate(file_names):
    name_file = file_name + "_prod.pkl"
    frame = pd.read_pickle(name_file)
    frames_prod.append(frame)
    if i%100 == 0:
        print(f"prod {i}")



print("reac concat starting")
result = pd.concat(frames)
print("prod concat starting")
result_prod = pd.concat(frames_prod)
print("concat done")


result["abs_prod"] = result_prod["max_abs"]
result["osc_prod"] = result_prod["osc_str"]

result.to_pickle('final_result.pkl')
result.to_csv("final_result.csv")
#result_prod.to_csv("final_result_prod.csv")
