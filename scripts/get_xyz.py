#!/groups/kemi/koerstz/anaconda3/envs/quantum_ml/bin/python

import pickle
import pandas as pd
import numpy as np
import sys
from rdkit import Chem
from rdkit.Chem import Descriptors
sys.path.append("/groups/kemi/brq616/speciale/opt/xTB/QMC/QMC")
import qmconf
import os
import time
'''
numbers = ['a6_b0_c5_d2','a6_b0_c2_d1','a6_b0_c5_d4','a6_b0_c8_d4','a6_b0_c5_d5','a6_b0_c8_d1','a6_b0_c8_d2','a6_b0_c6_d2','a6_b0_c8_d0',
'a6_b0_c5_d0','a6_b0_c4_d2','a6_b0_c0_d2','a6_b0_c6_d1','a6_b0_c8_d3','a3_b0_c8_d5','a6_b0_c2_d5','a6_b0_c1_d5',
'a6_b0_c7_d2','a6_b0_c2_d2','a6_b0_c5_d6','a6_b0_c1_d2','a4_b0_c8_d5','a6_b0_c8_d6','a6_b0_c5_d3','a6_b0_c0_d5'
]
'''
#comp_names = pd.read_csv("comp_names.csv")
#for name in comp_names["nr"]:
#    print(name)
#names_list = [x for x in comp_names["nr"]]

df = pd.read_pickle(sys.argv[1])
if sys.argv[2] == 'all':
	print(df.head())
	for x in df.itertuples():
		#if x.rep in names_list:
		print(x)
		x.prod.write_xyz(to_file=True)
		x.reac.write_xyz(to_file=True)
		print(x.reac.label)
		x.scan_TSa.write_xyz(to_file=True)
		x.scan_TSb.write_xyz(to_file=True)
		x.scan_TSc.write_xyz(to_file=True)
		#break
elif sys.argv[2] == 'ts':
	print(df.head())
	for x in df.itertuples():
		#if x.rep in names_list:
		print(x)
		# x.prod.write_xyz(to_file=True)
		# x.reac.write_xyz(to_file=True)
		# print(x.reac.label)
		x.scan_TSa.write_xyz(to_file=True)
		x.scan_TSb.write_xyz(to_file=True)
		x.scan_TSc.write_xyz(to_file=True)
		#break
else:
	print(df.head())
	for x in df.itertuples():
		#if x.rep in names_list:
		print(x)
		x.prod.write_xyz(to_file=True)
		x.reac.write_xyz(to_file=True)
		# print(x.reac.label)
		# x.scan_TSa.write_xyz(to_file=True)
		# x.scan_TSb.write_xyz(to_file=True)
		# x.scan_TSc.write_xyz(to_file=True)
		#break
#for x in range(10):
#	number = input("What is the number of your structure?")
#	numbers.append(str(number))
'''
for x in df.itertuples():
	#print(x.rep)
	if x.rep in numbers:
		x.reac.write_xyz(to_file=True)
		print("It was there")
	else:
		continue
'''
#com_path = "/groups/kemi/elholm/bin/g16_elholm"
#dft_path = "/groups/kemi/elholm/bin/elholm_g16"
#for x in numbers:
	#os.mkdir(f"/groups/kemi/elholm/tddft/{x}")
#	os.popen(f"cp /groups/kemi/elholm/xyz_files/{x}_r*.xyz /groups/kemi/elholm/tddft/.")
	#os.popen(f"{com_path} /groups/kemi/elholm/tddft/{x}/{x}*.xyz")
	#time.sleep(2)
	#os.popen(f"cp /groups/kemi/elholm/bin/elholm_g16 /groups/kemi/elholm/tddft/{x}/.")
	#time.sleep(2)
	#os.popen(f"/groups/kemi/elholm/tddft/{x}/elholm_g16 /groups/kemi/elholm/tddft/{x}/{x}*.com")

print("All done!")

