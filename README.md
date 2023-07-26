# Automated xTB Scripts

## Acknowledgement (not finsihed)
This work is based on the initial scripts and research conducted by [Initial Author/Contributors Names include links to theirs githubs and all of that]. 

Their fundamental contribution to this project is gratefully acknowledged.

## Description
This repository contains scripts for automating the calculation of storage energies, theremal back reaction barriers, and maximum absorption wavelengths using the xTB method. The results are saved as pickled pandas dataframes for efficient storage and retrieval. Each molecule is stored as an RDKit object within a QMmol object, allowing easy access to SMILES, XYZ coordinates, and more.

## File descriptions
`data_preparation.py`

`electronic_azo_conf_search.py`

`find_max_abs.py`

`find_ts.py`

`get_xyz.py`

`read_pickle.py`

`settings.py`

`submit_calculations.py`

`utils.py`

## Output File Structure
Output files are generated at each stage of the workflow and saved in pickled format. Each output file is a DataFrame where each row corresponds to a molecule. The columns contain the properties calculated in each step, and the molecule itself is stored as an RDKit Mol object within a QMmol object.

## Prerequisites (not finished)
Ensure the following dependencies are installed:
- Python 3.x
- Pandas
- Numpy
- Scipy
- RDKit
- and more to be determined

## Workflow
1. **Prerequisites (not finished):**
   
   Ensure that the dependencies are installed
   Also to be determined how to reliably set up the xTB and stda program to run for anyone
   Change paths and username in the settings.py file

2. **Input file:**

   The input file should be a CSV file with the following structure:

    - Column A, named `comp_name`, should include compound names (e.g., `nbd_a1_b2_c3_d4`).
    - Column B, named `smiles`, should include the corresponding SMILES string.
    - Column C, named `charge`, should include the corresponding charge (often 0).
    - Column D, named `multiplicity`, should include the corresponding multiplicity (often 1).
    - Column E, named `molecular_weight`, should include molecular weight based on smiles.

3. **Calculate Storage Energies:**

    Run: `python submit_calculations.py $input_file.csv storage`

    You'll be prompted to enter the batch size. The script will first create the product molecule based on the reactant SMILES input then calculate and store the storage energies.

4. **Find Transition States and Thermal Back Reaction Barriers:**

    Run the data preparation script: `python data_preparation.py tbr`

    Then, run: `python submit_calculations.py batch_list.csv tbr`

    The script will calculate and store the thermal back reaction barriers.

5. **Calculate Maximum Absorption Wavelength and Oscillator Strength:**

    Prepare the data by running: `python data_preparation.py abs`

    Then, run: `python submit_calculations.py $batch_list_abs.csv abs`

    The script will find the excitation energies then calculate and store the maximum absorption wavelengths and oscillator strengths.

7. **Merge Data:**

    To gather all the data into one file, run: `python data_preparation.py final`

    The final results will be saved to `final_results.csv`.

7. **Extract .xyz files:**

    To extract all xyz structures, run: `python get_xyz.py $final_results.pkl all`

    To extract only transistion state xyz structures, run: `python get_xyz.py $final_results.pkl ts`

