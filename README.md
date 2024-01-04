# Automated xTB Scripts

Python scripts to automate calculations in computational chemistry using the xTB method.

## Acknowledgement (Not Finished)

This project is built on foundational scripts and research by [Add more names and repositories here]. 

- The **QMC package**, originally by Mads Koerstz, is used here with some changes. It helps in storing molecular shapes and calculation outputs. Check out the original at [QMC Repository](https://github.com/koerstz/QMC).
  
- The **Semiempirical Extended Tight-Binding Program Package** by the Grimme lab is essential for our calculations. More details can be found at their [GitHub repository](https://github.com/grimme-lab/xtb). 

We deeply appreciate their foundational work.

## Description

- **What it does**: The scripts automate calculations of storage energies, thermal back reaction barriers, and maximum absorption wavelengths using the xTB method.
  
- **Storage**: All results are stored in an SQLite database. Molecules are stored using RDKit as QMmol objects, making it easy to access details like SMILES and XYZ coordinates.
  
- **Languages & Tools**: The scripts are in Python, using SQLite for the database. Calculations run in parallel using SLURM.

- **Order of Operations**: After setting up and populating the database, start by calculating storage energies and ground states for both reactants and products. Next, identify transition states, followed by determining excitation energies and solar conversion efficiency.

- **Submitting to ORCA**: There are specific scripts designed for submitting computations to ORCA, primarily for result comparison. Scripts that end in '_dft' are utilized for this purpose.


## File Descriptions

Each script has a specific role:

- `data_preparation.py`: Prepares the database.
  
- `electronic_conf_search.py`: Generates products from reactant SMILES, finds ground state conformers using xTB, and optimizes them.
  
- `find_max_abs.py`: Uses sTDA to determine max absorption wavelengths from ground state structures.
  
- `find_ts.py`: Finds transition states and back reaction barriers using a reaction path search.
  
- `settings.py`: Contains various settings and paths for the scripts.
  
- `submit_calculations.py`: Manages the submission of molecule batches to SLURM.
  
- `utils.py`: Contains functions for shell commands and extracting energies from xTB outputs.
  
- `db_utils.py`: Functions related to the database.
  
- `generate_azobenzenes.ipynb`: A notebook for creating azobenzene variations.

## Database Structure

SQLite database columns include:

- **HashedName**: A unique name derived from the hashed SMILES string.
- **Smiles**: The SMILES string representing the reactant molecule.
- **Multiplicity**: The molecular multiplicity, detailing the number of unpaired electrons in the molecule.
- **Charge**: The charge on the molecule.
- **BatchID**: A string identifier that groups molecules into batches for parallel calculations.
- **CalculationStage**: A status field that can be updated to indicate the current stage of calculation for a molecule.
- **ProductObject**: A serialized version of the QMmol object that contains data about the product molecule.
- **ReactantObject**: A serialized version of the QMmol object that contains data about the reactant molecule.
- **TransitionStateObject**: A serialized version of the QMmol object that contains data about the transition state.
- **StorageEnergy**: The energy difference between the reactant and the product, helping in the assessment of the storage potential.
- **BackReactionBarrier**: The energy barrier for the reverse reaction, i.e., the energy difference between the transition state and the product.
- **MaxAbsorptionProduct**: The absorption wavelength corresponding to the maximum oscillator strength for the product molecule.
- **MaxOscillatorStrengthProduct**: The maximum oscillator strength for transitions in the product molecule.
- **MaxAbsorptionReactant**: The absorption wavelength corresponding to the maximum oscillator strength for the reactant molecule.
- **MaxOscillatorStrengthReactant**: The maximum oscillator strength for transitions in the reactant molecule.
- **SolarConversionEfficiency**: The solar conversion efficiency calculated from absorptions, storage energy and backreaction barrier.
- **GroundStateStats**: A serialized dictionary containing statistical data (mean, median, standard deviation, variance, and range) of various energy calculations and RMSD values for ground state configurations.
- **TransitionStateStats**: A BLOB object to store statistical data for the transition state configurations (details on obtaining good statistics for this are under consideration).
- **ExcitationStats**: A serialized dictionary containing statistical data on several sTDA calculations, including mean, median, standard deviation, variance, and range of maximum absorptions and oscillator strengths.

Each column gets filled by specific scripts, as detailed in the **Workflow** section.

## Prerequisites

Before you begin, ensure the following dependencies are installed:
(Not sure about this section yet, the yml file is recommended.)
- Python 3.x
- Pandas
- Numpy
- Scipy
- RDKit
- (Others to be determined)

### Setting up the Python Environment

You can find a `screening_env.yml` file in the repository, which helps in setting up a compatible Python environment. To use it, follow these steps:

1. Open your terminal and navigate to the directory where the `screening_env.yml` file is located.
2. Create a new environment by running the command: `conda env create -f screening_env.yml`.
3. Once the environment is created, activate it with: `conda activate screening_env`.

### xTB Program

To run the calculations, you also need to have access to the xTB program. Make sure you have paths set up for both the xTB and the xTB sTDA programs. You will need to specify these paths in the `settings.py` file.

Additionally, ensure that the `.param_stda1.xtb` and `.param_stda2.xtb` files are placed in your home directory. These files are available in this repository and originally stem from the [xTB repository](https://github.com/grimme-lab/xtb). . The precise setup is still under exploration; any suggestions or contributions to solve this are highly appreciated.

## Workflow

1. **Setup**:
   - Install the prerequisites.
   - Adjust `settings.py` for your setup.
   
2. **Input File**:
   - Create a CSV file:
     - `comp_name`: Compound names.
     - `smiles`: SMILES strings.
     - `charge`: Molecule charges.
     - `multiplicity`: Molecule multiplicities.

3. **Database Initialization**:
   - Run: `python data_preparation.py create` to make the database.
   - Populate with: `python data_preparation.py fill $input_file.csv`. You'll be asked for a batch size.

4. **Storage Energies**:
   - Command: `python submit_calculations.py storage`
   - It calculates storage energies using xTB for ground states.
   - Populates: ProductObject, ReactantObject, StorageEnergy, GroundStateStats

5. **Transition States**:
   - Command: `python submit_calculations.py tbr`
   - It finds transition states and back reaction barriers using a reaction path search.
   - Populates: TransitionStateObject, BackReactionBarrier, TransitionStateStats

6. **Absorption & Oscillator Strength**:
   - Command: `python submit_calculations.py abs`
   - It calculates max absorption wavelengths and oscillator strengths using sTDA for excitations.
   - Populates: MaxAbsorptionProduct, MaxOscillatorStrengthProduct, and other related columns.

## Database Functionalities

### Choosing Which Table to Use

In the file `db_utils.py`, the functions `insert_data`, `retrieve_data`, and `update_data` have options to select the table within the database they operate on. By default, they all work on the "MoleculeData" table.

### Command Line Inputs
Manipulating the database is performed through the command line using the script `data_preparation.py`. By executing `python data_preparation.py {command}`, you can access the following functionalities:

**Fill**
- Command: `python data_preparation.py fill $input_file.csv`
- Result: Adds the contents of `input_file.csv` to the database table.
- Prompt: Requests batch size, determining the number of molecules in each calculation batch.

**Create**
- Command: `python data_preparation.py create`
- Result: Creates the database `molecule_data.db` along with the `MoleculeData` table, including the appropriate columns.

**ORCA**
- Command: `python data_preparation.py orca`
- Result: Creates a table named `ORCAData` within the database.

**Change**
- Command: `python data_preparation.py change`
- Result: Manually alters the value of the "CalculationStage" column.
- Prompt: Asks for "CalculationStage:", the input specifies the new value.

**Delete**
- Command: `python data_preparation.py delete`
- Result: Deletes a specified table within the database.
- Prompt: Asks "Delete Table:", input is the name of the table to be deleted.

**Rename**
- Command: `python data_preparation.py rename`
- Result: Renames a table in the database.
- Prompts:
    1. "Old Table Name:" – specify the name of the table to be renamed.
    2. "New Table Name:" – specify the new name for the table.

**Backup**
- Command: `python data_preparation.py backup`
- Result: Automatically creates a copy of the chosen table, appending the current time to the name (YYYYMMDD_HHMMSS).
- Prompt: "Which table:" – specify the table name for backup.

**Add Batch**
- Command: `python data_preparation.py add_batch $extra_input.csv $new_batch_id`
- Result: Adds an extra batch with a specific batch ID to the end of the database table.


### Retrieve Data to a Pandas DataFrame

The `db_utils.py` file contains a function `retrieve_data` that constructs a SQLite query based on two arguments. The first argument, `filters`, is a dictionary where the keys represent column names and the values specify filter criteria. The second argument determines which database table to source the data from. This function returns a pandas DataFrame containing all columns from the table, filtered to include rows that match the criteria. Additionally, it automatically deserializes columns containing objects and dictionaries, making them accessible in Python.

Here's an example retrieving molecules with the string "storage_completed" in the "CalculationStage" column from a table named "TestMolecules":

```python
filters = {'CalculationStage': '%storage_completed%'}
molecules_df = retrieve_data(filters, db_table='TestMolecules')
```

#### Getting XYZ Files Using the QMC Package

To obtain XYZ files using the QMC package, start by fetching the QMmol or QMconf object:

```python
reactant = molecules_df['ReactantObject']
reactant.write_xyz()
```

The resulting file will be named according to the label of the reactant (e.g., `{reactant.label}.xyz`), which is stored in the object.

#### Using RDKit

To convert the QMmol or QMconf object to an XYZ file using RDKit, follow these steps:

1. Convert the object to an RDKit molecule:

   ```python
   reactant_rdkit_mol = reactant.get_rdkit_mol()
   ```

2. Then, use RDKit to create an XYZ file, choosing your own filename:

   ```python
   from rdkit.Chem import rdmolfiles

   xyz_data = rdmolfiles.MolToXYZBlock(reactant_rdkit_mol)
   with open(filename, 'w') as f:
       f.write(xyz_data)
   ```


## Submitting Calculations

### Command Line Input

Generally, use `python submit_calculations.py {calculation_type}`. This fetches data from the database that has reached the necessary calculation stage for the specified type. 

**Resubmit a Batch**
- Command: `python submit_calculations.py {calculation_type} {batch_id}`
- Note: *{batch_id}* is optional. It can be added to the command if a specific batch needs to be resubmitted.

#### Calculation Types

**Storage Energy Calculations**
- Command: `python submit_calculations.py storage`
- Calculation Stage Needed: "storage"

**Transition State Search**
- Command: `python submit_calculations.py tbr`
- Calculation Stage Needed: "storage_completed"

**Excitation Energies and SCE**
- Command: `python submit_calculations.py abs`
- Calculation Stages Needed: "storage_completed, ts_completed"

### ORCA

To submit to ORCA, modify the settings by uncommenting `SLURM_TEMPLATE = SLURM_TEMPLATE_ORCA` and adding "_dft" to the calculation types. Note: This feature is not fully tested and may have some issues.







Please reach out for questions or feedback `:)`
