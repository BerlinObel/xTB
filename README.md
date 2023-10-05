# Automated xTB Scripts

Python scripts to automate calculations in computational chemistry using the xTB method.

## Acknowledgement (Not Finished)

This project is built on foundational scripts and research by [Add more names and repositories here]. 

- The **QMC package**, originally by Mads Koerstz, is used here with some changes. It helps in storing molecular shapes and calculation outputs. Check out the original at [QMC Repository](https://github.com/koerstz/QMC).
  
- The **Semiempirical Extended Tight-Binding Program Package** by the Grimme lab is essential for our calculations. More details can be found at their [GitHub repository](https://github.com/grimme-lab/xtb). 

We deeply appreciate their foundational work.

## Description

- **What it does**: The scripts automate calculations of storage energies, thermal back reaction barriers, and maximum absorption wavelengths using the xTB method.
  
- **Storage**: All results are stored in an SQLite database, which is efficient for both storage and retrieval. Molecules are stored using RDKit as QMmol objects, making it easy to access details like SMILES and XYZ coordinates.
  
- **Languages & Tools**: The scripts are in Python, using SQLite for the database. Calculations run in parallel using SLURM.

- **Progress Tracking**: Once the database is set up and filled, you can calculate storage energies and ground states for reactants and products. Subsequent calculations can be submitted after that. We're also working on gathering statistics from each calculation to understand uncertainties better.

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

Please reach out for questions or feedback `:)`
