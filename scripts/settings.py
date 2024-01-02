import numpy as np
import os
from typing import Dict, Optional

# Path to database
# DB_PATH = "/groups/kemi/obel/github/xTB/molecule_data.db"
DB_PATH = "/groups/kemi/brq616/speciale/molecule_data.db"

# Base path - replace with the path to the root directory of your project
# BASE_PATH = "/groups/kemi/obel/github/xTB/"
BASE_PATH = "/groups/kemi/brq616/speciale/opt/xTB/"

#USER = "obel"
USER = "brq616"

QUEUE = "kemi6"
# QUEUE = "chem"

QMC_PATH = os.path.join(BASE_PATH, "tQMC/QMC")

## Path to xtb programs

XTB4STDA_PATH = "/groups/kemi/brq616/speciale/opt/stda/xtb4stda"
STDA_PATH = "/groups/kemi/brq616/speciale/opt/stda/stda_v1.6.1"

## Calculation script paths

STORAGE_ENERGY_SCRIPT = os.path.join(BASE_PATH, "scripts/electronic_conf_search.py")
TRANSITION_STATE_SCRIPT = os.path.join(BASE_PATH, "scripts/find_ts.py")
ABSORPTION_SCRIPT = os.path.join(BASE_PATH, "scripts/find_max_abs.py")
TRANSITION_STATE_SCRIPT_ORCA = os.path.join(BASE_PATH, "scripts/find_ts_orca.py")
DFT_GEOMETRY_OPTIMIZATION_SCRIPT = os.path.join(BASE_PATH, "scripts/ground_state_search_dft.py")
DFT_TRANSISTION_STATE_SCRIPT = os.path.join(BASE_PATH, "scripts/find_ts_dft.py")
DFT_MAX_ABS_SCRIPT = os.path.join(BASE_PATH, "scripts/ground_state_search_dft.py")

## SLURM template and job types for submit_calculations.py

JOB_TYPE_CONFIG: Dict[str, Dict[str, Optional[str]]] = {
    'storage': {
        'python_script': STORAGE_ENERGY_SCRIPT,
        'job_suffix': '_energy',
    },
    'tbr': {
        'python_script': TRANSITION_STATE_SCRIPT,
        'job_suffix': '_TS',
    },
    'abs': {
        'python_script': ABSORPTION_SCRIPT,
        'job_suffix': '_abs',
    },
    'tbr_orca': {
        'python_script': TRANSITION_STATE_SCRIPT_ORCA,
        'job_suffix': '_ts_orca'
    },
    'storage_dft': {
        'python_script': DFT_GEOMETRY_OPTIMIZATION_SCRIPT,
        'job_suffix': '_energy_dft',
    },
    'tbr_dft': {
        'python_script': DFT_TRANSISTION_STATE_SCRIPT,
        'job_suffix': '_ts_dft'
    },
    'abs_dft': {
        'python_script': DFT_MAX_ABS_SCRIPT,
        'job_suffix': '_abs_dft',
    }
}

SLURM_TEMPLATE_XTB = '''\
#!/bin/sh
#SBATCH --job-name={job_name}
#SBATCH --cpus-per-task={cpus}
#SBATCH --tasks-per-node={cpus}
#SBATCH --mem={memory}
#SBATCH --ntasks=1
#SBATCH --error={working_dir}/{job_name}.stderr
#SBATCH --output={working_dir}/{job_name}.stdout
#SBATCH --time=1000:00:00
#SBATCH --partition={queue}
#SBATCH --no-requeue

cd /scratch/$SLURM_JOB_ID

echo "Job id: $SLURM_JOB_ID"

# run python code
python -u {python_script} {batch_name} {memory}

'''
SLURM_TEMPLATE_ORCA = '''\
#!/bin/sh
#SBATCH --job-name={job_name}
#SBATCH --cpus-per-task={cpus}
#SBATCH --tasks-per-node={cpus}
#SBATCH --mem={memory}
#SBATCH --ntasks=1
#SBATCH --error={working_dir}/{job_name}.stderr
#SBATCH --output={working_dir}/{job_name}.stdout
#SBATCH --time=1000:00:00
#SBATCH --partition={queue}
#SBATCH --no-requeue

# Set environment variables for ORCA
export PATH="/software/kemi/Orca/orca_5_0_1_linux_x86-64_openmpi411:/software/kemi/openmpi/openmpi-4.1.1/bin:$PATH"
export LD_LIBRARY_PATH="/software/kemi/openmpi/openmpi-4.1.1/lib"
export ORCA="/software/kemi/Orca/orca_5_0_1_linux_x86-64_openmpi411/orca"

cd /scratch/$SLURM_JOB_ID

echo "Job id: $SLURM_JOB_ID"

# run python code
python -u {python_script} {batch_name} {memory}

'''

# SLURM_TEMPLATE = SLURM_TEMPLATE_ORCA
SLURM_TEMPLATE = SLURM_TEMPLATE_XTB

REACTION_PATH_TEMPLATE = '''$path
nrun=6
npoint=35
anopt=10
kpush=0.003
kpull=-0.015
ppull=0.05
alp=0.7
$end
$md
 restart=true
$end
'''

# Constants
AVOGADROS_NUMBER = 6.02214199e23
SPEED_OF_LIGHT = 299792458
ELECTRON_CHARGE = 1.60217662e-19
ELECTRON_MASS = 9.10938e-31
PERMITTIVITY_OF_VACUUM = 8.8541878176e-12

# Sigma (width of the Gaussian) in electron volts is a parameter defining the width of 
# the Gaussian distribution that will be used to model the UV-Vis absorption spectrum.
SIGMA_ELECTRON_VOLT = 0.4

# Convert the sigma value from electron volts to cm^-1.
SIGMA_CM = SIGMA_ELECTRON_VOLT * 8065.544

# FACTOR converts the oscillator strengths to absorption values, enabling a comparison with experimental UV-Vis spectra.
FACTOR = (AVOGADROS_NUMBER * ELECTRON_CHARGE**2) / (np.log(10) * 2 * ELECTRON_MASS * SPEED_OF_LIGHT**2 * PERMITTIVITY_OF_VACUUM) * np.sqrt(np.log(2)/np.pi) * 1e-1

# Define a range of wavelengths over which to calculate the absorption spectrum. Here based on visible light
WAVELENGTH_RANGE = np.linspace(125, 525, 593)