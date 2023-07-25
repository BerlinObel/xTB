import numpy as np
from typing import Dict, Optional

USER = "brq616"
QUEUE = "kemi6"
QMC_PATH = "/groups/kemi/brq616/speciale/opt/xTB/tQMC/QMC"


## Path to xtb programs

XTB4STDA_PATH = "/groups/kemi/brq616/speciale/opt/stda/xtb4stda"
STDA_PATH = "/groups/kemi/brq616/speciale/opt/stda/stda_v1.6.1"

## Calculation script paths

STORAGE_ENERGY_SCRIPT = '/groups/kemi/brq616/speciale/opt/xTB/scripts/electronic_azo_conf_search.py'
TRANSITION_STATE_SCRIPT = '/groups/kemi/brq616/speciale/opt/xTB/scripts/find_ts.py'
ABSORPTION_SCRIPT = '/groups/kemi/brq616/speciale/opt/xTB/scripts/find_max_abs.py'


## SLURM template and job types for submit_calculations.py

JOB_TYPE_CONFIG: Dict[str, Dict[str, Optional[str]]] = {
    'storage': {
        'python_script': STORAGE_ENERGY_SCRIPT,
        'file_suffix': '',
        'chunk_size': None,  # This will be filled in based on user input
    },
    'tbr': {
        'python_script': TRANSITION_STATE_SCRIPT,
        'file_suffix': '_new',
    },
    'abs': {
        'python_script': ABSORPTION_SCRIPT,
        'file_suffix': '_abs',
    },
}

SLURM_TEMPLATE = '''\
#!/bin/sh
#SBATCH --job-name={job_name}
#SBATCH --cpus-per-task={cpus}
#SBATCH --mem={memory}
#SBATCH --ntasks=1
#SBATCH --error={working_dir}/{job_name}.stderr
#SBATCH --output={working_dir}/{job_name}.stdout
#SBATCH --time=1000:00:00
#SBATCH --partition=chem
#SBATCH --no-requeue

cd /scratch/$SLURM_JOB_ID

echo "Job id: $SLURM_JOB_ID"

# copy batch file
cp {working_dir}/{job_name}{job_suffix} .

# run python code
python {python_script} {job_name}{job_suffix} {memory}

# copy data back
cp *{job_name}{file_suffix}.pkl {working_dir}
'''

REACTION_PATH_TEMPLATE = '''$path
nrun=1
npoint=25
anopt=10
kpush=0.003
kpull=-0.015
ppull=0.05
alp=0.9
$end'''

## 2 runs
# '''$path
# nrun=2
# npoint=25
# anopt=10
# kpush=0.003
# kpull=-0.015
# ppull=0.05
# alp=0.9
# $end'''

## everything
# '''$path
# nrun=2
# npoint=30
# anopt=10
# kpush=0.001
# kpull=-0.01
# ppull=0.04
# alp=0.9
# $end'''

## 2 runs and 30 points
# '''$path
# nrun=2
# npoint=30
# anopt=10
# kpush=0.003
# kpull=-0.015
# ppull=0.05
# alp=0.9
# $end'''

## Lower push pull and 30 points
# '''$path
# nrun=1
# npoint=30
# anopt=10
# kpush=0.001
# kpull=-0.01
# ppull=0.04
# alp=0.9
# $end'''

## Lower push pull plus two runs
# '''$path
# nrun=2
# npoint=25
# anopt=10
# kpush=0.001
# kpull=-0.01
# ppull=0.04
# alp=0.9
# $end'''

## Lower push pull
# '''$path
# nrun=1
# npoint=25
# anopt=10
# kpush=0.001
# kpull=-0.01
# ppull=0.04
# alp=0.9
# $end'''

## 30 points
# '''$path
# nrun=1
# npoint=30
# anopt=10
# kpush=0.003
# kpull=-0.015
# ppull=0.05
# alp=0.9
# $end'''

## Original
# '''$path
# nrun=1
# npoint=25
# anopt=10
# kpush=0.003
# kpull=-0.015
# ppull=0.05
# alp=0.9
# $end'''

## Constants and settings for absorbtion calculations in "find_max_abs.py"

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