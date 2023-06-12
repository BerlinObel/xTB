import os
import sys
import time
import pandas as pd
import numpy as np
import textwrap

STORAGE_ENERGY_SCRIPT = '/groups/kemi/brq616/speciale/opt/xTB/scripts/electronic_azo_conf_search.py'
TRANSITION_STATE_SCRIPT = '/groups/kemi/brq616/speciale/opt/xTB/scripts/find_ts.py'
ABSORPTION_SCRIPT = '/groups/kemi/brq616/speciale/opt/xTB/scripts/find_max_abs.py'

JOB_TYPE_CONFIG = {
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

def prepare_slurm_script(job_name, python_script, cpus, memory, working_dir, file_suffix):
    """Creates a SLURM script to submit a calculation to the queue."""
    job_suffix = '.csv' if file_suffix == '' else '.pkl'
    slurm_script = SLURM_TEMPLATE.format(
        job_name=job_name, cpus=cpus, memory=memory, working_dir=working_dir,
        job_suffix=job_suffix, python_script=python_script, file_suffix=file_suffix,
    )

    slurm_script_filename = f"{job_name}_qsub.tmp"

    with open(slurm_script_filename, 'w') as slurm_file:
        slurm_file.write(textwrap.dedent(slurm_script))

    return slurm_script_filename

def submit_slurm_job(slurm_script):
    """Submits a job to SLURM and returns the job ID."""
    batch_id = os.popen("sbatch " + slurm_script).read()
    return batch_id.strip().split()[-1]

def run_calculations(batch_files, python_script, memory, cpus, max_running_jobs, file_suffix):
    """Runs a batch of calculations on a SLURM cluster."""
    submitted_jobs = set()
    for batch_file in batch_files:
        slurm_script = prepare_slurm_script(batch_file, python_script, cpus, memory, os.getcwd(), file_suffix)
        batch_id = submit_slurm_job(slurm_script)
        submitted_jobs.add(batch_id)

        while len(submitted_jobs) >= max_running_jobs:
            output = os.popen("squeue -u brq616 -p chem").readlines()[1:]
            running_jobs = {job.split()[0] for job in output}

            # remove finished jobs
            submitted_jobs -= running_jobs

            if len(submitted_jobs) < max_running_jobs:
                break
            time.sleep(10)

if __name__ == "__main__":
    cpus = 2
    memory = "8GB"
    max_running_jobs = 200

    data_file = sys.argv[1]
    job_type = sys.argv[2]

    if job_type not in JOB_TYPE_CONFIG:
        print(f"Unknown job type: {job_type}")
        sys.exit(1)

    config = JOB_TYPE_CONFIG[job_type]

    if job_type == "storage":
        config['chunk_size'] = int(input("Batch size: "))
        data = pd.read_csv(data_file)
        chunks = [data[i:i+config['chunk_size']] for i in range(0, data.shape[0], config['chunk_size'])]
        batch_files = []
        for idx, chunk in enumerate(chunks):
            chunk_name = f"smiles_batch-{idx}"
            chunk.to_csv(chunk_name + ".csv", index=False)
            batch_files.append(chunk_name)
    else:
        batch_files = np.loadtxt(data_file, dtype=str).tolist()

    run_calculations(batch_files, config['python_script'], memory, cpus, max_running_jobs, config['file_suffix'])


    
