import os
import queue
import sys
import time
import pandas as pd
import numpy as np
import textwrap
from typing import List, Set


# Add the directory of the script to the Python path
sys.path.append(os.path.dirname(os.path.realpath(__file__)))
from settings import SLURM_TEMPLATE, JOB_TYPE_CONFIG, USER, QUEUE

def prepare_slurm_script(job_name, python_script, cpus, memory, working_dir, file_suffix):
    """Creates a SLURM script to submit a calculation to the queue."""
    job_suffix = '.csv' if file_suffix == '' else '.pkl'
    slurm_script = SLURM_TEMPLATE.format(
        job_name=job_name, cpus=cpus, memory=memory, working_dir=working_dir,
        job_suffix=job_suffix, python_script=python_script, file_suffix=file_suffix, queue=QUEUE,
    )

    slurm_script_filename = f"{job_name}_qsub.tmp"

    with open(slurm_script_filename, 'w') as slurm_file:
        slurm_file.write(textwrap.dedent(slurm_script))

    return slurm_script_filename

def submit_slurm_job(slurm_script):
    """Submits a job to SLURM and returns the job ID."""
    batch_id = os.popen("sbatch " + slurm_script).read()
    return batch_id.strip().split()[-1] 

def run_calculations(batch_files: List[str], python_script: str, memory: int, cpus: int, max_running_jobs: int, file_suffix: str):
    """Runs a batch of calculations on a SLURM cluster."""
    submitted_jobs = set()
    for batch_file in batch_files:
        slurm_script = prepare_slurm_script(batch_file, python_script, cpus, memory, os.getcwd(), file_suffix)
        batch_id = submit_slurm_job(slurm_script)
        submitted_jobs.add(batch_id)

        manage_job_submission(max_running_jobs, submitted_jobs)

def manage_job_submission(max_running_jobs: int, submitted_jobs: Set[str]) -> None:
    """
    Manages job submissions based on pre-defined constraints.  
    This function ensures that the number of submitted jobs does 
    not exceed the pre-defined limit on nodes ('max_running_jobs'). 
    New jobs are only submitted when there are sufficient system resources 
    available as per this limit.
    """
    while len(submitted_jobs) >= max_running_jobs:
        output = os.popen(f"squeue -u {USER} -p {QUEUE}").readlines()[1:]
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


    
