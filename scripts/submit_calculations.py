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
from db_utils import insert_data, connect_db, close_db
from settings import SLURM_TEMPLATE, JOB_TYPE_CONFIG, USER, QUEUE


def prepare_slurm_script(batch_name, job_name, python_script, cpus, memory, working_dir):
    """Creates a SLURM script to submit a calculation to the queue."""
    slurm_script = SLURM_TEMPLATE.format(batch_name = batch_name,
        job_name=job_name, cpus=cpus, memory=memory, working_dir=working_dir, python_script=python_script, queue=QUEUE,
    )

    slurm_script_filename = f"{job_name}_qsub.tmp"

    with open(slurm_script_filename, 'w') as slurm_file:
        slurm_file.write(textwrap.dedent(slurm_script))

    return slurm_script_filename

def submit_slurm_job(slurm_script):
    """Submits a job to SLURM and returns the job ID."""
    batch_id = os.popen("sbatch " + slurm_script).read()
    return batch_id.strip().split()[-1] 

def run_calculations(batch_files: List[str], python_script: str, memory: int, cpus: int, max_running_jobs: int, job_suffix: str):
    """Runs a batch of calculations on a SLURM cluster."""
    submitted_jobs = set()
    for batch_file in batch_files:
        job_name = f"{batch_file}{job_suffix}"
        batch_name = batch_file
        slurm_script = prepare_slurm_script(batch_name, job_name, python_script, cpus, memory, os.getcwd())
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
    cpus = 8
    memory = "8GB"
    max_running_jobs = 200

    job_type = sys.argv[1]
    
    if job_type not in JOB_TYPE_CONFIG:
        print(f"Unknown job type: {job_type}")
        sys.exit(1)
        
    config = JOB_TYPE_CONFIG[job_type]
    # Mapping of job types to the calculation stages that should be used to find batches ready for processing
    stage_mapping = {
        'storage': 'storage',
        'tbr': 'storage_completed',
        'abs': 'storage_completed',
        # ... (add other job types as needed)
    }
    # Query the database to get a list of distinct batch IDs that are at the appropriate stage for this job type
    conn = connect_db()
    cursor = conn.cursor()
    cursor.execute("SELECT DISTINCT BatchID FROM MoleculeData WHERE CalculationStage = ?", (stage_mapping[job_type],))
    batch_files = [row[0] for row in cursor.fetchall()]
    close_db(conn)
    
    run_calculations(batch_files, config['python_script'], memory, cpus, max_running_jobs, config['job_suffix'])
    
