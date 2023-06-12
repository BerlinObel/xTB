import time
import os
import sys
import textwrap

import pandas as pd
import numpy as np

storage_energy_script = '/groups/kemi/brq616/speciale/opt/xTB/scripts/electronic_azo_conf_search.py'
transistion_state_script = '/groups/kemi/brq616/speciale/opt/xTB/scripts/find_ts.py'
absorption_script = '/groups/kemi/brq616/speciale/opt/xTB/scripts/find_max_abs.py'

def prepare_slurm_script(job_name, python_script, cpus, memory, working_dir, file_suffix, job_type):
    """Creates a SLURM script to submit a calculation to the queue."""
    
    if job_type == "storage":
        slurm_script = '''\
        #!/bin/sh
        #SBATCH --job-name={0}
        #SBATCH --cpus-per-task={1}
        #SBATCH --mem={2}
        #SBATCH --ntasks=1
        #SBATCH --error={3}/{0}.stderr
        #SBATCH --output={3}/{0}.stdout
        #SBATCH --time=1000:00:00
        #SBATCH --partition=chem
        #SBATCH --no-requeue

        cd /scratch/$SLURM_JOB_ID

        echo "Job id: $SLURM_JOB_ID"

        # copy batch file
        cp {3}/{0}.csv .

        # run python code
        python {4} {0}.csv {2}

        # copy data back
        cp {0}.pkl {3}

    '''.format(job_name, cpus, memory, working_dir, python_script)
    
    else:
        slurm_script = '''\
        #!/bin/sh
        #SBATCH --job-name={0}
        #SBATCH --cpus-per-task={1}
        #SBATCH --mem={2}
        #SBATCH --ntasks=1
        #SBATCH --error={3}/{0}.stderr
        #SBATCH --output={3}/{0}.stdout
        #SBATCH --time=1000:00:00
        #SBATCH --partition=chem
        #SBATCH --no-requeue

        cd /scratch/$SLURM_JOB_ID

        echo "Job id: $SLURM_JOB_ID"

        # copy batch file
        cp {3}/{0}.pkl .

        # run python code
        python {4} {0}.pkl {2}

        # copy data back
        cp *{0}{5}.pkl {3}

        '''.format(job_name, cpus, memory, working_dir, python_script, file_suffix)

    slurm_script_filename = f"{job_name}_qsub.tmp"

    with open(slurm_script_filename, 'w') as slurm_file:
        slurm_file.write(textwrap.dedent(slurm_script))

    return slurm_script_filename

def submit_slurm_job(slurm_script):
    """Submits a job to SLURM and returns the job ID."""
    batch_id = os.popen("sbatch " + slurm_script).read()
    return batch_id.strip().split()[-1]

def run_calculations(batch_files, python_script, memory, cpus, max_running_jobs, file_suffix, job_type):
    """Runs a batch of calculations on a SLURM cluster."""
    submitted_jobs = set()
    for batch_file in batch_files:
        slurm_script = prepare_slurm_script(batch_file, python_script, cpus, memory, os.getcwd(), file_suffix, job_type)
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
    

    if sys.argv[2] == "storage":
        python_script = storage_energy_script
        file_suffix = ""
        chunk_size=int(input("Batch size: "))
        data = pd.read_csv(data_file)
        chunks = [data[i:i+chunk_size] for i in range(0, data.shape[0], chunk_size)]
        batch_files = []
        for idx, chunk in enumerate(chunks):
            chunk_name = f"smiles_batch-{idx}"
            chunk.to_csv(chunk_name + ".csv", index=False)
            batch_files.append(chunk_name)

    if sys.argv[2] == "tbr":
        python_script = transistion_state_script
        file_suffix = "_new"
        batch_files = np.loadtxt(data_file, dtype=str).tolist()

    if sys.argv[2] == "abs":
        python_script = absorption_script
        file_suffix = "_abs"
        batch_files = np.loadtxt(data_file, dtype=str).tolist()
        
    
    run_calculations(batch_files, python_script, memory, cpus, max_running_jobs, file_suffix, sys.argv[2])

        

    
