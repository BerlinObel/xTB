import time
import os
import sys
import textwrap

import pandas as pd
import numpy as np

sys.path.append("/groups/kemi/brq616/speciale/opt/xTB/QMC/QMC")
from qmconf import QMConf


def qsub_prep(batchname, script_path, cpus, mem):
    """ """
    pwd = os.getcwd()

    qsub_file = '''\
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
    /groups/kemi/ree/anaconda3/envs/my-rdkit-env/bin/python {4} {0}.pkl {1}
    #/groups/kemi/koerstz/anaconda3/envs/quantum_ml/bin/python {4} {0}.pkl {1}
    # copy data back
    cp *_abs_prod.pkl {3}

    '''.format(batchname, cpus, mem, pwd, script_path)

    with open(batchname + "_qsub.tmp", 'w') as qsub:
        qsub.write(textwrap.dedent(qsub_file))

    return batchname + "_qsub.tmp"

def submit_job(pkl, script, mem, cpus, nodes):

    qsub_name = qsub_prep(pkl, script, cpus, mem)

    batch_id = os.popen("sbatch " + qsub_name).read()
    batch_id = batch_id.strip().split()[-1]

    return int(batch_id)


def run_calculations(pkl_names, script, mem, cpus, nodes):
    ''' '''

    submitted_jobs = set()
    for pkl in pkl_names:
        batch_id = submit_job(pkl, script, mem, cpus, nodes)
        submitted_jobs.add(batch_id)

        if len(submitted_jobs) >= nodes:

            while True:
                output = os.popen("squeue -u obel -p chem").readlines()[1:]
                all_running_jobs = set([int(job.split()[0]) for job in output])

                if len(all_running_jobs & submitted_jobs) >= nodes: # intersect
                    time.sleep(10)
                else:
                    # remove finished jobs
                    finished_jobs = submitted_jobs - all_running_jobs
                    for job in finished_jobs:
                        submitted_jobs.remove(job)
                    break


if __name__ == "__main__":

    # input params
    cpus = 2
    mem = "8GB"

    nodes = 200

    script = '/groups/kemi/brq616/speciale/opt/xTB/scripts/find_max_abs_prod.py'

    data_file = sys.argv[1] #a csv file with all the batch names

    ##########################################################################
    #
    ##########################################################################
    
    file_names = np.loadtxt(sys.argv[1], dtype=str) #a csv file with all the batch names without .pkl
    file_names = file_names.tolist()

    # run calculations on nodes
    run_calculations(file_names, script, mem, cpus, nodes)
