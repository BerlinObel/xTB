import time
import os
import sys
import textwrap

import pandas as pd
import numpy as np

sys.path.append("./tQMC/QMC")
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
    echo "{0}{1}{2}{3}{4}"
    
    # copy batch file
    cp {3}/{0}.csv .

    # run python code
    /groups/kemi/koerstz/anaconda3/envs/quantum_ml/bin/python {4} {0}.csv {1}
    #/groups/kemi/ree/anaconda3/envs/my-rdkit-env/bin/python {4} {0}.csv {1}
    echo "/groups/kemi/koerstz/anaconda3/envs/quantum_ml/bin/python {4} {0}.csv {1}"
    # copy data back
    #tar -czf {0}_output.tar.gz *xyz
    #cp {0}_output.tar.gz *csvl {3}
    #cp *out {3}

    cp *pkl {3}

    '''.format(batchname, cpus, mem, pwd, script_path)

    with open(batchname + "_qsub.tmp", 'w') as qsub:
        qsub.write(textwrap.dedent(qsub_file))

    return batchname + "_qsub.tmp"

def submit_job(csv, script, mem, cpus, nodes):

    qsub_name = qsub_prep(csv, script, cpus, mem)

    batch_id = os.popen("sbatch " + qsub_name).read()
    batch_id = batch_id.strip().split()[-1]

    return int(batch_id)


def run_calculations(csv_names, script, mem, cpus, nodes):
    ''' '''

    submitted_jobs = set()
    for csv in csv_names:
        batch_id = submit_job(csv, script, mem, cpus, nodes)
        submitted_jobs.add(batch_id)

        if len(submitted_jobs) >= nodes:

            while True:
                output = os.popen("squeue -u brq616 -p chem").readlines()[1:]
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

    nodes = 100
    
    if sys.argv[2]:
        chunk_size=int(sys.argv[2])
    else:
        chunk_size = 1000
    test = False
    if sys.argv[3] == "test":
        test = True
    if test:
        script = './xTB/scripts/electronic_azo_conf_search.py'
    else:
        script = './xTB/scripts/electronic_azo_conf_search.py'

    data_file = sys.argv[1]

    ##########################################################################
    #
    ##########################################################################

    # import data
    data = pd.read_csv(data_file)

    # split data into chunks
    chunked_data = [data[i:i+chunk_size] for i in range(0, data.shape[0], chunk_size)]

    chunk_names = list()
    for idx, chunk in enumerate(chunked_data):
        chunk_name = "smiles_batch-{}".format(idx)
        chunk.to_csv(chunk_name + ".csv", index=False)

        chunk_names.append(chunk_name)

    # run calculations on nodes
    run_calculations(chunk_names, script, mem, cpus, nodes)
