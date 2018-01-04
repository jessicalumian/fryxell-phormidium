#!/bin/bash

#SBATCH -t 4:00:00
#SBATCH --mem=1G

source modules.pel

logdir=log
mkdir -p $logdir

SBATCH="sbatch -l nodes=1:ppn={threads} "
SBATCH="$SBATCH -l walltime={cluster.time} -l mem={cluster.mem}"
SBATCH="$SBATCH -o $logdir -e $logdir peloton/log"

snakemake --unlock

snakemake                                  \
    -j 8                                   \
    --local-cores 4                        \
    --cluster-config peloton/cluster.yaml  \
    --js peloton/jobscript.sh              \
    --rerun-incomplete                     \
    --cluster "$SBATCH"                    \
    --latency-wait 10                      \
    --use-conda                            \
    >>peloton/log/snakemake.log 2>&1
