#!/bin/bash

logdir=peloton/log
mkdir -p $logdir

source modules.pel

SBATCH="sbatch -N 1 -n {threads}"
SBATCH="$SBATCH -t {cluster.time} --mem={cluster.mem}"
SBATCH="$SBATCH -o peloton/log/o-%j -e peloton/log/e-%j"

snakemake                                    \
    -j 1000                                  \
    --cluster-config peloton/cluster.yaml    \
    --js peloton/jobscript.sh                \
    --rerun-incomplete                       \
    --keep-going                             \
    --latency-wait 10                        \
    --use-conda                              \
    --cluster "$SBATCH" $@
