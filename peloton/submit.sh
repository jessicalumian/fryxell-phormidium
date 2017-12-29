logdir=phormidium/log
mkdir -p $logdir

SUB="sbatch -p shared -A m342 -N 1 -n {threads} "
SUB="$SUB -t {cluster.time} --mem={cluster.mem}"
SUB="$SUB -o $logdir -e $logdir"

snakemake                                    \
    -j 3000                                  \
    --cluster-config phormidium/cluster.yaml \
    --js phormidium/jobscript.sh             \
    --rerun-incomplete                       \
    --keep-going                             \
    --latency-wait 10                        \
    --use-conda                              \
    --cluster "$SUB" $@
