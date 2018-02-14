#!/bin/bash -l
#
#SBATCH --job-name=humann2-bin20-4-map-genes.sh
#SBATCH -t 5:00:00
#SBATCH --mem=6G
#SBATCH -o /home/jemizzi/phormidium/output/humann2_scripts/humann2-bin20-4-map-genes.out
#SBATCH -e /home/jemizzi/phormidium/output/humann2_scripts/humann2-bin20-4-map-genes.err
#SBATCH --mail-user=jemizzi@ucdavis.edu
#SBATCH --mail-type=ALL

set -e

module load bio

# map gene family CPMs to pathways

humann2_regroup_table --input /home/jemizzi/phormidium/output/humann2_scripts/Bin_20-contigs_genefamilies.tsv --groups uniref50_ko --output /home/jemizzi/phormidium/output/humann2_scripts/Bin_20-contigs_genefamilies_uniref50_ko.tsv

