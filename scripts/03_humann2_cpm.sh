#!/bin/bash -l
#
#SBATCH --job-name=humann2-bin2-3-cpm.sh
#SBATCH -t 03:00:00
#SBATCH --mem=6G
#SBATCH -o /home/jemizzi/phormidium/output/humann2_scripts/03_cpm_bin20.out
#SBATCH -e /home/jemizzi/phormidium/output/humann2_scripts/03_cpm_bin20.err
#SBATCH --mail-user=jemizzi@ucdavis.edu
#SBATCH --mail-type=ALL

set -e

module load bio

# convert output to CPM

humann2_renorm_table --input /home/jemizzi/phormidium/output/humann2_scripts/Bin_20-contigs_genefamilies.tsv --units "cpm" --output /home/jemizzi/phormidium/output/humann2_scripts/Bin_20-contigs_genefamilies_cpm.tsv

