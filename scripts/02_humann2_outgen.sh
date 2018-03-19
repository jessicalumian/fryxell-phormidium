#!/bin/bash -l
#
#SBATCH --job-name=02_humann2_outgen_bin20
#SBATCH -t 10:00:00
#SBATCH --mem=12G
#SBATCH -o /home/jemizzi/phormidium/output/humann2_scripts/02_outgen_bin20.out
#SBATCH -e /home/jemizzi/phormidium/output/humann2_scripts/02_outgen_bin20.err
#SBATCH --mail-user=jemizzi@ucdavis.edu
#SBATCH --mail-type=ALL

set -e

module load python usearch samtools perlbrew bowtie2 metaphlan/2.0

humann2 --input /home/jemizzi/phormidium/output/04_anvio/merged_mat_lab_coassembly/SAMPLES-SUMMARY/bin_by_bin/Bin_20/Bin_20-contigs.fa --nucleotide-database /home/jemizzi/rotation-project/databases/humann2-database/chocophlan --protein-database /home/jemizzi/rotation-project/databases/humann2-database/uniref --output /home/jemizzi/phormidium/output/humann2_scripts

