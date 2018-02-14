#!/bin/bash -l
#
#SBATCH --job-name=checkm.sh
#SBATCH -t 4:00:00
#SBATCH --mem=57344
#SBATCH -o /home/jemizzi/phormidium/output/checkm_coassembly_bins/checkm.out
#SBATCH -e /home/jemizzi/phormidium/output/checkm_coassembly_bins/checkm.err
#SBATCH --mail-user=jemizzi@ucdavis.edu
#SBATCH --mail-type=ALL

module load bio
module load pplacer

# folder which results go in must be empty

rm -fr /home/jemizzi/phormidium/output/checkm_coassembly_bins/checkm_lineage

# reduced trees saves memory

checkm lineage_wf --reduced_tree -x fa /home/jemizzi/phormidium/output/checkm_coassembly_bins/bins /home/jemizzi/phormidium/output/checkm_coassembly_bins/checkm_lineage
