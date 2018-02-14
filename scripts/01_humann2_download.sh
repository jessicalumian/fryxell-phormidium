#!/bin/bash -l
#
#SBATCH --job-name=humann2-bin2-1_download_dbs.sh
#SBATCH -t 10:00:00
#SBATCH --mem=6G
#SBATCH -o /home/jemizzi/phormidium/output/humann2_scripts/01_download.out
#SBATCH -e /home/jemizzi/phormidium/output/humann2_scripts/01_download.err 
#SBATCH --mail-user=jemizzi@ucdavis.edu
#SBATCH --mail-type=ALL

set -e

module load bio

# download databases and utility mapping (KEGG legacy mapping)

humann2_databases --download chocophlan full /home/jemizzi/rotation-project/databases/humann2-database/
humann2_databases --download uniref uniref50_diamond /home/jemizzi/rotation-project/databases/humann2-database/
humann2_databases --download utility_mapping full /home/jemizzi/rotation-project/databases/humann2-database/

