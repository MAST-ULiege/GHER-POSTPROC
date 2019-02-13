#!/bin/bash
# Submission script for NIC4 
#SBATCH --job-name=PyAtWork
#SBATCH --time=05:00:00 # hh:mm:ss
#
#SBATCH --ntasks=1 
#SBATCH --mem-per-cpu=20000 # megabytes 
#SBATCH --partition=defq 
#
#SBATCH --mail-user=acapet@ulg.ac.be
#SBATCH --mail-type=ALL
#SBATCH --workdir=/home/ulg/mast/acapet/NEMO/azote/

#module purge 
source /home/ulg/mast/acapet/pyload
#module load  EasyBuild Python/3.5.1-foss-2016a

echo 'Running '$1
python $1
