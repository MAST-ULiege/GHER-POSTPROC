#!/bin/bash
#PBS -N py_diags
#PBS -q main
#PBS -r y
#PBS -W group_list=bsmfc
#PBS -l walltime=00:20:00
#PBS -l select=1:ncpus=1:mem=15Gb
#PBS -l model=haswell_fat


#pyenv
export PATH="/home/acad/ulg-mast/acapet/.pyenv/bin:$PATH"
eval "$(pyenv init -)"
eval "$(pyenv virtualenv-init -)"
eval "$(pyenv local 2.7.15)"


cd /home/acad/ulg-mast/acapet/NEMO_ANALYSIS

echo 'Running ' $scr
python $scr
