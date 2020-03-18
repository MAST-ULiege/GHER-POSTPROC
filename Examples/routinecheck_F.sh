#!/bin/bash 

# qsub -v scr=OSR5_F.py,y=local_NEMO_OSR5d.yml                       -N OSR5   -l walltime=03:00:00 -l select=1:ncpus=1:mem=20Gb fat_submit_zenobe_args.sh      
#qsub -v scr=TimeSeriesComp_R_F.py,y=local_NEMO_OSR5d.yml           -N TSC    -l walltime=01:30:00 -l select=1:ncpus=1:mem=1Gb submit_zenobe_args.sh  
#qsub -v scr=TimeSeriesComp_R_Benth_F.py,y=local_NEMO_OSR5d.yml     -N TSC_B  -l walltime=02:00:00 -l select=1:ncpus=1:mem=1Gb submit_zenobe_args.sh  
qsub -v scr=SurfInsight2_OSR5_F.py,y=local_NEMO_OSR5d.yml          -N SI     -l walltime=06:00:00 -l select=1:ncpus=1:mem=40Gb fat_submit_zenobe_args.sh      
qsub -v scr=BenthInsight_OSR5_F.py,y=local_NEMO_OSR5d.yml          -N BI     -l walltime=04:00:00 -l select=1:ncpus=1:mem=50Gb fat_submit_zenobe_args.sh      
qsub -v scr=gNTRACX_OSR5_F.py,y=local_NEMO_OSR5d.yml               -N NT     -l walltime=06:00:00 -l select=1:ncpus=1:mem=40Gb fat_submit_zenobe_args.sh  

