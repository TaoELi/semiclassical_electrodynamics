#!/bin/bash

#SBATCH --job-name=TLS_sigma
#SBATCH --time=03:00:00
#SBATCH --mem-per-cpu=2G
#SBATCH --cpus-per-task=1
#SBATCH --partition=taoeli
#SBATCH --array=0-17
#SBATCH --output=array_example_%A_%a.out
#SBATCH --error=array_example_%A_%a.err

i=$SLURM_ARRAY_TASK_ID
echo "Performing jobs $i"
python job_sigmas.py $i  
