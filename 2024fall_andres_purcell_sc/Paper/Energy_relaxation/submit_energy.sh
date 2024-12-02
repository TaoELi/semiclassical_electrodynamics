#!/bin/bash

#SBATCH --job-name=maxwell
#SBATCH --time=01:00:00
#SBATCH --mem-per-cpu=2G
#SBATCH --cpus-per-task=1
#SBATCH --partition=taoeli
#SBATCH --array=0-7
#SBATCH --output=array_example_%A_%a.out
#SBATCH --error=array_example_%A_%a.err

i=$SLURM_ARRAY_TASK_ID
echo "Performing jobs $i"
python job_energy.py $i  
