#!/usr/bin/env bash
#SBATCH -n 1					# Each job runs on 1 core
#SBATCH -N 1
#SBATCH -J suppfig4
#SBATCH -o suppfig4.out
#SBATCH -e suppfig4.err
#SBATCH -p short						# You will need to change this value to the name of the queue used by your system. 
#SBATCH -t 2-00:00					# Jobs should finish in under 2 days
#SBATCH --array=1-40:1				# You must run the job array from 1-40 to capture all the information.
#SBATCH --mem=400


Rscript gen_data_suppfig456.R ${SLURM_ARRAY_TASK_ID} 1 >suppfig4.${SLURM_ARRAY_TASK_ID}.out 2>suppfig4.${SLURM_ARRAY_TASK_ID}.err


