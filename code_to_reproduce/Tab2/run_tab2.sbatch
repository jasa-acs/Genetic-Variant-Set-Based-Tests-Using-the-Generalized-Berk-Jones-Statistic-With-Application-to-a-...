#!/usr/bin/env bash
#SBATCH -n 1					# Each job runs on 1 core
#SBATCH -N 1
#SBATCH -J tab2
#SBATCH -o tab2.out
#SBATCH -e tab2.err
#SBATCH -p short						# You will need to change this value to the name of the queue used by your system. 
#SBATCH -t 5-00:00					# Jobs should finish in under 5 days
#SBATCH --array=1-181:1				# You must run the job array from 1-181 to capture all the information.
#SBATCH --mem=2000


Rscript gen_data_tab2.R ${SLURM_ARRAY_TASK_ID} 1 >tab2.${SLURM_ARRAY_TASK_ID}.out 2>tab2.${SLURM_ARRAY_TASK_ID}.err


