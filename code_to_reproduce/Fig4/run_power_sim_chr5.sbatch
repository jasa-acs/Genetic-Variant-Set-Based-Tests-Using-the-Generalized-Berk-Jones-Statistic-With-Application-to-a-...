#!/usr/bin/env bash
#SBATCH -n 1					# Each job runs on 1 core
#SBATCH -N 1
#SBATCH -J power_sim_chr5
#SBATCH -o power_sim_chr5.out
#SBATCH -e power_sim_chr5.err
#SBATCH -p short						# You will need to change this value to the name of the queue used by your system. 
#SBATCH -t 2-00:00					# Jobs should finish in under 2 days
#SBATCH --array=1-320:1				# You must run the job array from 1-320 to capture all the information.
#SBATCH --mem=400


Rscript gen_data_fig4.R ${SLURM_ARRAY_TASK_ID} 4 >power_sim_chr5.${SLURM_ARRAY_TASK_ID}.out 2>power_sim_chr5.${SLURM_ARRAY_TASK_ID}.err


