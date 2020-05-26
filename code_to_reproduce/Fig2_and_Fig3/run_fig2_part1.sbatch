#!/usr/bin/env bash
#SBATCH -n 1					# Each job runs on 1 core
#SBATCH -N 1
#SBATCH -J fig2_part1
#SBATCH -o fig2_part1.out
#SBATCH -e fig2_part1.err
#SBATCH -p short						# You will need to change this value to the name of the queue used by your system. 
#SBATCH -t 2-00:00					# Jobs should finish in under 2 days
#SBATCH --array=1-50:1				# You must run the job array from 1-50 to capture all the information.
#SBATCH --mem=200


Rscript gen_sim_data_fig2.R ${SLURM_ARRAY_TASK_ID} 1 >fig2_part1.${SLURM_ARRAY_TASK_ID}.out 2>fig2_part1.${SLURM_ARRAY_TASK_ID}.err


