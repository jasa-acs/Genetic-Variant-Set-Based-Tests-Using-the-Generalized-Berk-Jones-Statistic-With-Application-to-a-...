#!/usr/bin/env bash
#SBATCH -n 1					# Each job runs on 1 core
#SBATCH -N 1
#SBATCH -J tab1_part2
#SBATCH -o tab1_part2.out
#SBATCH -e tab1_part2.err
#SBATCH -p short						# You will need to change this value to the name of the queue used by your system. 
#SBATCH -t 8-00:00					# Jobs should finish in under 8 days
#SBATCH --array=1-2400:1				# You must run the job array from 1-2400 to capture all the information.
#SBATCH --mem=500


Rscript gen_data_tab1_part2.R ${SLURM_ARRAY_TASK_ID} 1 >tab1_part2.${SLURM_ARRAY_TASK_ID}.out 2>tab1_part2.${SLURM_ARRAY_TASK_ID}.err


