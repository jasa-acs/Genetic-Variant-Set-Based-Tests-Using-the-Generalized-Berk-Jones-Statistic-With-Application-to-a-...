#!/usr/bin/env bash
#SBATCH -n 1					# Each job runs on 1 core
#SBATCH -N 1
#SBATCH -J fig2_part6
#SBATCH -o fig2_part6.out
#SBATCH -e fig2_part6.err
#SBATCH -p short						# You will need to change this value to the name of the queue used by your system. 
#SBATCH -t 2-00:00					# Jobs should finish in under 2 days
#SBATCH --array=1-10:1				# You must run the job array from 1-10 to capture all the information.
#SBATCH --mem=200


Rscript gen_omni_data_fig2.R ${SLURM_ARRAY_TASK_ID} 2 >fig2_part6.${SLURM_ARRAY_TASK_ID}.out 2>fig2_part6.${SLURM_ARRAY_TASK_ID}.err


