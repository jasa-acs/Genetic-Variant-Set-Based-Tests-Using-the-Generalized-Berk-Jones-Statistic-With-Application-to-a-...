#!/usr/bin/env bash
#SBATCH -n 1					# Each job runs on 1 core
#SBATCH -N 1
#SBATCH -J supp_tab4_part4
#SBATCH -o supp_tab4_part4.out
#SBATCH -e supp_tab4_part4.err
#SBATCH -p short						# You will need to change this value to the name of the queue used by your system. 
#SBATCH -t 0-12:00					# Jobs should finish in under 12 hours
#SBATCH --array=1-50:1				# You must run the job array from 1-50 to capture all the information.
#SBATCH --mem=200


Rscript gen_sim_data_supptab4.R ${SLURM_ARRAY_TASK_ID} 8 >supp_tab4_part4.${SLURM_ARRAY_TASK_ID}.out 2>supp_tab4_part4.${SLURM_ARRAY_TASK_ID}.err


