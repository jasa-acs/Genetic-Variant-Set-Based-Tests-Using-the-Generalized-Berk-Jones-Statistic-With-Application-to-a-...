#!/usr/bin/env bash
#SBATCH -n 1					# Each job runs on 1 core
#SBATCH -N 1
#SBATCH -J supp_fig1_part2.R
#SBATCH -o supp_fig1_part2.out
#SBATCH -e supp_fig1_part2.err
#SBATCH -p short						# You will need to change this value to the name of the queue used by your system. 
#SBATCH -t 1-00:00					# Jobs should finish in under 12 hours
#SBATCH --array=1-600:1				
#SBATCH --mem=500


Rscript gen_data_suppfig1.R ${SLURM_ARRAY_TASK_ID} 2 >supp_fig1_part2.${SLURM_ARRAY_TASK_ID}.out 2>supp_fig1_part2.${SLURM_ARRAY_TASK_ID}.err


