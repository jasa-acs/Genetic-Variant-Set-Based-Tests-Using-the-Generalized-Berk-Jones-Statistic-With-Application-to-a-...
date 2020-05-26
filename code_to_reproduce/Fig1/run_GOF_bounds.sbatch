#!/usr/bin/env bash
#SBATCH -n 1					# Each job runs on 1 core
#SBATCH -N 1
#SBATCH -J GOF_bounds.R
#SBATCH -o GOF_bounds.out
#SBATCH -e GOF_bounds.err
#SBATCH -p short						# You will need to change this value to the name of the queue used by your system. 
#SBATCH -t 1-00:00					# Jobs should finish in under 12 hours
#SBATCH --array=3-21:1				# This will produce all the necessary data for Fig 1, plus some extra data
#SBATCH --mem=5000


Rscript gen_data_fig1.R ${SLURM_ARRAY_TASK_ID} 1 >GOF_bounds.${SLURM_ARRAY_TASK_ID}.out 2>GOF_bounds.${SLURM_ARRAY_TASK_ID}.err


