#!/bin/bash
#SBATCH --job-name sim_8shape	# Name for your job
#SBATCH --ntasks 1		# Number of task
#SBATCH --time 5		# Runtime in minutes.
#SBATCH --mem=2000		# Reserve 2 GB RAM for the job
#SBATCH --partition small 	# Partition to submit
#SBATCH --output sim_8shape.out	# Standard out goes to this file
#SBATCH --error sim_8shape.err	# Standard err goes to this file
#SBATCH --mail-user jchoi177@jhu.edu	# this is the email you wish to be notified at
#SBATCH --mail-type ALL		# ALL will alert you of job beginning, completion, failure etc

module load R

Rscript beta_infer_setting.R # Execute the script