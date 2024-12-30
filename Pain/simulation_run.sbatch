#!/bin/bash
#SBATCH --job-name=sim_8shape_smo    # Name for your job
#SBATCH --array=1-10             # Create 10 parallel jobs (array indices 1 to 10)
#SBATCH --ntasks=1               # Number of tasks per job
#SBATCH --mem=64000              # Memory per job (64 GB shared across all tasks)
#SBATCH --output=sim_8shape_smo_%A_%a.out  # Output file (JobID and array index in filename)
#SBATCH --error=sim_8shape_smo_%A_%a.err   # Error file (JobID and array index in filename)
#SBATCH --mail-user=jchoi177@jhu.edu   # Email for notifications
#SBATCH --mail-type=END,FAIL     # Notify only for job end or failure

module load R

Rscript --vanilla sim_8shape/sim_8shape_smo.R ${SLURM_ARRAY_TASK_ID}
