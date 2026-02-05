#!/usr/bin/env bash
#SBATCH --job-name=scTransientEval
#SBATCH --output=logs/%x_%A_%a.out
#SBATCH --error=logs/%x_%A_%a.err
#SBATCH --array=0-14999
#SBATCH --time=0:30:00
#SBATCH --partition=defq
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=16G

# Create log directory if not exists
mkdir -p logs

# Load required modules
## NOTE: This will be different on your HPC; find the appropriate Singularity module.
module load singularity-apptainer/1.1.6

# Echo basic context
echo "Job ID: ${SLURM_JOB_ID}"
echo "Array Task ID: ${SLURM_ARRAY_TASK_ID}"
echo "Host: $(hostname)"
echo "Started at: $(date)"


singularity exec method_evaluations/method_eval.sif python3.11 method_evaluations/method_execution.py --index ${SLURM_ARRAY_TASK_ID} --output_dir results --parameter_file method_evaluations/params.py

echo "Finished at: $(date)"
