#!/bin/bash
#sbatch --mem=200G
#SBATCH -c 25
#SBATCH -t 30:00:00

start_time=$(date +%s)  # Get start time in seconds

Rscript ./20250205_run_model.R

end_time=$(date +%s)  # Get end time in seconds
elapsed_time=$((end_time - start_time))

echo "Execution time: $elapsed_time seconds"