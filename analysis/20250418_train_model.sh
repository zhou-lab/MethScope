#!/bin/bash
#SBATCH --mem=400G
#SBATCH -c 6
#SBATCH -t 30:00:00

start_time=$(date +%s)  # Get start time in seconds

Rscript-4.5.0.devel ./20250418_train_model.R

end_time=$(date +%s)  # Get end time in seconds
elapsed_time=$((end_time - start_time))

echo "Execution time: $elapsed_time seconds"