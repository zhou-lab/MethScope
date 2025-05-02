#!/bin/bash
#sbatch --mem=420G
#SBATCH -c 24
#SBATCH -t 30:00:00

start_time=$(date +%s)  # Get start time in seconds

# Define input/output paths
IDX_FILE="/home/fuh1/tmp/20240322_SparseM/Zhou2025_test.cg.idx"
CG_FILE="/home/fuh1/tmp/20240322_SparseM/Zhou2025_test.cg"
OUTPUT_DIR="/home/fuh1/intermediate/20250429_Zhou2024"
N=40  # Number of splits

# Ensure output directory exists
mkdir -p "$OUTPUT_DIR"

# Split the index file into N roughly equal chunks
split -d -n l/$N "$IDX_FILE" "$OUTPUT_DIR/split_"

# Use parallel to process each chunk and generate corresponding .cg file
ls "$OUTPUT_DIR/split_"* | parallel -j "$N" 'yame subset '"$CG_FILE"' -l {} -o '"$OUTPUT_DIR"'/subset_{#}.cg'

echo "Splitting and subset creation complete!"

Rscript-4.5.0.devel ./20250205_run_model.R

end_time=$(date +%s)  # Get end time in seconds
elapsed_time=$((end_time - start_time))

echo "Execution time: $elapsed_time seconds"