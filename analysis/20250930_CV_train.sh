#!/bin/bash
#SBATCH --job-name=trainCV
#SBATCH --output=logs/trainCV_%j.out
#SBATCH --error=logs/trainCV_%j.err
#SBATCH --mem=420G
#SBATCH -c 24
#SBATCH -t 120:00:00

#mkdir -p ./train_test_sample/

start_time=$(date +%s) 
# raw_cg="/mnt/isilon/zhou_lab/projects/20220324_SingleCellMeth/mm10/2021_Liu.cg"
raw_cg="/mnt/isilon/zhou_lab/projects/20220324_SingleCellMeth/hg38/2025_Zhou.cg"
#Rscript-4.5.0.devel ./20250930_CV_train.R
#for k in {1..5}; do
for k in {4..5}; do
  sample_train="/home/fuh1/intermediate/MethScope/analysis/train_test_sample/Zhou2025_training_sample_${k}.txt"
  sample_test="/home/fuh1/intermediate/MethScope/analysis/train_test_sample/Zhou2025_testing_sample_${k}.txt"
  tsv="/home/fuh1/intermediate/MethScope/analysis/train_test_sample/Zhou2025_training_sample_${k}.tsv"
  cg="/home/fuh1/intermediate/MethScope/analysis/train_test_sample/Zhou2025_training_sample_${k}.cg"
  IDX_FILE="/home/fuh1/intermediate/MethScope/analysis/train_test_sample/Zhou2025_testing_sample_${k}.cg.idx"
  CG_FILE="/home/fuh1/intermediate/MethScope/analysis/train_test_sample/Zhou2025_testing_sample_${k}.cg"
  yame subset "$raw_cg" -l "$sample_train" -o "$cg"
  yame subset "$raw_cg" -l "$sample_test" -o "$CG_FILE"

  echo "[INFO] Running GenerateReference.sh on fold $k"
  ./GenerateReference.sh "$tsv" "$cg"

  reference_pattern="/home/fuh1/cell_type/Tool/MethScope/analysis/patterns.cm"
  save_path="/home/fuh1/intermediate/MethScope/analysis/models_Zhou2025/trained_model_fold${k}.rds"
  mkdir -p "/home/fuh1/intermediate/MethScope/analysis/models_Zhou2025/"
  echo "[INFO] Training Fold $k"
  Rscript-4.5.0.devel ./20250418_train_model.R "$reference_pattern" "$cg" "$tsv" "$save_path"

  OUTPUT_DIR="/home/fuh1/intermediate/MethScope/analysis/models_output_Zhou2025/${k}/"
  mkdir -p "$OUTPUT_DIR"
  N=10
  split -d -n l/$N "$IDX_FILE" "$OUTPUT_DIR/split_"
  ls "$OUTPUT_DIR/split_"* | parallel -j "$N" 'yame subset '"$CG_FILE"' -l {} -o '"$OUTPUT_DIR"'/subset_{#}.cg'
  echo "[INFO] Splitting complete Fold $k"

  echo "[INFO] Modeling Fold $k"
  Rscript-4.5.0.devel ./20250205_run_model.R \ "$reference_pattern" \
    "$OUTPUT_DIR" \ "$save_path"
done

end_time=$(date +%s)
elapsed_time=$((end_time - start_time))
echo "Execution time: $elapsed_time seconds"
