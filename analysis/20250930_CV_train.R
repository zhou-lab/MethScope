library(devtools)
library(readr)
library(tidyr)
library(dplyr)
library(scales)
library(readxl)
library(data.table)
setwd("~/cell_type/Tool/MethScope/")
devtools::load_all()

celltype <- read.table("~/cell_type/20250425_eckerhuman.tsv",sep="\t")
#celltype <- read_xlsx("/home/fuh1/cell_type/20220630_Liu2021_mouseBrainWGBS.xlsx")
set.seed(20250930)
celltype_shuffled <- celltype[sample(nrow(celltype)), ]
K <- 5
fold_ids <- sample(rep(1:K, length.out = nrow(celltype_shuffled)))

setwd("/home/fuh1/intermediate/MethScope/analysis/")
for (k in 1:K) {
  test_idx <- which(fold_ids == k)
  train_idx <- which(fold_ids != k)
  test  <- celltype_shuffled[test_idx, ]
  train <- celltype_shuffled[train_idx, ]
  output_train_sample <- train$cell
  output_test_sample <- test$cell
  write.table(output_train_sample,paste0("./train_test_sample/Zhou2025_training_sample_",k,".txt"),quote=FALSE,col.names = FALSE,row.names = FALSE)
  write.table(output_test_sample,paste0("./train_test_sample/Zhou2025_testing_sample_",k,".txt"),quote=FALSE,col.names = FALSE,row.names = FALSE)
  fwrite(train, file = paste0("./train_test_sample/Zhou2025_training_sample_",k,".tsv"), sep = "\t")
  fwrite(test, file = paste0("./train_test_sample/Zhou2025_testing_sample_",k,".tsv"), sep = "\t")
}

