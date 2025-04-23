library(devtools)
library(readr)
library(tidyr)
library(dplyr)
library(scales)

setwd("~/cell_type/Tool/MethScope/")
devtools::load_all()

reference_pattern <- system.file("extdata", "Zhou2025_HumanAtlas.cm", package = "MethScope")
input_pattern <- GenerateInput("/home/fuh1/tmp/20240322_SparseM/Zhou2025_train.cg",reference_pattern)
#input_pattern <- readRDS("/home/fuh1/intermediate/20250419_Zhou2025_training_model_input.rds")
celltype <- read.table("~/cell_type/20250418_eckerhuman.tsv",sep="\t")
sampled_celltype <- celltype$subtype[match(rownames(input_pattern), celltype$cell)]
#saveRDS(input_pattern,"/home/fuh1/intermediate/20250419_Zhou2025_training_model_input.rds")
trained_model <- Input_training(input_pattern,sampled_celltype,cross_validation=F,number_patterns=3000)
saveRDS(trained_model,"/home/fuh1/intermediate/20250419_Zhou2025_training_model.rds")