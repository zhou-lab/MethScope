library(devtools)
library(readr)
library(tidyr)
library(dplyr)
library(scales)

setwd("/home/fuh1/cell_type/Tool/MethScope/analysis")
devtools::load_all()

args <- commandArgs(trailingOnly = TRUE)
reference_pattern <- args[1]
input_pattern_path <- args[2]
cell_type_path <- args[3]
save_path <- args[4]
#print(commandArgs(trailingOnly = TRUE))
input_pattern <- GenerateInput(input_pattern_path,reference_pattern)
celltype <- read.table(cell_type_path,sep="\t",header=T)
# sampled_celltype <- celltype$MajorType[match(rownames(input_pattern), celltype$SampleID)]
sampled_celltype <- celltype$label[match(rownames(input_pattern), celltype$cell)]
trained_model <- Input_training(input_pattern,sampled_celltype,cross_validation=F)
saveRDS(trained_model,save_path)



#reference_pattern <- system.file("extdata", "Zhou2025_HumanAtlas_updated.cm", package = "MethScope")
#input_pattern <- GenerateInput("/home/fuh1/tmp/20240322_SparseM/Zhou2025_train.cg",reference_pattern)
#input_pattern <- readRDS("/home/fuh1/intermediate/20250507_Zhou2025_training_model_input.rds")
# celltype <- read.table("~/cell_type/20250425_eckerhuman.tsv",sep="\t")
# sampled_celltype <- celltype$label[match(rownames(input_pattern), celltype$cell)]
#saveRDS(input_pattern,"/home/fuh1/intermediate/20250419_Zhou2025_training_model_input.rds")
# trained_model <- Input_training(input_pattern,sampled_celltype,cross_validation=T)
#saveRDS(trained_model,"/home/fuh1/intermediate/20250724_Zhou2025_training_model.rds")