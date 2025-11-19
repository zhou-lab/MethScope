library(devtools)
library(doParallel)
library(MethScope)

setwd("/home/fuh1/cell_type/Tool/MethScope/analysis")
devtools::load_all()

args <- commandArgs(trailingOnly = TRUE)
reference_pattern <- args[1]
output_dir <- args[2]
model_trained <- args[3]
# Define parameters
#reference_pattern <- system.file("extdata", "Liu2021_MouseBrain.cm", package = "MethScope")
#reference_pattern <- system.file("extdata", "Zhou2025_HumanAtlas.cm", package = "MethScope")
#output_dir <- "/home/fuh1/intermediate/20250725_Zhou2024"
file_list <- list.files(output_dir, pattern = "*.cg$", full.names = TRUE)  # Get all .cg files

num_cores <- 20
cl <- makeCluster(num_cores)
registerDoParallel(cl)

# Parallel loop over split .cg files
results <- foreach(file = file_list, .combine = 'c', .packages = c("MethScope")) %dopar% {
  input_pattern <- GenerateInput(file, reference_pattern)
  #Liu2021_MouseBrain_P1000 <- getFromNamespace("Liu2021_MouseBrain_P1000", "MethScope")
  #prediction_result <- PredictCellType(Liu2021_MouseBrain_P1000, input_pattern, smooth = FALSE)
  #Zhou2025_HumanAtlas_P1000 <- readRDS("/home/fuh1/intermediate/20250724_Zhou2025_training_model.rds")
  #Zhou2025_HumanAtlas_P1000 <- getFromNamespace("Zhou2025_HumanAtlas_P1000", "MethScope")
  model_trained <- readRDS(model_trained)
  prediction_result <- PredictCellType(model_trained, input_pattern, smooth = FALSE)
  # Generate unique output filenames
  file_base <- tools::file_path_sans_ext(basename(file))
  saveRDS(prediction_result, file.path(output_dir, paste0(file_base, "_prediction.rds")))
  saveRDS(input_pattern, file.path(output_dir, paste0(file_base, "_patterns.rds")))
  return(file_base)  # Return the processed file name for logging
}
stopCluster(cl)
cat("Processing completed")
