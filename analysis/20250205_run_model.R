library(devtools)
setwd("/home/fuh1/cell_type/Tool/MethScope/")
devtools::load_all()

reference_pattern <- system.file("extdata", "Liu2021_MouseBrain.cm", package = "MethScope")
input_pattern <- GenerateInput("/home/fuh1/intermediate/20250130_Liu2021.cg",reference_pattern)
prediction_result <- PredictCellType(Liu2021_MouseBrain_P1000,input_pattern,smooth = F)
saveRDS(prediction_result,"/home/fuh1//intermediate/20250130_Liu2021_prediction.rds")
saveRDS(input_pattern,"/home/fuh1//intermediate/20250130_Liu2021_prediction_patterns.rds")