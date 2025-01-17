#' Predict cell type annotation from the trained model
#'
#' @param bst_model boosting model trained from ModelTrain
#' @param predictMatrix a wide cell by pattern matrix generated from GenerateInput function
#' @return A cell by pattern matrix.
#' @import xgboost
#' @importFrom stats setNames
#' @export
PredictCellType <- function(bst_model, predictMatrix) {
  numberOfClasses <- bst_model$params$num_class
  cell_type_factor <- bst_model$cell_type
  number_patterns <- bst_model$npattern
  predictMatrix = do.call(cbind, lapply(predictMatrix[,1:number_patterns], as.numeric))
  dtest <- xgboost::xgb.DMatrix(data = predictMatrix)
  pred_result <- predict(bst_model, newdata = dtest)
  pred_result <- matrix(pred_result, nrow = numberOfClasses,
                            ncol=length(pred_result)/numberOfClasses) %>%
                     t() %>% data.frame() %>%
                     mutate(max_prob = max.col(., "last"))
  num_to_factor <- stats::setNames(cell_type_factor, 1:length(cell_type_factor))
  pred_result$prediction_label <- factor(sapply(pred_result$max_prob, function(x) num_to_factor[as.character(x)]), levels = cell_type_factor)
  confiscore <- apply(pred_result[,1:numberOfClasses], 1, confidence_score)
  pred_result$confidence_score <- confiscore
  pred_result
}

#' Produce confidence score for XGBoost prediction
#'
#' @param vec A vector of predicted probability for each cell type
#' @return A numeric confidence score from 0 to 1.
#' @importFrom stats quantile
#' @export
confidence_score <- function(vec){
  percentile_95 <- stats::quantile(vec, 0.95)
  values_above_95th <- vec[vec >= percentile_95]
  max_value <- max(vec)
  metric <- max_value / sum(values_above_95th)
  return(metric)
}


