#' Predict cell type annotation from the trained model
#'
#' @param bst_model The boosting model trained from ModelTrain
#' @param predictMatrix A wide cell by pattern matrix generated from GenerateInput function
#' @param smooth A Boolean variable to indicate whether smooth the matrix (Default: FALSE)
#' @param KNeighbor number of knn neighbors to use for smoothing (Default: 5)
#' @return A cell by cell type matrix with confidence score and labeled cell type.
#' @import xgboost
#' @importFrom stats setNames
#' @importFrom dplyr mutate
#' @examples
#' # Use the example.cg file included in the package
#' reference_pattern <- system.file("extdata", "Liu2021_MouseBrain.cm", package = "MethScope")
#' example_file <- system.file("extdata", "example.cg", package = "MethScope")
#' result <- GenerateInput(example_file,reference_pattern)
#' prediction_result <- PredictCellType(MethScope:::Liu2021_MouseBrain_P1000,result)
#' @export
PredictCellType <- function(bst_model, predictMatrix,smooth=FALSE,KNeighbor=5) {
  numberOfClasses <- bst_model$params$num_class
  cell_type_factor <- bst_model$cell_type
  number_patterns <- bst_model$npattern
  sample_names <- rownames(predictMatrix)
  predictMatrix = do.call(cbind, lapply(predictMatrix[,1:number_patterns], as.numeric))
  if (smooth){
    Matrix_smooth <- smooth_matrix(predictMatrix,KNeighbor)
    predictMatrix_smooth <- Matrix_smooth[[1]]
    rownames(predictMatrix_smooth) <- rownames(predictMatrix)
    colnames(predictMatrix_smooth) <- colnames(predictMatrix)
    na_loc <- is.na(predictMatrix[,1:number_patterns])
    predictMatrix_smooth[na_loc] <- NA
    predictMatrix <- as.matrix(predictMatrix_smooth)
  }
  dtest <- xgboost::xgb.DMatrix(data = predictMatrix)
  pred_result <- predict(bst_model, newdata = dtest)
  pred_result <- matrix(pred_result, nrow = numberOfClasses,
                            ncol=length(pred_result)/numberOfClasses) %>%
                     t() %>% data.frame() %>%
                     dplyr::mutate(max_prob = max.col(., "last"))
  num_to_factor <- stats::setNames(cell_type_factor, 1:length(cell_type_factor))
  pred_result$prediction_label <- factor(sapply(pred_result$max_prob, function(x) num_to_factor[as.character(x)]), levels = cell_type_factor)
  confiscore <- apply(pred_result[,1:numberOfClasses], 1, confidence_score)
  pred_result$confidence_score <- confiscore
  if(smooth){
  pred_result <- filter_cell(pred_result,Matrix_smooth[[2]])
  }
  rownames(pred_result) <- sample_names
  pred_result
}

#' Produce confidence score for XGBoost prediction
#'
#' @param vec A vector of predicted probability for each cell type
#' @return A numeric confidence score from 0 to 1.
#' @export
confidence_score <- function(vec){
  K <- length(vec)  # Number of classes
  entropy <- -sum(vec * log(vec + 1e-10))  # Compute Shannon entropy
  max_entropy <- log(K)  # Maximum possible entropy
  # Normalize confidence between 0 and 1
  confidence <- 1 - (entropy / max_entropy)
  # Ensure confidence is within the valid range [0,1]
  confidence <- max(0, min(confidence, 1))
  return(confidence)
}

#' Produce confidence score based on top 95 percent for XGBoost prediction
#'
#' @param vec A vector of predicted probability for each cell type
#' @return A numeric confidence score from 0 to 1.
#' @export
confidence_score_top95 <- function(vec){
  percentile_95 <- stats::quantile(vec, 0.95)
  values_above_95th <- vec[vec >= percentile_95]
  max_value <- max(vec)
  metric <- max_value / sum(values_above_95th)
  return(metric)
}

#' Smooth cell by pattern matrix to reduce noise
#'
#' @param predictMatrix A wide cell by pattern matrix generated from GenerateInput function
#' @param KNeighbor Number of knn neighbors to use for smoothing (Default: 5)
#' @importFrom FNN get.knn
#' @return A wide cell by pattern matrix after smoothing and knn graph
#' @export
#' 
smooth_matrix <- function(predictMatrix,KNeighbor = 5){
  all_na_cols <- apply(predictMatrix, 2, function(x) all(is.na(x)))
  predictMatrix[, all_na_cols] <- 1
  k <- which(is.na(predictMatrix), arr.ind=TRUE)
  predictMatrix[k] <- colMeans(predictMatrix, na.rm=TRUE)[k[,2]]
  knn_res <- FNN::get.knn(predictMatrix, k = KNeighbor)
  smooth_methylation <- function(data, knn_res) {
    smoothed_matrix <- matrix(0, nrow = nrow(data), ncol = ncol(data))
    for (i in 1:nrow(data)) {
      neighbors <- knn_res$nn.index[i, ]  # Get indices of k nearest neighbors
      smoothed_matrix[i, ] <- colMeans(data[neighbors, , drop = FALSE])  # Compute mean across neighbors
    }
    return(smoothed_matrix)
  }
  predictMatrix <- smooth_methylation(predictMatrix, knn_res)
  predictMatrix <- as.data.frame(predictMatrix)
  list(predictMatrix,knn_res)
}

#' Filter final prediction to reduce noise
#'
#' @param pred_result The prediction result from XGBoost
#' @param knn_res knn graph from smooth_matrix
#' @param KNeighbor Number of knn neighbors to use for smoothing (Default: 5)
#' @return The final prediction result after dropping few cell types
#' @export
#' 
filter_cell <- function(pred_result,knn_res,KNeighbor = 5){
  cell_type_counts <- table(pred_result$prediction_label)
  low_confidence_idx <- which(cell_type_counts[pred_result$prediction_label] < KNeighbor & pred_result$confidence_score < 0.5)
  majority_vote <- function(neighbors) {
    labels <- pred_result$prediction_label[neighbors]  # Get labels of nearest neighbors
    labels <- labels[!is.na(labels)]  # Remove NAs (if any)
    return(names(sort(table(labels), decreasing = TRUE)[1]))  # Return most frequent label
  }
  for (i in low_confidence_idx) {
    pred_result$prediction_label[i] <- majority_vote(knn_res$nn.index[i,])
  }
  pred_result
}


#' Estimate cell type relative proportion
#'
#' @param ref An imputed wide cell by pattern matrix generated from GenerateInput function using reference Pseudobulk
#' @param mixture_matrix An imputed wide cell by pattern matrix generated from GenerateInput function
#' @param number_patterns a numeric value to indicate number of patterns to be used (Default: 1000)
#' @param var_threshold a numeric value to indicate variance that should filter the patterns (Default: 0.1)
#' @return A cell type by cell matrix showing the relative cell type proportion estimate for each cells 
#' @import nnls
#' @export
#' 
nnls_deconv <- function(ref, mixture_matrix,number_patterns= 1000,var_threshold=0.01) {
  ref <- t(ref[,1:number_patterns])
  mixture_matrix <- t(mixture_matrix[,1:number_patterns])
  common_rows <- intersect(rownames(ref), rownames(mixture_matrix))
  ref <- ref[common_rows, , drop = FALSE]
  mixture_matrix <- mixture_matrix[common_rows, , drop = FALSE]
  mixture_matrix <- as.matrix(mixture_matrix[rownames(ref),])
  row_vars_ref <- apply(ref, 1, var)
  high_var_rows <- names(row_vars_ref[row_vars_ref > var_threshold])
  ref <- ref[high_var_rows, , drop = FALSE]
  mixture_matrix <- mixture_matrix[high_var_rows, , drop = FALSE]
  result <- apply(mixture_matrix, 2, function(sample) {
    fit <- nnls(ref, sample)
    prop <- fit$x
    prop / sum(prop)
  })
  colnames(result) <- colnames(mixture_matrix)
  rownames(result) <- colnames(ref)
  return(result)
}
