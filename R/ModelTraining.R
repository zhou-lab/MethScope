#' Train XGBoost model to predict cell type
#'
#' @param summary_results a wide cell by pattern matrix generated from GenerateInput function
#' @param cell_type_label a vector of the corresponding cell type label for each row of the summary results
#' @param number_patterns a numeric value to indicate number of patterns to be used (Default: 1000)
#' @param cross_validation a boolean varaible whether to perform cross_validation to obtain the best hyper parameters for the model
#' @param xgb_parameters an optional list for xgb model parameters provided by the user
#' @return A MethScopeModel object containing the trained XGBoost booster and metadata.
#' @import xgboost
#' @import caret
#' @import doParallel
#' @import parallel
#' @export

Input_training <- function(summary_results,cell_type_label,number_patterns=1000,
                           cross_validation = FALSE,xgb_parameters = list()){
  # Make sure the input matrix is numeric
  train = do.call(cbind, lapply(summary_results[,1:number_patterns], as.numeric))
  cell_type_label_factor <- as.factor(cell_type_label)
  cell_type_label <- as.numeric(as.factor(cell_type_label)) - 1
  numberOfClasses <- length(unique(cell_type_label))
  # Obtain the XGBoost model paramters
  if(cross_validation == TRUE){
    num_cores <- parallel::detectCores() - 1  # Leave 1 core free for the OS
    # Create and register the parallel backend
    cl <- parallel::makeCluster(num_cores)
    doParallel::registerDoParallel(cl)
    xgb_grid <- expand.grid(
      nrounds = round(sqrt(nrow(train))),
      max_depth = c(4, 6, 8),             # Maximum depth of the trees
      eta = c(0.01, 0.1, 0.3),            # Learning rate
      gamma = c(0, 1),                    # Minimum loss reduction
      colsample_bytree = c(0.5, 0.7, 1),  # Subsample ratio of columns
      min_child_weight = c(1, 3, 5),     # Minimum sum of instance weight
      subsample = c(0.7, 1)         # Subsample ratio of the training instances
    )
    control <- trainControl(method = "cv", number = 5, verboseIter = TRUE)
    train_matrix <- train
    train_labels <- cell_type_label
    xgb_tune <- train(
      x = train_matrix, y = train_labels,
      method = "xgbTree",
      tuneGrid = xgb_grid,
      trControl = control
    )
    parallel::stopCluster(cl)
    best_params <- xgb_tune$bestTune
    xgb_params <- list("objective" = "multi:softprob",
                       "eval_metric" = "mlogloss",
                       "num_class" = numberOfClasses,
                       booster = "gbtree",
                       max_depth = best_params$max_depth,
                       eta = best_params$eta,
                       gamma = best_params$gamma,
                       colsample_bytree = best_params$colsample_bytree,
                       min_child_weight = best_params$min_child_weight,
                       subsample = best_params$subsample)
  } else{
    if(length(xgb_parameters) == 0){
      xgb_params <- list("objective" = "multi:softprob",
                         "eval_metric" = "mlogloss",
                         "num_class" = numberOfClasses,
                         booster = "gbtree")
    } else{xgb_params <- xgb_parameters}
  }

  dtrain <- xgb.DMatrix(data = train, label= cell_type_label)
  train_args <- list(
    params = xgb_params,
    data = dtrain,
    nrounds = round(sqrt(nrow(train))),
    print_every_n = 20
  )
  eval_arg <- if ("evals" %in% names(formals(xgboost::xgb.train))) "evals" else "watchlist"
  train_args[[eval_arg]] <- list(train = dtrain)
  bst_model <- do.call(xgboost::xgb.train, train_args)
  structure(
    list(
      booster = bst_model,
      cell_type = levels(cell_type_label_factor),
      npattern = number_patterns,
      num_class = numberOfClasses
    ),
    class = "MethScopeModel"
  )
}
