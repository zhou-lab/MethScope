test_that("pretrained models load and predict with the current xgboost", {
  model_loaders <- list(
    Liu2021_MouseBrain_P1000,
    Zhou2025_HumanAtlas_P1000
  )

  for (loader in model_loaders) {
    model <- loader()
    expect_s3_class(model, "MethScopeModel")
    expect_true(any(inherits(model$booster, c("xgb.Booster", "xgb.Booster.handle"))))
    expect_equal(model$npattern, 1000)
    expect_length(model$cell_type, model$num_class)

    dummy <- as.data.frame(matrix(0, nrow = 2, ncol = model$npattern))
    rownames(dummy) <- c("cell1", "cell2")
    colnames(dummy) <- paste0("P", seq_len(model$npattern))

    pred <- PredictCellType(model, dummy)
    expect_equal(nrow(pred), 2)
    expect_true("prediction_label" %in% colnames(pred))
    expect_true("confidence_score" %in% colnames(pred))
    expect_equal(
      unname(as.numeric(pred[1, seq_len(model$num_class)])),
      unname(as.numeric(pred[2, seq_len(model$num_class)]))
    )
  }
})

test_that("Input_training models work with PredictCellType", {
  set.seed(1)
  input <- as.data.frame(matrix(runif(40), nrow = 8, ncol = 5))
  colnames(input) <- paste0("P", seq_len(5))
  rownames(input) <- paste0("cell", seq_len(8))
  labels <- rep(c("A", "B"), each = 4)

  model <- Input_training(input, labels, number_patterns = 5)
  expect_s3_class(model, "MethScopeModel")
  expect_true(any(inherits(model$booster, c("xgb.Booster", "xgb.Booster.handle"))))
  expect_equal(model$npattern, 5)
  expect_equal(model$cell_type, c("A", "B"))

  pred <- PredictCellType(model, input)
  expect_equal(nrow(pred), nrow(input))
  expect_true("prediction_label" %in% colnames(pred))
  expect_true("confidence_score" %in% colnames(pred))
})
