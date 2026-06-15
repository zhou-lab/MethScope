.load_MethScope_model <- function(model_name) {
  meta <- MethScope_model_metadata[[model_name]]
  if (is.null(meta)) {
    stop("Unknown MethScope model: ", model_name, call. = FALSE)
  }

  booster_path <- system.file("extdata", meta$file, package = "MethScope", mustWork = TRUE)
  con <- gzfile(booster_path, "rb")
  on.exit(close(con), add = TRUE)
  booster_raw <- readBin(con, what = "raw", n = 1e8)
  booster <- xgboost::xgb.load.raw(booster_raw)

  structure(
    list(
      booster = booster,
      cell_type = meta$cell_type,
      npattern = meta$npattern,
      num_class = meta$num_class
    ),
    class = "MethScopeModel"
  )
}

#' Load the Liu et al. mouse brain pretrained model
#'
#' @return A MethScopeModel object for use with \code{\link{PredictCellType}}.
#' @export
Liu2021_MouseBrain_P1000 <- function() {
  .load_MethScope_model("Liu2021_MouseBrain_P1000")
}

#' Load the Zhou et al. human atlas pretrained model
#'
#' @return A MethScopeModel object for use with \code{\link{PredictCellType}}.
#' @export
Zhou2025_HumanAtlas_P1000 <- function() {
  .load_MethScope_model("Zhou2025_HumanAtlas_P1000")
}
