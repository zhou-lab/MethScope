#' Generate UMAP for the final prediction based on 100kb bin widows
#'
#' @param query_fn File path to query .cg
#' @param knowledge_fn File path to 100bk bins window
#' @param prediction_result Prediction result from PredictCellType
#' @param n_component Number of PCA components to use (Default: 30)
#' @param seed A number for random seed (Default: 123)
#' @param ... Additional arguments passed to `uwot::umap` (e.g., `n_neighbors`, `metric`).
#' @return A ggplot2 UMAP object.
#' @useDynLib MethScope, .registration = TRUE
#' @import ggplot2
#' @import uwot
#' @importFrom stringr str_extract
#' @importFrom tidyr spread
#' @importFrom utils read.table
#' @importFrom stats prcomp
#' @export
PlotUMAP <- function(query_fn, knowledge_fn,prediction_result,n_component=30,seed=123,...) {
  set.seed(seed)
  stopifnot(is.character(query_fn), is.character(knowledge_fn))
  if (.Platform$OS.type == "windows") {
    stop("Testing sequencing data does not support Windows.")
  }
  yame_result <- .Call("yame_summary_cfunc", query_fn, knowledge_fn)

  summary_results <- utils::read.table(text = paste(yame_result, collapse = "\n"), header = TRUE)

  summary_results <- summary_results %>%
    dplyr::select('Query','Mask','Beta') %>%
    tidyr::spread(key=Mask,value=Beta,fill = NA)

  summary_results$cell_type <- prediction_result$prediction_label
  summary_results_value <- summary_results %>% dplyr::select(-"Query",-"cell_type")
  summary_results_value <- do.call(cbind, lapply(summary_results_value[,1:ncol(summary_results_value)], as.numeric))
  summary_results_value <- imputeRowMean(as.data.frame(summary_results_value))

  pc <- stats::prcomp(as.data.frame(summary_results_value),scale=TRUE)
  pc_data <- pc$x
  UMAP_results = uwot::umap(pc_data[,1:n_component], n_neighbors = 10L, n_components = 2L, metric = "cosine",
                      n_epochs = NULL,learning_rate = 1, min_dist = 0.3,spread = 1,set_op_mix_ratio = 1,
                      local_connectivity = 1L,...)
  UMAP_results = as.data.frame(UMAP_results)
  colnames(UMAP_results) = c("UMAP1", "UMAP2")
  df_plot <- cbind(UMAP_results$UMAP1,UMAP_results$UMAP2)
  df_plot <- as.data.frame(df_plot)
  df_plot$col <- summary_results$cell_type
  df_plot$col <- as.factor(df_plot$col)
  colnames(df_plot) <- c("UMAP1","UMAP2","color")
  plot <- ggplot(df_plot) + geom_point(aes(UMAP1,UMAP2,fill=color),pch = 21,size = 1, stroke=NA) +
                          theme_bw()+labs(fill="")+theme(panel.grid.major = element_blank(),
                                                         panel.grid.minor = element_blank(),
                                                         panel.background = element_blank())
  plot
}

#' Impute missing value for 100K window matrix
#'
#' @param mtx A cell by 100K window data frame with missing values
#' @return A cell by 100K window data frame with imputed values
#' @importFrom stats quantile
#' @export
#'
imputeRowMean <- function(mtx) {
  na_percentage <- sapply(mtx, function(x) sum(is.na(x)) / length(x)) * 100
  mtx <- mtx[, na_percentage < 30]
  k <- which(is.na(mtx), arr.ind=TRUE)
  mtx[k] <- rowMeans(mtx, na.rm=TRUE)[k[,1]]
  mtx
}
