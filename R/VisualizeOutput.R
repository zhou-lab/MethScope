#' Generate UMAP for the final prediction based on cell patterns
#'
#' @param predictMatrix a wide cell by pattern matrix generated from GenerateInput function
#' @param prediction_result Prediction result from PredictCellType
#' @param n_component Number of PCA components to use (Default: 30)
#' @param seed A number for random seed (Default: 123)
#' @param ... Additional arguments passed to `uwot::umap` (e.g., `n_neighbors`, `metric`).
#' @return A list of two ggplot2 UMAP object.
#' @useDynLib MethScope, .registration = TRUE
#' @import ggplot2
#' @import uwot
#' @importFrom stringr str_extract
#' @importFrom tidyr spread
#' @importFrom utils read.table
#' @importFrom stats prcomp
#' @export
PlotUMAP <- function(predictMatrix,prediction_result,n_component=30,seed=123,...) {
  set.seed(seed)
  summary_results <- predictMatrix
  summary_results$cell_type <- prediction_result$prediction_label
  summary_results_value <- summary_results %>% dplyr::select(-"cell_type")
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
  df_plot$conf <- prediction_result$confidence_score
  colnames(df_plot) <- c("UMAP1","UMAP2","color","conf")
  plot1 <- ggplot(df_plot) + geom_point(aes(UMAP1,UMAP2,fill=color),pch = 21,size = 1, stroke=NA) +
    theme_bw()+labs(fill="")+theme(panel.grid.major = element_blank(),
                                   panel.grid.minor = element_blank(),
                                   panel.background = element_blank())+
    guides(fill = guide_legend(override.aes = list(size = 4)))
  plot2 <- ggplot(df_plot) + geom_point(aes(UMAP1,UMAP2,fill=conf),pch = 21,size = 1, stroke=NA) +
    theme_bw()+labs(fill="")+theme(panel.grid.major = element_blank(),
                                   panel.grid.minor = element_blank(),
                                   panel.background = element_blank())
  list(plot1,plot2)
}

#' Generate UMAP for the final prediction based on fixed window eg.100kb bin widows
#'
#' @param query_fn File path to query .cg
#' @param knowledge_fn File path to 100bk bins window or reference pattern 
#' @param prediction_result Prediction result from PredictCellType
#' @param n_component Number of PCA components to use (Default: 30)
#' @param seed A number for random seed (Default: 123)
#' @param ... Additional arguments passed to `uwot::umap` (e.g., `n_neighbors`, `metric`).
#' @return A list of two ggplot2 UMAP object.
#' @useDynLib MethScope, .registration = TRUE
#' @import ggplot2
#' @import uwot
#' @importFrom stringr str_extract
#' @importFrom tidyr spread
#' @importFrom utils read.table
#' @importFrom stats prcomp
#' @export
PlotUMAP_fixedwindow <- function(query_fn, knowledge_fn,prediction_result,n_component=30,seed=123,...) {
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
  df_plot$conf <- prediction_result$confidence_score
  colnames(df_plot) <- c("UMAP1","UMAP2","color","conf")
  plot1 <- ggplot(df_plot) + geom_point(aes(UMAP1,UMAP2,fill=color),pch = 21,size = 1, stroke=NA) +
                          theme_bw()+labs(fill="")+theme(panel.grid.major = element_blank(),
                                                         panel.grid.minor = element_blank(),
                                                         panel.background = element_blank())
  plot2 <- ggplot(df_plot) + geom_point(aes(UMAP1,UMAP2,fill=conf),pch = 21,size = 1, stroke=NA) +
                          theme_bw()+labs(fill="")+theme(panel.grid.major = element_blank(),
                                   panel.grid.minor = element_blank(),
                                   panel.background = element_blank())
  list(plot1,plot2)
}

#' Generate confusion table for the final prediction 
#'
#' @param prediction_result Prediction result from PredictCellType
#' @param actual_label Ground truth cell label
#' @param log2 Log scale count (Default: False)
#' @return A ggplot2 confusion table object.
#' @import ggplot2
#' @importFrom caret confusionMatrix
#' @export
PlotConfusion <- function(prediction_result,actual_label,log2=F) {
  confusion_matrix <- caret::confusionMatrix(factor(prediction_result$prediction_label),
                                      factor(actual_label),mode = "everything")
  cm_table <- as.data.frame(confusion_matrix$table)
  
  # Compute accuracy
  accuracy_rate <- sum(diag(confusion_matrix$table)) / sum(confusion_matrix$table)
  
  # Rename columns for ggplot
  colnames(cm_table) <- c("Actual", "Predicted", "Freq")
  legend_title <- "Count"
  if(log2){cm_table$Freq <- log2(cm_table$Freq+1)
  legend_title <- "Log2(Count)"}
  
  # Plot using ggplot2
  p1 <- ggplot(cm_table, aes(x = Predicted, y = Actual, fill = Freq)) +
    geom_tile(color = "white") +
    scale_fill_distiller(direction=1)+
    labs(title = paste0("Confusion Matrix (Accuracy: ", round(accuracy_rate, 3), ")"),
         x = "Predicted Label", y = "Actual Label", fill = legend_title) +
    theme_minimal() +
    theme(plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
          axis.text.x = element_text(size = 8, angle = 90, vjust = 0.5, hjust = 1),  # Rotated and smaller
          axis.text.y = element_text(size = 8),legend.position = "right")
  p1
}

#' Generate F1 score barplot for each class
#'
#' @param prediction_result Prediction result from PredictCellType
#' @param actual_label Ground truth cell label
#' @return A ggplot2 object.
#' @import ggplot2
#' @importFrom caret confusionMatrix
#' @importFrom stringr str_remove
#' @export
PlotF1 <- function(prediction_result,actual_label) {
  confusion_matrix <- caret::confusionMatrix(factor(prediction_result$prediction_label),
                                             factor(actual_label),mode = "everything")
  f1_scores <- confusion_matrix$byClass[, "F1"]
  f1_df <- data.frame(Class = rownames(confusion_matrix$byClass),
                      F1_Score = f1_scores)
  f1_df$F1_Score[is.na(f1_df$F1_Score)] <- 0
  f1_df$Class <- stringr::str_remove(f1_df$Class, "^Class: ")
  ggplot(f1_df, aes(x = reorder(Class, F1_Score), y = F1_Score, fill = F1_Score)) +
    geom_bar(stat = "identity") +
    coord_flip() +  # Flip to horizontal for better readability
    scale_fill_distiller(palette = "YlOrBr",direction=1) +  # Better contrast
    labs(title = "",
         x = "",
         y = "F1-score",
         fill = "") +
    theme_bw() +
    theme(plot.title = element_text(hjust = 0.5, face = "bold"),
          axis.text.y = element_text(size = 10),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.background = element_blank(),
          legend.position = "None") 
}

#' Impute missing value for 100K window matrix
#'
#' @param mtx A cell by 100K window data frame with missing values
#' @param na_percent A na percent threshold to be filterd (Default: 30)
#' @return A cell by 100K window data frame with imputed values
#' @importFrom stats quantile
#' @export
#'
imputeRowMean <- function(mtx,na_percent = 30) {
  na_percentage <- sapply(mtx, function(x) sum(is.na(x)) / length(x)) * 100
  mtx <- mtx[, na_percentage < na_percent]
  k <- which(is.na(mtx), arr.ind=TRUE)
  mtx[k] <- rowMeans(mtx, na.rm=TRUE)[k[,1]]
  mtx
}
