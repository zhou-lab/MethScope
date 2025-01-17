#' Generate pattern level data for cell type annotation
#'
#' @param query_fn File path to query
#' @param knowledge_fn File path to pattern file
#' @return A cell by pattern matrix.
#' @useDynLib MethScope, .registration = TRUE
#' @import dplyr
#' @importFrom stringr str_extract
#' @importFrom tidyr spread
#' @importFrom utils read.table
#' @export
GenerateInput <- function(query_fn, knowledge_fn) {

  stopifnot(is.character(query_fn), is.character(knowledge_fn))
  if (.Platform$OS.type == "windows") {
    stop("Testing sequencing data does not support Windows.")
  }
  yame_result <- .Call("yame_summary_cfunc", query_fn, knowledge_fn)

  summary_results <- utils::read.table(text = paste(yame_result, collapse = "\n"), header = TRUE)

  summary_results <- summary_results %>%
    dplyr::select('Query','Mask','Beta') %>%
    tidyr::spread(key=Mask,value=Beta,fill = NA)

  column_order <- order(as.numeric(stringr::str_extract(colnames(summary_results), "\\d+")))
  summary_results <- summary_results[, column_order]
  rownames(summary_results) <- summary_results$Query
  summary_results <- summary_results %>% dplyr::select(-"Query")
  summary_results
}
