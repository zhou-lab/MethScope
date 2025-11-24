#' Generate pattern level data for cell type annotation
#'
#' @param query_fn File path to query .cg
#' @param knowledge_fn File path to pattern file .cm
#' @return A cell by pattern matrix.
#' @useDynLib MethScope, .registration = TRUE
#' @importFrom dplyr select
#' @importFrom stringr str_extract
#' @importFrom tidyr spread
#' @importFrom utils read.table
#' @importFrom magrittr %>%
#' @importFrom data.table fread
#' @export
GenerateInput <- function(query_fn, knowledge_fn) {

  stopifnot(is.character(query_fn), is.character(knowledge_fn))
  if (.Platform$OS.type == "windows") {
    stop("Testing sequencing data does not support Windows. Please directly use yame to generate inputs.")
  }
  #yame_result <- .Call("yame_summary_cfunc", query_fn, knowledge_fn)
  
  temp_file <- tempfile(fileext = ".txt")
  .Call("yame_summary_cfunc", query_fn, knowledge_fn, temp_file)
  summary_results <- data.table::fread(temp_file, header = TRUE)
  on.exit(unlink(temp_file), add = TRUE)
  summary_results <- summary_results %>%
    dplyr::select('Query','Mask','Beta') %>%
    tidyr::spread(key=Mask,value=Beta,fill = NA)

  column_order <- order(as.numeric(stringr::str_extract(colnames(summary_results), "\\d+")))
  summary_results <- summary_results[, column_order]
  rownames(summary_results) <- summary_results$Query
  summary_results <- summary_results %>% dplyr::select(-"Query")
  summary_results
}


#' Generate reference pattern file
#'
#' @param binary_file File path from the output of the bash script 
#' @param min_CG The minimum CpG # a pattern must have (Default: 50)
#' @param output_path File path to store the output txt file (Default: current path)
#' @return NULL, write out a patterns.txt file
#' @export

GenerateReference <- function(binary_file,min_CG = 50,output_path="./patterns.txt"){
  reference_set <- readLines(binary_file)
  # Calculate the frequency of each string
  stringFrequency <- table(reference_set)
  frequencyDF <- as.data.frame(stringFrequency)
  # Sort the data frame by frequency
  sortedFrequency <- frequencyDF[order(-frequencyDF$Freq),]
  sortedFrequency$pattern <- paste0("P",seq(1:nrow(sortedFrequency)))
  # If a pattern has less than certain number of CpG covered, group them into one Pna
  if(length(which(sortedFrequency$Freq <= min_CG)) > 0){
    sortedFrequency[which(sortedFrequency$Freq <= min_CG),]$pattern <- "Pna"}
  patterns_out <- sortedFrequency$pattern[match(reference_set, sortedFrequency$reference_set)]
  write.table(patterns_out,output_path,quote=FALSE,col.names = FALSE,row.names = FALSE)
}
