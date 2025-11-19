#' Generate pattern reference
#'
#' @param binary_file File path from the output of the bash script 
#' @param min_CG The minimum CpG # a pattern must have (Default: 50)

args <- commandArgs(trailingOnly = TRUE)

binary_file <- args[1]
min_CG <- as.numeric(args[2])

GenerateReference <- function(binary_file,min_CG = 50){
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
  write.table(patterns_out,"./patterns.txt",quote=FALSE,col.names = FALSE,row.names = FALSE)
}

GenerateReference(binary_file,min_CG)