#' Filter pattern reference
#'
#' @param binary_file File path from the output of the bash script 
#' @param min_CG The minimum CpG # a pattern must have (Default: 50)

args <- commandArgs(trailingOnly = TRUE)

binary_file <- args[1]
min_CG <- as.numeric(args[2])

FilterReference <- function(binary_file,min_CG = 50){
  reference_set <- readLines(binary_file)
  # Calculate the frequency of each string
  stringFrequency <- table(reference_set)
  frequencyDF <- as.data.frame(stringFrequency)
  # Sort the data frame by frequency
  sortedFrequency <- frequencyDF[order(-frequencyDF$Freq),]
  sortedFrequency$pattern <- paste0("P",seq(1:nrow(sortedFrequency)))
  # If a pattern has less than certain number of CpG covered, group them into one Pna
  sortedFrequency[which(sortedFrequency$Freq <= min_CG),]$pattern <- "Pna"
  patterns_out <- sortedFrequency$pattern[match(reference_set, sortedFrequency$reference_set)]
  idx_file <- read.table("tmp.cg.idx",sep="\t")
  a = sortedFrequency[sortedFrequency$pattern != "Pna",]
  a1 = cbind(a[,c("pattern","Freq")], a[,"reference_set",drop=F] %>% separate(reference_set, into = idx_file$V1, sep=seq_len(length(idx_file$V1)-1)))
  write.table(a1,"./pattern_definition.txt",quote=FALSE,col.names = TRUE,row.names = FALSE)
  dfdef = a1
  mxdef = as.matrix(dfdef[,3:ncol(dfdef)])
  rownames(mxdef) = dfdef[,1]
  metx_na = do.call(c, mclapply(0:24, function(gchunk) {
    mtx = as.matrix(read.table(text=system(sprintf("yame rowsub -1 -R ~/references/hg38/KYCGKB_hg38/cpg_nocontig.cr -I %d_1225000 tmp.cg | yame unpack -f 3 -C -r 2 -l tmp.cg.idx -a -", gchunk), intern=TRUE), header=T, row.names=1, check.names=FALSE))
    pat = read.table(text=system(sprintf("yame rowsub -I %d_1225000 patterns.cm | yame unpack -", gchunk), intern=T), header=F, check.names=FALSE)$V1
    sapply(seq_along(pat), function(i) { if (pat[i] == "Pna") { NA } else { v0 = mtx[i, mxdef[pat[i],] == 0]; v1 = mtx[i, mxdef[pat[i],] == 1];  if (length(v0) == 0 || length(v1) == 0) { 0 } else { max(sum(is.na(v0))/length(v0), sum(is.na(v1))/length(v1)) }}})
  }, mc.cores=12))
  metx_na_df <- as.data.frame(cbind(as.numeric(metx_na),reference_set))
  metx_na_df <- metx_na_df %>% mutate(reference_set = ifelse(is.na(V1) | V1 > 0, strrep("2", ncol(mxdef)), reference_set))
  reference_set <- metx_na_df$reference_set
  write.table(reference_set,"./tmp_filter.txt",quote=FALSE,col.names = FALSE,row.names = FALSE)
  stringFrequency <- table(reference_set)
  frequencyDF <- as.data.frame(stringFrequency)
  # Sort the data frame by frequency
  sortedFrequency <- frequencyDF[order(-frequencyDF$Freq),]
  sortedFrequency$pattern <- paste0("P",seq(1:(nrow(sortedFrequency)))-1)
  # If a pattern has less than certain number of CpG covered, group them into one Pna
  sortedFrequency[which(sortedFrequency$Freq <= min_CG),]$pattern <- "Pna"
  sortedFrequency$pattern <- ifelse(sortedFrequency$pattern == "P0","Pna",sortedFrequency$pattern)
  patterns_out <- sortedFrequency$pattern[match(reference_set, sortedFrequency$reference_set)]
  write.table(patterns_out,"./patterns_updated.txt",quote=FALSE,col.names = FALSE,row.names = FALSE)
}

FilterReference(binary_file,min_CG)