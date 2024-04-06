library(TFBSTools)
library(Biostrings)
library(seqinr)
library(dplyr)
library(arrow)
library(GenomicFeatures)
library(foreach)
library(doParallel)

# Usage function
usage <- function() {
  cat("Usage: Rscript TFBStools_prediction.R tf_jaspar_file fasta_file thresholds_file output_dir\n")
}

# Check command line arguments
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 4) {
  usage()
  quit(status = 1)
}

create_output_df <- function(rows_num, cols, df){
  output_df <- data.frame(matrix(nrow = rows_num, ncol = length(cols)))
  colnames(output_df) <- cols
  output_df$chr <- sapply(combined_df$seqname, function(x) unlist(strsplit(x, ":"))[1])
  output_df$coords <- sapply(combined_df$seqname, function(x) unlist(strsplit(x, ":"))[2])
  output_df$intStart <- sapply(output_df$coords, function(x) unlist(strsplit(x, "-"))[1])
  output_df$rest <- sapply(output_df$coords, function(x) unlist(strsplit(x, "-"))[2])
  output_df$intEnd <- sapply(output_df$rest, function(x) unlist(strsplit(x, "(", fixed = TRUE))[1])

  output_df$start <- as.integer(output_df$intStart) + combined_df$start
  output_df$end <- as.integer(output_df$intStart) + combined_df$end
  output_df$ID <- combined_df$TF_ID
  output_df$name <- combined_df$TF_name
  output_df$score <- combined_df$relative_score
  output_df$strand <- combined_df$strand
  output_df$source <- "TFBSTools"

  output_df <- output_df[cols]
  return(output_df)
}
handle_combined_df <- function(combined_df){
  combined_df <- cbind(combined_df, do.call(rbind, strsplit(combined_df$attributes, ";")))
  combined_df$sequence <- gsub("sequence=", "", combined_df[["3"]])
  combined_df$TF <- gsub("TF=", "", combined_df[["1"]])
  tf_split <- strsplit(combined_df$TF, "\\.")
  combined_df$TF_ID <- sapply(tf_split, function(x) paste(x[1], x[2], sep = "."))
  combined_df$TF_name <- sapply(tf_split, function(x) x[3])
  return(combined_df)
}
find_mtx <- function(matrix, seqs, threshold, output_dir){
  sitesList <- list()

  sitesList <- foreach(seq_name = names(seqs), .packages = c("TFBSTools", "Biostrings")) %dopar% {
    seq <- DNAString(seqs[[seq_name]][[1]])
    sites <- searchSeq(matrix, seq, seqname = seq_name, min.score = threshold, strand = "*")
    scores_list <- unlist(relScore(sites))
    gff <- writeGFF3(sites)
    gff$relative_score <- scores_list
    gff
  }

  combined_df <- bind_rows(sitesList)
  rownames(combined_df) <- NULL

  cols <- c("chr", "start", "end", "ID", "score", "strand", "name", "source")
  combined_df <- handle_combined_df(combined_df)
  output_df <- create_output_df(nrow(combined_df), cols, combined_df)
  output_name <- paste0(output_dir,"/", name(matrix), ".prq")
  write_parquet(output_df, output_name)
}
# Parse command line arguments
tf_jaspar_file <- args[1]
fasta_file <- args[2]
thresholds_file <- args[3]
output_dir <- args[4]

matrices <- readJASPARMatrix(tf_jaspar_file, matrixClass="PWM")
seqs <- read.fasta(fasta_file, seqtype = "DNA", as.string = TRUE)
thresholds <- read.csv(thresholds_file, header = TRUE)

# Parallel setup
cores <- detectCores()
cl <- makeCluster(cores)
registerDoParallel(cl)

# Process matrices in parallel
foreach(matrix = matrices, .packages = c("TFBSTools", "Biostrings")) %dopar% {
  mtx_name <- name(matrix)
  mtx_name <- strsplit(mtx_name, "\\.")[[1]]
  id <- paste0(mtx_name[1], '.', mtx_name[2])
  thresh <- thresholds[thresholds$matrixID == id, "threshold"]
  find_mtx(matrix, seqs, thresh, output_dir)
}

# Clean up parallel processing
stopCluster(cl)
