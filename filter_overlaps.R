## This script filters out overlapping TFBS of the same TF and keeps the one with the highest score.

# Load command-line arguments
args <- commandArgs(trailingOnly = TRUE)

# Check if any arguments were passed
if (length(args) == 0) {
  cat("No arguments provided.\n 
      Usage: Rscript filter_output.R path/to/tfbs.prq path/to/output.prq")
}
library(arrow)
library(GenomicRanges)

# filter each transcription factor separately
process_tf <- function(tf_df) {
  gr <- makeGRangesFromDataFrame(tf_df, 
                                 seqnames.field = SEQNAMES_FILEDS,
                                 start.field = START_FIELD, 
                                 end.field = END_FIELD)
  
  overlaps <- findOverlaps(gr, gr, type = "any")
  
  query_hits <- queryHits(overlaps)
  subject_hits <- subjectHits(overlaps)
  
  rows_to_remove <- ifelse(query_hits == subject_hits, 
                           NA,
                           ifelse(tf_df$score[query_hits] >= tf_df$score[subject_hits],
                                  subject_hits, query_hits)
  )
  rows_to_remove <- na.omit(rows_to_remove)
  rows_to_remove <- unique(rows_to_remove)
  
  tf_df <- tf_df[-rows_to_remove, ]
  return(tf_df)
}

## CONSTANTS
SEQNAMES_FILEDS <- c("V1", "seqname", "chromosome", "chrom", "chr", "chromosome_name", "seqid", "CHROM")
START_FIELD <- c("V2", "POS", "start")
END_FIELD <- c("V3", "stop", "end")

input_file <- args[1]
output_file <- args[2]

df <- read_parquet(input_file)
print(paste("original nrow:", nrow(df)))
df_list <- list()
tf_list <- unique(df$name)

# Apply the function to each transcription factor's data
df_list <- lapply(tf_list, function(tf) {
  tf_df <- df[df$name == tf, ]
  df <- df[df$name != tf, ]
  return(process_tf(tf_df))
})

final_df <- do.call(rbind, df_list)
print(paste("final nrow:", nrow(final_df)))
write_parquet(final_df, output_file)
