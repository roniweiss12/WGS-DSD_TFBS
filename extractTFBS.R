# Load command-line arguments
args <- commandArgs(trailingOnly = TRUE)

# Check if any arguments were passed
if (length(args) == 0) {
  cat("No arguments provided.\n 
      Usage: Rscript extractTFBS.R dataFile intervalFile outputFile.tsv")
}
library(GenomicRanges)
library(rtracklayer)
library(stringr)
library(readr)
library(arrow)
library(data.table)

## CONSTANTS
SEQNAMES_FILEDS <- c("V1", "seqname", "chromosome", "chrom", "chr", "chromosome_name", "seqid", "CHROM")
START_FIELD <- c("V2", "POS", "start")
END_FIELD <- c("V3", "stop", "end")
QUERY_COL <- 1
SUBJ_COL <- 2
INTERVAL_COL <- "V4"

data_file <- args[1]
query_file <- args[2]
output_file <- args[3]

## FUNCTIONS
overlap_regions <- function(query_object, subject_object) {
  overlaps <- findOverlaps(query_object, subject_object)
  query_hits <- queryHits(overlaps)
  subject_hits <- subjectHits(overlaps)
  len <- length(subject_hits)
  return(list(query_hits, subject_hits, len))
}
## WORKFLOW
data <- read_parquet(data_file)
data_obj <- makeGRangesFromDataFrame(data, seqnames.field = SEQNAMES_FILEDS,
                                 start.field = START_FIELD, end.field = END_FIELD)

query <- read.table(query_file, sep = "\t")
query_obj <- makeGRangesFromDataFrame(query, seqnames.field = SEQNAMES_FILEDS,
                                       start.field = START_FIELD, end.field = END_FIELD)


ovlp <- overlap_regions(query_obj, data_obj)
result_data <- data[ovlp[[SUBJ_COL]], ]
colnames(result_data) <- c("chr",	"start",	"end",	"ID_or_sequence",	"score",	"strand",	"name",	"source")
result_data$interval_name <- query[ovlp[[QUERY_COL]], INTERVAL_COL]

if (length(ovlp[[QUERY_COL]]) == 0){
  print("The coordinates you entered don't match any relevant TFBS in our ATAC intervals.")
} else {
  write.table(result_data, output_file, row.names = FALSE, quote = FALSE, sep = "\t")
}
