# Load command-line arguments
args <- commandArgs(trailingOnly = TRUE)

# Check if any arguments were passed
if (length(args) == 0) {
  cat("No arguments provided.\n 
      Usage: Rscript mutate_fasta.R variants.csv atac_intervals.fasta output_file.fasta interval_id")
}
library(Biostrings)
library(seqinr)
library(dplyr)

filter_seqs <- function(interval, seqs) {
    # Vectorized approach to check if the name of each sequence is in the interval
    seq_names <- names(seqs)
    
    # Extract the part of the name after ")"
    split_names <- sapply(seq_names, function(name) {
        split_name <- strsplit(name, split = ")")[[1]][2]
        return(split_name)
    })
    
    # Determine which indices match the interval
    matching_indices <- which(split_names %in% interval)
    
    # Return the indices that match the interval
    return(matching_indices)
}
get_seq <- function(x){
    #convert list of chars to string
  seq <- paste(as.character(getSequence(x)), collapse = "")
  seq <- paste(as.character(seq), collapse = "")
  return(seq)
}
change_sequence <- function(orig_seq, var_index, ref, alt){
    #This function changes the sequence according to the variant, including all use cases
  ref_len <- nchar(ref)
  alt_len <- nchar(alt)
  gap <- ref_len - alt_len
  if (gap == 0) {
    for (i in 1:ref_len){
      orig_seq[[1]][var_index + i - 1] <- alt[i]
    }
  } else if (gap > 0){
    replace_seq <- orig_seq
    rest <- orig_seq[[1]][(var_index + ref_len):length(orig_seq[[1]])]
    if (alt == ""){
      replace_seq[[1]][var_index: (length(orig_seq[[1]]) - gap - 1)] <- rest
      orig_seq[[1]] <- replace_seq[[1]][1:(length(orig_seq[[1]]) - gap)]
    } else {
      for (i in 1:alt_len){
        replace_seq[[1]][(var_index + i - 1)] <- substr(alt, i, i)
      }
      replace_seq[[1]][(var_index + alt_len): (length(orig_seq[[1]]))] <- rest
      orig_seq[[1]] <- replace_seq[[1]][1:(length(orig_seq[[1]]) - gap)]
    }
    
  } else {
    rest <- orig_seq[[1]][var_index + ref_len:length(orig_seq[[1]])]
    for (i in 1:alt_len){
      orig_seq[[1]][var_index + i - 1] <- substr(alt, i, i)
    }
    orig_seq[[1]][(var_index + alt_len): (length(orig_seq[[1]]) - gap)] <- rest
  }
  return(orig_seq)
}
get_var_data <- function(row, seqs){
    #This function gets all relevant data for a given variant
  data <- c()
  data[["pos"]] <- row[["POS"]]
  data[["ref"]] <- row[["REF"]]
  data[["alt"]] <- row[["ALT"]]
  if (data[["alt"]] == "*"){
    data[["alt"]] = ""
  }
  data[["interval"]] <- row[["INTERVAL_ID"]]
  rex <- paste0(".*", data[["interval"]], "$")
  ref_seq <- seqs[grepl(rex, names(seqs))]
  alt_seq <- ref_seq
  full_name <- names(ref_seq)
  data[["chr"]] <- unlist(strsplit(full_name, split = ":"))[1]
  rest <- unlist(strsplit(full_name, split = ":"))[2]
  data[["start_interval"]] <- as.numeric(unlist(strsplit(rest, split = "-"))[1])
  rest <- unlist(strsplit(rest, split = "-"))[2]
  data[["end_interval"]] <- as.numeric(unlist(strsplit(rest, split = "\\("))[1])
  var_index <- data[["pos"]] - data[["start_interval"]]
  alt_seq <- change_sequence(alt_seq, var_index, data[["ref"]], data[["alt"]])
  names(alt_seq) <- paste0(names(alt_seq), "_mt")
  data[["seqs_to_scan"]] <- c(ref_seq, alt_seq)
  return(data)
}

variants_file <- args[1]
fasta_file <- args[2]
output_file <- args[3]
interval_id <- args[4]

seqs <- read.fasta(fasta_file)
vars <- read.csv(variants_file, header = TRUE)
vars <- vars[vars$INTERVAL_ID == interval_id,]

# # Filter sequences - keep only requested interval
keep_indices <- filter_seqs(interval_id, seqs)
seq <- seqs[keep_indices]

output_seqs <- list()
seqs_names <- list()
for (i in 1:nrow(vars)){
  data <- get_var_data(vars[i,], seq)
  mt_seq <- data[["seqs_to_scan"]][[2]]
  seq_name <- paste0(data[["interval"]], "_", i, "_", data[["chr"]], ":", 
                    data[["pos"]], "_", data[["ref"]], "_", data[["alt"]])
  seqs_names <- c(seqs_names, seq_name)
  mt_seq <- paste(mt_seq, collapse = "")
  output_seqs[[i]] <- mt_seq
}
write.fasta(output_seqs, seqs_names, file.out = output_file, open = "w")
