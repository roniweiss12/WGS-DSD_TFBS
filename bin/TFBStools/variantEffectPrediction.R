# Load command-line arguments
args <- commandArgs(trailingOnly = TRUE)

# Check if any arguments were passed
if (length(args) == 0) {
  cat("No arguments provided.\n 
      Usage: Rscript variantEffectPrediction.R selected_matrices.jaspar variants.csv atac_intervals.fasta matrices_thresholds.csv")
}
library(TFBSTools)
library(Biostrings)
library(seqinr)
library(dplyr)
library(arrow)
library(GenomicFeatures)

## CONSTANTS
SEQNAMES_FILEDS <- c("V1", "seqname", "chromosome", "chrom", "chr", "chromosome_name", "seqid", "CHROM")
START_FIELD <- c("V2", "POS", "start")
END_FIELD <- c("V3", "stop", "end")
HUMAN_CELLS <- c("PG", "PS", "SE", "GR")
THRESHOLD <- 0

filter_seqs <- function(intervals_set, seqs) {
    # Vectorized approach to check if the name of each sequence is in the interval
    seq_names <- names(seqs)
    
    # Extract the part of the name after ")"
    split_names <- sapply(seq_names, function(name) {
        split_name <- strsplit(name, split = ")")[[1]][2]
        return(split_name)
    })
    
    # Determine which indices match the interval
    matching_indices <- which(split_names %in% intervals_set)
    
    # Return the indices that match the intervals_set
    return(matching_indices)
}
get_seq <- function(x){
  seq <- paste(as.character(getSequence(x)), collapse = "")
  seq <- paste(as.character(seq), collapse = "")
  return(seq)
}

split_att <- function(df){
  df <- cbind(df, do.call(rbind, strsplit(df$attributes, ";")))
  df$sequence <- gsub("sequence=", "", df[["3"]])
  df$TF <- gsub("TF=", "", df[["1"]])
  tf_split <- strsplit(df$TF, "\\.")
  df$TF_ID <- sapply(tf_split, function(x) paste(x[1], x[2], sep = "."))
  df$TF_name <- sapply(tf_split, function(x) x[3])
  cols_to_remove <- c("attributes", "1", "2", "3", "TF")
  df <- df[, -which(names(df) %in% cols_to_remove)]
  return(df)
}
  
find_mtx <- function(matrices, seqs, threshold){
  sitesList <- c()
  
  for (i in 1:length(seqs)) {
    seq <- seqs[i]
    seq_name <- names(seq)
    seq <- as(vapply(seq, get_seq, character(1)), "DNAStringSet")
    seq <- paste(as.character(seq), collapse = "")
    sites <- searchSeq(matrices, seq, seqname = seq_name, min.score = threshold, strand = "*")
    scores_list <- unlist(relScore(sites))
    gff <- writeGFF3(sites)
    gff$relative_score <- scores_list
    gff <- gff[,-1]
    gff <- split_att(gff)
    rownames(gff) <- NULL
    sitesList[[seq]] <- gff
  }
  #combined_df <- bind_rows(sitesList)
  return (sitesList)
}
change_sequence <- function(orig_seq, var_index, ref, alt){
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
  data <- c()
  data[["pos"]] <- as.numeric(row[["POS"]])
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
# Define a function to apply row-wise
find_shift <- function(row, df) {
  # Check if relative_score.x is not NA
  if (!is.na(row["relative_score.x"])) {
    # Find the index of the row in df where relative_score.y matches relative_score.x from the input row
    matching_index <- which(df["relative_score.y"] == row["relative_score.x"][[1]])
    # Check if a matching index is found
    if (length(matching_index) > 0) {
      # Return the first matching row from df
      return(c(rownames(row), rownames(df[matching_index,])))
    }
  }
  # Return NULL if no match is found or relative_score.x is NA
  return(NULL)
}

filter_by_threshold <- function(thresholds, df){
  combined_df <- merge(df, thresholds, by.x = "TF_ID", by.y = "matrixID")
  combined_df <- subset(combined_df, relative_score.x >= threshold | relative_score.y >= threshold)
  cols_to_remove <- c("TF", "X")
  combined_df <- combined_df[, -which(names(combined_df) %in% cols_to_remove)]
  rows_to_remove <- lapply(unique(combined_df$TF_ID), function(id) {
    ss <- subset(combined_df, TF_ID == id)
    double_rows <- lapply(seq_len(nrow(ss)), function(row_idx) {
      find_shift(ss[row_idx, ], ss)
    })
    double_rows
  })
  rows_to_remove <- as.numeric(unlist(rows_to_remove))
  filtered_df <- combined_df[!(rownames(combined_df) %in% rows_to_remove), ]
  return(filtered_df)
}
find_overlap <- function(df){
  overlaps <- data.frame(pair1 = integer(), pair2 = integer())
  rows_to_remove <- c()
  # Loop through each pair of intervals
  for (i in 1:(nrow(df) - 1)) {
    for (j in (i + 1):nrow(df)) {
      # Check if intervals overlap
      if ((df$start[i] <= df$end[j]) && (df$end[i] >= df$start[j])) {
        overlaps <- rbind(overlaps, data.frame(pair1 = i, pair2 = j))
      }
    }
  }
  if(nrow(overlaps) > 0){
    # Identify rows to remove efficiently
    rows_to_remove_x <- ifelse(df$relative_score.x[overlaps$pair1] >= df$relative_score.x[overlaps$pair2],
                             rownames(df[overlaps$pair2,]), rownames(df[overlaps$pair1,]))
    # Remove duplicate rows to avoid unnecessary multiple removals
    rows_to_remove_x <- unique(rows_to_remove_x)
    df <- df[!(rownames(df) %in% rows_to_remove_x), ]
    if(nrow(df) > 1){
      overlaps <- data.frame(pair1 = integer(), pair2 = integer())
      
      # Loop through each pair of intervals
      for (i in 1:(nrow(df) - 1)) {
        for (j in (i + 1):nrow(df)) {
          # Check if intervals overlap
          if ((df$start[i] <= df$end[j]) && (df$end[i] >= df$start[j])) {
            overlaps <- rbind(overlaps, data.frame(pair1 = i, pair2 = j))
          }
        }
      }
      rows_to_remove_y <- ifelse(df$relative_score.y[overlaps$pair1] >= df$relative_score.y[overlaps$pair2],
                                 rownames(df[overlaps$pair2,]), rownames(df[overlaps$pair1,]))
      # Remove duplicate rows to avoid unnecessary multiple removals
      rows_to_remove_y <- unique(rows_to_remove_y)
      rows_to_remove <- c(rows_to_remove_x, rows_to_remove_y)
    } else {
    rows_to_remove <- rows_to_remove_x
    }
  }
  return(rows_to_remove)
}
filter_overlapping_tf <- function(df) {
  rows_to_remove <- integer(0)
  # Iterate over unique TF_IDs
  for (id in unique(df$TF_ID)) {
    ss <- subset(df, TF_ID == id)
    if(nrow(ss) > 1){
      rows <- find_overlap(ss)
      rows_to_remove <- c(rows_to_remove, rows)
    }
  }
  return(rows_to_remove)
}
get_abs_difference <- function(row){
  if (is.na(row[["relative_score.x"]])) {
    value <- as.numeric(row[["relative_score.y"]]) - THRESHOLD
  } else if (is.na(row[["relative_score.y"]])){
    value <- as.numeric(row[["relative_score.x"]]) - THRESHOLD
  } else {
    value <- abs(as.numeric(row[["relative_score.y"]]) - as.numeric(row[["relative_score.x"]]))
  }
  return(value)
}
get_delta_cols <- function(col_names){
  cols <- c()
  for (col_name in col_names) {
    # Check if the column name contains 'delta'
    if (grepl("delta", col_name)) {
      cols <- c(cols, col_name)
    }
  }
  return(cols)
}
tfbs_delta <- function(row, cols){
  # Initialize a vector to store delta values for each row
  row_deltas <- c()
  # Loop through each column
  for (col_name in cols) {
    # Split the string by "_"
    col_values <- unlist(strsplit(as.character(row[col_name]), "_"))
    # Convert values to numeric
    col_values <- as.numeric(col_values)
    # Calculate delta and append to row_deltas
    row_deltas <- c(row_deltas, col_values)
  }  
  # Take the maximum absolute delta value for the row
  max_deltas <- ifelse(length(row_deltas) > 1, max(abs(row_deltas)), 0)
  return(max_deltas)
}
# function to check if any values in a string are contained in the list
contains_human_cells <- function(cell_string) {
  any(sapply(HUMAN_CELLS, grepl, cell_string))
}
blend_rows <- function(df, gap){
  x_rows <- df[is.na(df$relative_score.y), ]
  y_rows <- df[is.na(df$relative_score.x), ]
  for(i in 1:nrow(x_rows)){
    for (j in 1:nrow(y_rows)) {
      if(((x_rows$start[i] - abs(gap)) == y_rows$start[j]) & ((x_rows$end[i] - abs(gap)) == y_rows$end[j])){
        x_rows$relative_score.y[i] <- y_rows$relative_score.y[j]
        x_rows$score.y[i] = y_rows$score.y[j]
        x_rows$abs_diff[i] = abs(x_rows$relative_score.y[i] - x_rows$relative_score.x[i])
        df[row.names(x_rows[i]),] <- x_rows[i,]
        df <- df[!(row.names(df) %in% row.names(y_rows[j])), ]
      }
    }
  }
  return(df)
}
main_pipeline <- function(row){
  #get variant data
  data <- get_var_data(row, seqs)
  #scan matrices for wt and mt
  res <- find_mtx(matrices, data[["seqs_to_scan"]], THRESHOLD)
  #res[[2]] <- apply_shift(res[[2]], data)
  #find differential rows
  cols <- c("start", "end", "score", "strand", "relative_score", "TF_ID", "TF_name")
  ref_non_alt <- anti_join(res[[1]], res[[2]], by = colnames(res[[1]]))[,cols]
  alt_non_ref <- anti_join(res[[2]], res[[1]], by = colnames(res[[1]]))[,cols]
  #alt_non_ref <- apply_shift(alt_non_ref, data)
  changes <- merge(x = ref_non_alt, y = alt_non_ref, by = c("start", "end", "TF_ID", "TF_name", "strand"), all = TRUE)
  #calculate difference
  changes$abs_diff <- apply(changes, 1, get_abs_difference)
  #keep only rows with difference of over 1%
  changes <- subset(changes, abs_diff > 0.01 | is.na(abs_diff))
  #remove rows where neither wt or mt were above threshold
  filtered_df <- filter_by_threshold(thresholds, changes)
  if(nrow(filtered_df) > 0){
    for (id in unique(filtered_df$TF_ID)){
      tf_row <- subset(filtered_df, TF_ID == id)
      #remove overlapping TFBS of same TF
      rows_to_remove <- unlist(filter_overlapping_tf(tf_row))
      if (length(rows_to_remove) > 0){
        tf_row <- tf_row[!(rownames(tf_row) %in% rows_to_remove), ]
      }
      if (nrow(tf_row) == 2){
        if (tf_row$start[1] == tf_row$start[2]){
          tf_row <- tf_row[1, ]
        }
      }
      tf <- unique(tf_row$TF_name)
      func <- ifelse(is.na(tf_row$relative_score.x), "gain", 
                     ifelse(is.na(tf_row$relative_score.y), "loss",
                            ifelse(tf_row$relative_score.x > tf_row$relative_score.y, "loss", "gain")))
      if(length(func) > 1){
        result <- character()  # Initialize an empty character vector
        result_orig <- character()
        result_mut <- character()
        delta <- character()
        # Iterate through each row of the data frame
        for (j in 1:nrow(tf_row)) {
          # Concatenate func[j], tf_row$start[j], tf_row$end[j] with the specified separators
          # and append to the result vector
          result <- c(result, paste0(func[j], ";", tf_row$start[j], ";", tf_row$end[j]))
          result_orig <- c(result_orig, tf_row$relative_score.x[j])
          result_mut <- c(result_mut, tf_row$relative_score.y[j])
          delta <- c(delta, tf_row$abs_diff[j])
        }
        row[paste0(tf, "_", id)] <- paste0(result, collapse = "_")
        row[paste0(tf, "_", id, "_original")] <- paste0(result_orig, collapse = "_")
        row[paste0(tf, "_", id, "_mutated")] <- paste0(result_mut, collapse = "_")
        row[paste0(tf, "_", id, "_delta")] <- paste0(delta, collapse = "_")
      } else {
        row[paste0(tf, "_", id)] <- paste0(func, ";", tf_row$start, ";", tf_row$end, "_")
        row[paste0(tf, "_", id, "_original")] <- tf_row$relative_score.x
        row[paste0(tf, "_", id, "_mutated")] <- tf_row$relative_score.y
        row[paste0(tf, "_", id, "_delta")] <- tf_row$abs_diff
      }
    }
  }
  return(row)
}

tf_jaspar_file <- args[1]
variants_file <- args[2]
fasta_file <- args[3]
threshold_file <- args[4]
output_file <- strsplit(variants_file, "\\.")[[1]][1]
output_file <- paste0(output_file, "_effect_prediction.tsv")

matrices <- readJASPARMatrix(tf_jaspar_file, matrixClass="PWM")
seqs <- read.fasta(fasta_file)
vars <- read.csv(variants_file, header = TRUE)
thresholds <- read.csv(threshold_file)

# Filter sequences - keep only those that have a variant
intervals_set <- unique(vars$INTERVAL_ID)
keep_indices <- filter_seqs(intervals_set, seqs)
seqs <- seqs[keep_indices]

thresholds <- thresholds[order(thresholds$TF),]
for (id in thresholds$matrixID){
  tf <- thresholds[thresholds$matrixID == id, "TF"]
  vars[[paste0(tf, "_", id)]] <- ""
  vars[[paste0(tf, "_", id, "_original")]] <- ""
  vars[[paste0(tf, "_", id, "_mutated")]] <- ""
  vars[[paste0(tf, "_", id, "_delta")]] <- ""
}

# Apply the main_pipeline function to each row of the dataframe
modified_rows <- apply(vars, 1, main_pipeline)
# Assign the modified rows back to the original dataframe
vars <- as.data.frame(t(modified_rows))
                       
delta_cols <- get_delta_cols(names(vars))
# find each row's maximal delta
vars$TFBS_delta <- apply(vars, 1, tfbs_delta, cols = delta_cols)
#order by highest delta
vars <- vars[order(vars$TFBS_delta, decreasing = TRUE), ]
#order by interval source
vars$contains_human_cells <- sapply(vars$INTERVAL_ID, contains_human_cells)
vars <- vars[order(vars$contains_human_cells, decreasing = TRUE), ]
# Remove the temporary column
vars$contains_human_cells <- NULL
#output_file
write.table(vars, "tmp.tsv", sep = '\t', row.names = FALSE, quote = FALSE, na = "")
