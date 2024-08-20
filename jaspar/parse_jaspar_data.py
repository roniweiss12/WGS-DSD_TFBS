# This script parses the data extracted from JASPAR2024.bb
import pandas as pd
import ast
from pandarallel import pandarallel
import sys
import os
import numpy as np
import pyarrow

LIST_COL = 0
DATA_COL = 'Results'

def convert_to_tuple(val):
    return ast.literal_eval(val)

def process_data(row):
    '''
    Convert the input unparsed data to a tuple
    '''
    val = row[DATA_COL]
    row_tuple = ast.literal_eval(val)
    row[DATA_COL] = row_tuple
    return row

def parse_jaspar_data(df):
    '''
    This function converts each row, seperates the different entries to multiple rows and splits the details of each entry to multiple columns
    '''
    # Initialize pandarallel
    pandarallel.initialize(nb_workers=8)
    # Assuming 'df' is your DataFrame
    df_parallel = df.parallel_apply(process_data, axis=1)
    # Explode the DataFrame based on the tuples in the DATA_COL column
    df_exploded = df_parallel.explode(DATA_COL)
    # Repeat 'chr' column values for each exploded row
    chr_column_expanded = np.repeat(df_parallel['chr'], df_parallel[DATA_COL].str.len())
    # Split the tuples into separate columns in the exploded DataFrame
    exploded_details = pd.DataFrame(df_exploded[DATA_COL].tolist(), columns=['start', 'end', 'details'])
    # Concatenate 'chr' column and exploded details horizontally
    df_exploded = pd.concat([pd.DataFrame(chr_column_expanded.reset_index(drop=True), columns=['chr']), exploded_details.reset_index(drop=True)], axis=1)
    # Further split 'Value3' by '\t'
    value3_split = df_exploded['details'].str.split('\t', expand=True)
    value3_split.columns = ['ID', 'score', 'strand', 'name']
    # Concatenate the split columns with the original DataFrame
    df_final = pd.concat([df_exploded, value3_split], axis=1)
    df_final.drop('details', axis=1, inplace=True)  # Drop the original 'Value3' column
    return df_final

def filter_tf(df, tf_list, score_threshold):
    '''
    Check if the entry is in our list of TF families, and the p-value is over the threshold
    '''
    print(len(df))
    df_name_series = df['name'] if isinstance(df['name'], pd.Series) else pd.Series(df['name'])
    # Filter out rows where 'name' is not in tf_list
    df_filtered = df[df_name_series.isin(tf_list)]
    print(len(df_filtered))
    df_filtered = df_filtered[df_filtered['score'].astype(int) >= score_threshold]
    return df_filtered

def main(tf_file, jaspar_non_parsed, output_file, score_threshold=400):
    tf_list = pd.read_csv(tf_file, header=None)

    chunk_size = 10000
    parsed_dfs = []
    for chunk in pd.read_csv(jaspar_non_parsed, chunksize=chunk_size):
        parsed_df = parse_jaspar_data(chunk)
        filtered_df = filter_tf(parsed_df, tf_list[LIST_COL], score_threshold)
        parsed_dfs.append(filtered_df)

    result_df = pd.concat(parsed_dfs)
    result_df['source'] = 'jaspar'
    result_df.to_parquet(output_file, index=False)

def check_file_exists(file_path):
    if not os.path.exists(file_path):
        print(f"Error: File '{file_path}' does not exist.")
        return False
    return True  

if __name__ == "__main__":
    if len(sys.argv) < 4:
        print("Usage: python parse_jaspar_data.py tf_file jaspar_non_parsed output_file score_threshold(optional)")
        exit(0)
    print(len(sys.argv))
    tf_file = sys.argv[1]
    jaspar_non_parsed = sys.argv[2]
    output_file = sys.argv[3]

    if not (check_file_exists(tf_file) and check_file_exists(jaspar_non_parsed)):
        print("Some or all of the files don't exist")
        exit(0)

    if len(sys.argv) == 5:
        score_threshold = int(sys.argv[4])
        main(tf_file, jaspar_non_parsed, output_file, score_threshold)
    else:
        main(tf_file, jaspar_non_parsed, output_file)
