import pyBigWig
import pandas as pd
import sys
import os

CHR_COL = 0
START_COL = 1
END_COL = 2

def process_bed_entry(row, bb):
    """
    This function extracts the coordinates of each interval and gets the relevant entries from the bigbed file
    """
    chromosome = row[CHR_COL]
    start = int(row[START_COL])
    end = int(row[END_COL])
    result = bb.entries(chromosome, start, end)
    return result

def main(bigbed_file, intervals_file, output_file):
    #read files
    bb = pyBigWig.open(bigbed_file)
    intervals = pd.read_table(intervals_file, header=None)
    #apply processing function
    results = intervals.apply(process_bed_entry, args=(bb,), axis=1)
    # Creating a DataFrame for 'results'
    results_df = pd.DataFrame({'Results': results})
    # Adding the 'intervals'['chr'] column to the 'results_df'
    results_df['chr'] = intervals[CHR_COL]
    # Writing the combined data to a CSV file
    results_df.to_csv(output_file, index=False)


def check_file_exists(file_path):
    if not os.path.exists(file_path):
        print(f"Error: File '{file_path}' does not exist.")
        return False
    return True  

if __name__ == "__main__":
    if len(sys.argv) != 4:
        print("Usage: python intersect_intervals_to_bigBed.py bigbed_file intervals_file output_file")
        exit(0)
    
    bigbed_file = sys.argv[1]
    intervals_file = sys.argv[2]
    output_file = sys.argv[3]

    if not (check_file_exists(bigbed_file) and check_file_exists(intervals_file)):
        print("Some or all of the files don't exist")
        exit(0)

    main(bigbed_file, intervals_file, output_file)