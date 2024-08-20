import pandas as pd
import sys
## This script converts variants in a CSV to a 0-based BED file.
#usage
if len(sys.argv) != 1:
    print('Usage: python vars_csv_to_bed.py <vars_file.csv>')
    sys.exit(0)

vars_file = sys.argv[1]
df = pd.read_csv(vars_file)
df['start'] = df['POS'] - 1
df['end'] = df['POS']
#df['rsID'] = '-'

new_df = df[['CHROM', 'start', 'end', 'ALT', 'REF', 'rsID', 'AF_popmax']]

vars_file = vars_file.replace('.csv', '.bed')
new_df.to_csv(vars_file, sep='\t', index=False, header=False)
