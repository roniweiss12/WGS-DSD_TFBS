#!/usr/bin/env python
import pandas as pd
import sys

keyTFpath = sys.argv[1]
peaksMotifsPath = sys.argv[2]

tf = pd.read_csv(keyTFpath, header=None)[0].to_list()

peaks_motif = pd.read_csv(peaksMotifsPath, sep='\t')

relevant_motifs = peaks_motif[peaks_motif['Motif Name'].str.contains('|'.join(f'^{t}' for t in tf), case=False)]

relevant_motifs.to_csv("keyTF_motifs_in_peaks.csv")
