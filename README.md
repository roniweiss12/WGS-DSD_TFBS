# WGS-DSD_TFBS
This repository contains all the scripts used to extract TFBS data from various DBs (JASPAR, Homer, etc.) or predict them (TFBStools).

## Commonly Used Scripts:
Want to predict TFBS for a specific sequence? use `TFBStools/TFBStools_prediction.R`

Want to predict variant effect? use `TFBStools/variantEffectPrediction.R`

Do you already have a TFBS file and want to extract the results within specific genomic intervals? use `extractTFBS.R`

Are you looking for a script that changes a fasta file according to given variants? use `mutate_fasta.R`

The rest of the scripts have explanations within the code.
