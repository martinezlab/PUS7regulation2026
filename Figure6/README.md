# Figure 6: cell type specificity

This folder has code and relevant files for the cell type specificity analysis, specifically those that are unique to Figure 6.

Cell type specific site calling code is located in the Figure 3 folder.
Plotting of various sites uses code highly similar to that used to plot sites in Figure 3.

Sequence and structure analysis uses the same pipelines as established in the Figure 4 folder, just with different filtering criteria.

6G was create with STREME (https://meme-suite.org/meme/doc/streme.html). S6H was created with the GEPIA web server (http://gepia.cancer-pku.cn/).

## Read Count Normalization

To generate plots S6C and S6D, run `read_norm.R`, which has data paths hard-coded into it.

## Imaging Quantification

To quantify the nuclear to cytoplasmic ratios, Cellpose (https://www.cellpose.org/) was used to draw boundaries based on DAPI and brightfield. 
PUS7 (FITC) intensity was then quantified within these boundaries and ratio calculated.
The Python script `extract_ratios.py` is used to do this; adjust directory paths within the script itself.
The script processes one image at a time and then aggregates the results from multiple images to a single CSV file.

The above quantification pipeline was adapted from Ruan et al, 2026 (https://github.com/mlruan/Imaging_data_quantify).


To generate the plot in S6F, run `plot_stats.R`, with hard-coded data. This shows data from 100 representative cells in one of the three biological replicates we performed (rep2_data.csv)