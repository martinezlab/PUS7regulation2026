# Figure 3: MPRA Nano-BID-Amp

This folder has code and files relevant to the MPRA Nano-BID-Amp pipeline and analysis.

## FASTQ to BAM prep

This pipeline is optimized to run on Stanford's HPC Sherlock system. Make sure you adjust have the following installed and adjust package loading appropriately:

- samtools 1.16.1
- cutadapt 1.18 (requires Python 3.6)
- java 11
- minimap2: https://github.com/lh3/minimap2
- UMICollapse: https://github.com/Daniel-Liu-c0deb0t/UMICollapse 

The point of is to go from a compressed fastq file, named as a barcode from Nanopore demultiplexing, into a mapped bam file that is ready for downstream processing. All scripts are intended to work on Sherlock, Stanford's HPC.

For the sample map file, it is **strongly** recommended to set up you file names to follow the scheme of celltype_vector_rep_treat for integration with the rest of the pipeline.

Example sample map file:

``` bash
HepG2_WT_in_1:barcode66
HepG2_WT_BS_1:barcode67
HepG2_WT_BS_2:barcode68
HepG2_KD_BS_2:barcode70
```

### submit_file_prep.sh

This Bash script is a general-purpose SLURM job array submitter to launch a batch of jobs on an HPC cluster. This script does not perform the actual data processing itself.

Its steps are:
1. Parse Configuration: It reads all settings from command-line arguments. These include SLURM resource requests (CPU, memory, time) and the paths required by the pipeline script it will run (e.g., input/output directories, sample map).
2. Validate Inputs: It ensures all required arguments are provided and that the specified script, files, and directories exist and are accessible.
3. Submit Job Array: It counts the number of lines in the sample map file to determine the required number of jobs. It then uses the sbatch command to submit a job array to SLURM, where each task in the array is configured with the user-specified resources.
4. Launch Pipeline: Each job in the submitted array will execute the user-provided pipeline script. This submitter script passes the necessary file paths and settings on to that pipeline script, which then processes one sample per job.

How to use:

```bash
bash ~/Figure3/file_prep/submit_file_prep.sh \
    -u you@institution.edu -t 1:00:00 \
    -s ~/file_prep/trim_map_dedup_mpra.sh \
    -a <barcodes.txt> \
    -i <input_fastq> \
    -o <output_directory> \
    -r <reference_fasta>
```

### trim_map_dedup_mpra.sh

This script processes a single Nano-BID-Amp sample as part of a SLURM array job. This script is appropriate for MPRA where there are two sets of adapters, the Nanopore adapters and the MPRA adapters, to trim from the reads. Nanopore adapters are trimmed prior to UMI extraction, MPRA adapters are trimmed following UMI extraction.

Its steps are:
1. Setup and Initialization: It parses command-line arguments, uses the SLURM array task ID to identify a specific sample from a map file, sets up output directories, and redirects all output to a log file.
2. Filter Reads: It uses cutadapt with linked adapters to process the compressed FASTQ file twice, creating two separate files: one containing "sense" reads and another with "antisense" reads, based on the Nanopore adapter sequences.
3. Unify Read Orientation: It reverse-complements the "antisense" reads and then concatenates them with the "sense" reads. This creates a single FASTQ file where all reads are oriented in the sense direction.
4. Extract UMIs: It uses umi_tools extract to pull 10-nucleotide Unique Molecular Identifiers (UMIs) from the 3' end of each read, embedding the UMI into the read's header for later use.
5. Second Trim: Performs a second cutadapt pass on the UMI-tagged reads to remove adapter sequences for the MPRA.
6. Align Reads: It maps the UMI-tagged reads to a reference genome using minimap2 with default settings for short genomice reads and converts the resulting SAM output into a sorted BAM file with samtools.
7. Remove Multi-Aligners: Removes any reads that map to multiple locations within the reference, retaining only those that have a single primary alignment and a MAPQ score greater than or equal to 30.
8. Deduplicate Reads: It uses UMICollapse dedup to identify and remove PCR duplicates from the sorted BAM file, using the UMI and alignment position of each read.
9. Finalize and Clean Up: It copies the final, deduplicated BAM file to a common output directory, verifies the copy's integrity with a checksum, indexes the file, and removes temporary intermediate files.
10. Reporting: Throughout the process, it tracks read counts and file sizes at each major step, printing a final summary table to the log. For BAM files, it reports the number of alignments. For fastq files, it reports the number of reads.

## BIDdetect

These are the same scripts as in Figure2/BIDdetect. Please refer there for documentation on usage for counting and aggregating counts

### Site Calling - In Vitro

This pipeline is optimized to run on Stanford's HPC Sherlock system. Make sure you adjust have R version 4.3 or higher installed.

Perform analysis for statistical difference and equivalence of sites. At least **two replicates** are required for the statistical analysis to work properly, and more is of course better. 

This script performs statistical analysis on processed Nano-BID-Amp data to identify significantly modified sites. It takes a counts table (typically generated by BIDdetect.sh) with deletion fractions as input and applies a mixed-effects model and equivalence testing to categorize each site.

Required Inputs:
- -i, --input: Path to the input data table containing per-replicate deletion fractions.
- -o, --outdir: Parent directory where a new, prefixed output folder will be created.
- -p, --prefix: A unique name for the analysis run, used for the output folder and all file names.

Optional Inputs:
- --cores: Number of CPU cores to use for parallel processing. (Default: -1, all available cores).
- --sesoi: The Smallest Effect Size of Interest (SESOI) for equivalence testing. (Default: 0.05).
- --color: The hex code for coloring "modified" sites in plots.
- --plot_all_sites: A flag to generate a large PDF containing a barplot for every individual site. Include the flag to generate this plot.

Key Outputs:
- Final Results Table (data_summary/<prefix>_modification_significance.tsv): Contains the final category ('Modified', 'Unmodified', 'Inconclusive') and all statistical results for every site.
- Data Subsets (data_raw/ and data_summary/): Includes raw and summary tables for all sites in each category, as well as for quartiles of the modified sites.
- Plots (plots/): A folder containing summary plots, including a boxplot of average deletion fractions, heatmaps of modified sites (with and without labels), and an optional PDF of all individual site plots.
- Log File (<prefix>_analysis_log.txt): A detailed log capturing all parameters, progress messages, and summary tables printed to the console.

Workflow:

1. Set Up: Before running, the script verifies that the R version is 4.3.0 or newer is available. If not, it stops with a helpful error message explaining how to load the correct modules on Sherlock. It parses all command-line arguments and creates the nested output directory structure.
2. Data Loading & Validation: It loads the input data table and verifies that all required columns are present. Required columns should be generated by the counting script if used properly. Required columns are:
    - chr, pos: extracted from BED file during counting, indicates the positon to analyze
    - treat: values should be in / input and BS / BID; column the statistics are performed between, **requires at least two data points in each condition** for the analysis to run, otherwise that site will be skipped
    - rep: indicates replicate number, need to have at least two for each treatment condition (BS, input)
    - delrate: deletion rate at a given sites
    - totalReads: total number of reads for a site
    - vector: (important) placeholder, indicates noPUS / PUS or WT / KD, more of a tool we'll use later but still required to have
3. Statistical Analysis (Parallelized): This is the core computational step.
    - It first pre-filters sites where all deletion fractions are below the SESOI, categorizing them as 'Unmodified'.
    - It then sets up a parallel processing backend using the number of cores specified by the --cores argument. More cores means faster processing.
    - For each site, it fits a Bayesian mixed-effects model (bglmer) to test for significant differences and performs an equivalence test (TOST) to check for a lack of meaningful change.
4. Site Categorization: It uses the results from the statistical tests to assign a final category ('Modified', 'Unmodified', or 'Inconclusive') to every site.
    - Modified: adjusted p-value indicates a significant difference between BS-treated and input
    - Unmodified: equivalence between sites, or every data point for a given site is less than SESOI
    - Inconclusive: Did not meet statistical significance or equivalence, likely due to high variance in underlying data.
5. Generate Outputs: The script saves the final summary table, creates and saves the data subsets (including a quartile analysis for modified sites), and generates the summary plots (box plot, heat maps).

Use example:

```bash
Rscript ~/Figure3/mpra_sites/modification_analysis.R \
  --input <BIDdetect_data.txt> \
  --outdir <output_dir> \
  --prefix <invitro> \
  --plot_all_sites
```

### Site Calling - In Cellulo

```Rscript incell_analysis.R```

This pipeline is optimized to run on Stanford's HPC Sherlock system. Make sure you adjust have R version 4.3 or higher installed.

Perform analysis for bisulfite-treatment-adjusted statistical difference of sites. This is for assessing the difference between WT and KD, two cell lines, or any other case where an additional factor is added on top of the bisulfite treatment. At least **two replicates** are required for the statistical analysis to work properly, and more is of course better. 

This script performs statistical analysis on processed Nano-BID-Amp data to identify significantly modified sites. It takes a counts table (typically generated by BIDdetect.sh) with deletion fractions as input and applies a mixed-effects model and equivalence testing to categorize each site. Analysis works as a combination of the in vitro site calling, described above, and the endogenous site calling, described in the Figure 2 README.

This is a rather hard-coded script. Please edit the "user configuration" section to access most of the same features as in the in vitro pipeline. If you desire to run a different type of comparison, I would recommend using this script as a template and building from there.

The script does the following three analyses:

1. **WT modification**: Identifies sites, using treatment analysis, that are modified, unmodified, or inconclusive in the WT state. This does NOT assign a PUS to the site.
    - It performs this separately for each cell type, and then for the cell types combined.
    - It saves a list of raw data, summarized data, modified sites, and unmodified sites.
    - Generates individual bar plots of every site (if selected) and summary box plot of all modified sites comparing BS to input values.

2. **PUS7 dependency**: Identifies sites, using factor analysis, that are modulated by PUS7 KD or OE.
    - It performs this separately for each cell type, and then for the cell types combined.
    - For each condition, it saves a list of raw data, summarized data, modified sites, and unmodified sites.
    - It takes that union of all PUS7 dependency sites and splits them into quartiles based on the average deletion deletion rate between all conditions.
    - As output, generates individual bar plots of every site (if selected) and summary box plot and heatmap of all PUS7-dependent sites.

2. **Cell Line Specificity**: Identifies sites, using factor analysis, that are modulated by cell type.
    - For each condition, it saves a list of raw data, summarized data, modified sites, and unmodified sites.
    - It takes finds the intersection of PUS7 dependent sites with cell type specific sites.
    - As output, generates individual bar plots of every site (if selected) and summary box plot of all cell type specific sites.

### Individual Sites

```Rscript indiv_sites.R```

This script generates plots 3C, 3E, S3E, and S3F, wherein individual sites from the MPRA are featured.

## In Vitro versus In Cellulo Comparisons

```Rscript invitro_incell_comp.R```

This script generates plots 3H and 3I, comparing deletion fractions between in cellulo and in vitro MPRA experiments.

## Threshold Parameter Sweep

```Rscript invitro_thresholds.R```

This script performs an exhaustive parameter sweep to determine how different filtering and statistical thresholds affect site calling using the factor analysis described above.
It iterates through thresholds of coverage, delta deletion fraction (smallest effect size of interest / SESOI), and p-value from a hard-coded input file.

It generates sequence logos and modified sites lists for each parameter combination considered. It produces a summary file comparing each combination to the baseline case.

```Rscript plots.R```

This script generates the Pareto front and metaplots displayed in S3B. Sequence logos and generated in the script prior.