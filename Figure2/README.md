# Figure 2: endogenous Nano-BID-Amp

This folder has code and files relevant to the endogenous Nano-BID-Amp pipeline and analysis.

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
bash ~/file_prep/submit_file_prep.sh \
    -u you@institution.edu -t 1:00:00 \
    -s ~/file_prep/trim_map_dedup.sh \
    -a <barcodes.txt> \
    -i <input_fastq> \
    -o <output_directory> \
    -r <reference_fasta>
```

### trimp_map_dedup.sh

This script processes a single Nano-BID-Amp sample as part of a SLURM array job. This script is appropriate for amplicons where there is only one set of adapters, the Nanopore adapters, to trim from the reads.

Its steps are:
1. Setup and Initialization: It parses command-line arguments, uses the SLURM array task ID to identify a specific sample from a map file, sets up output directories, and redirects all output to a log file.
2. Filter Reads: It uses cutadapt with linked adapters to process the compressed FASTQ file twice, creating two separate files: one containing "sense" reads and another with "antisense" reads, based on the Nanopore adapter sequences.
3. Unify Read Orientation: It reverse-complements the "antisense" reads and then concatenates them with the "sense" reads. This creates a single FASTQ file where all reads are oriented in the sense direction.
4. Extract UMIs: It uses umi_tools extract to pull 10-nucleotide Unique Molecular Identifiers (UMIs) from the 3' end of each read, embedding the UMI into the read's header for later use.
5. Align Reads: It maps the UMI-tagged reads to a reference genome using minimap2 with default settings for long read splice alignment and converts the resulting SAM output into a sorted BAM file with samtools.
6. Filtering: It filters the mapped file to only contain primary alignments with a MAPQ > 30.
7. Deduplicate Reads: It uses UMICollapse dedup to identify and remove PCR duplicates from the sorted BAM file, using the UMI and alignment position of each read.
8. Finalize and Clean Up: It copies the final, deduplicated BAM file to a common output directory, verifies the copy's integrity with a checksum, indexes the file, and removes temporary intermediate files.
9. Reporting: Throughout the process, it tracks read counts and file sizes at each major step, printing a final summary table to the log. For BAM files, it reports the number of alignments. For fastq files, it reports the number of reads.

## BIDdetect

Aggregates the nucleotide counts for all samples. All scripts are intended to work on Sherlock, Stanford's HPC.

## BIDdetect.sh

This script automates the process of calculating and aggregating nucleotide counts from multiple bam files.

Required Inputs:
- -b, --bam_dir: A directory containing one or more sorted and indexed BAM files (.bam).
- -o, --output_dir: A path to a directory where all results and logs will be stored.
- -r, --ref_fasta: The path to the reference genome FASTA file.
- -e, --bed_file: The path to a BED file defining the genomic regions of interest. Requires 6 columns.

Optional Input:
- -n, --col_names: A string defining column names for the final formatting step (Default: 'celltype_vector_rep_treat').

Key Outputs:

- Final Data File (${DEST_DIR}/BIDdetect_data.txt): The main result of the pipeline, containing the fully processed and formatted count data after all steps. This should be taken into further analyses.
- Aggregated Counts File (${DEST_DIR}/BIDdetect_counts.txt): A master table containing the raw, combined counts from all input BAM files before the final formatting step.
- Log File (${DEST_DIR}/logs/BIDdetect_...log): A timestamped log file that captures all screen output for debugging and record-keeping.
- Intermediate Files (${DEST_DIR}/intermediate_counts/): A directory holding the temporary count files generated for each individual sample before they are aggregated.

Workflow:

1. Argument Parsing & Validation: makes sure all required variables are present and valid
2. Environment Setup: loads appropriate modules, sets up directories, creates logging system
3. Main Processing Loop (Per-Sample Analysis): iterates through every (indexed) .bam file in the input directory
    1. Extracts a clean sample name from the filename. 
    2. Executes the bam_counts_fast.R script that performs a pileup on the BAM file at the positions defined in the BED file to generate nucleotide, deletion, and insertion counts. Calculates the deletion rate (deletions / total reads) for each site. 
    3. Saves the per-sample results to a temporary file in the intermediate_counts directory. 
    4. Appends the results from the temporary file to a single master file (BIDdetect_counts.txt), adding a new column containing the sample name.
4. Variable Extraction from File Name: calls sample_name.R to split the file name into 4 different columns, as defined in the --col_names argument.
    - default is celltype_vector_rep_treat, which is the recommended naming of bam files to ensure integration with analysis scripts
5. Final Processing & Cleanup: saves final data to BIDdetect_data.txt, outputs locations of final files

Usage example:

```
bash $HOME/BIDamplicon/pipeline/BIDdetect.sh \
  --bam_dir <prep_output/deduplicated_bam> \
  --output_dir <output_dir> \
  --ref_fasta <reference_fasta> \
  --bed_file <reference_bed>
  ```

## Site Calling

```Rscript analysis_endo.R```

Generic script used for preliminary analysis and plotting of endogenous Nano-BID-Amp. 
Will not produce the exact plot in the paper, but should get close enough. 
Hard-coded with input file, which contains generic BIDdetect output file name.

The script does the following:
1. Creates a bar plot for each WT site with coverage, compared between BS and input.
2. Runs standard "treatment analysis" to determine sites which are statistically significant, indicating pseudouridine presence, in WT.
3. Create a bar plot for each WT and KD sites with coverage, plotting only the BS points.
4. Performs "factor analysis" to determine which sites have a bisulfite-treatment-adjusted statistically significant difference between WT and KD.
5. Performs "factor analysis" between two cell types to determine which sites have a bisulfite-treatment-adjusted statistically significant difference.

## Deletion Fraction Comparisons

```Rscript delfrac_comparison.R```

This script aggregates data from endogenous Nano-BID-Amp, BID WT, and PRAISE WT experiments.
It compares the deletion fractions between sites with coverage and generates plot S2D.