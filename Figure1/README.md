# Figure 1: direct RNA sequencing

This folder has code and files relevant to the direct RNA sequencing pipeline and analysis.

## Basecalling
The version of dorado used for basecalling can influence mismatch rates. We found that older models produced higher mismatch rates than newer models. All direct RNA sequencing in this manuscript was basecalled with **dorado 0.7.3**.

These scripts are optimized to run on Stanford's HPC, Sherlock, as multiple high quality GPUs are required for basecalling to finish in a timely manner. Scripts can be adapted for use on other clusters with the help of your favorite AI agent.

The pipeline can be called as follows:

```bash
bash run_basecaller_rna_pipeline.sh \
  -o <output_directory> \
  -s <script_directory> \
  -p <pod5_directory> \
  -b <number_of_batches> \
  -m <max_number_of_jobs>
```
Runs two scripts in sequence: basecaller_rna_array.sbatch and merge.sbatch

This defaults to running high-accuracy basecalling with dorado/0.7.3, as set in the basecaller_rna_array.sbatch script.

To specifiy the GPU to basecall on, modify the SBATCH part of the basecaller_sup_array.sbatch. Default is as so:

```#SBATCH -C '[GPU_GEN:HPR|GPU_GEN:LOV|GPU_GEN:VLT]' ```



To change any parameters around basecalling (kit, basecalling model, modified bases, etc), update the dorado command in basecaller_rna_array.sbatch. The current command is:

```bash
dorado basecaller hac "$INPUT_DIR" > "$BAM_OUTPUT_DIR/basecalled.bam"
```

Basecalling outputs one bam file per batch, which is then merged with samtools into a single bam file merge.sbatch. 

### Basecalling Job Fails

Sometimes, a GPU can go kaput halfway through your job. Never fear, you do not have to re-run the full script from above. Instead, you need to look at resume_failed_pipeline.sh.

Identify the tasks numbers that failed and update the user configuration settings accordingly. This will re-submit basecalling jobs only for the batches that intitially failed. It will then carry through merging both the basecalled batches that completed on the first attempt and those that did not. You will end up at the same spot: one bam file in a given directory.

## ModDetect

ModDetect was created based off of work done in the Rouhanifard lab:
- github.com/RouhanifardLab/PsiNanopore
- https://www.nature.com/articles/s41467-023-35858-w

ModDetect is a versatile RNA modification detection pipeline that can be used to identify candidate mod sites in a test sample when compared to a control sample. The script's default functionality searches for U-to-C mismatch signatures, but can be modified using any combination of flags listed below to look for deletion signatures or to perform a mod-agnostic search. Two samples are required: test (WT) and control (KD).

### Required Inputs: 

1. Bam file for test sample (along with bam index file that ends with file extension .bai)
2. Bam file for control sample (along with bam index file that ends with file extension .bai)
3. Fasta file for reference genome
4. BED file associated with the reference genome which will be used to define the genomic search space intervals (10 columns)

### Options: 
```
-f	Path to bam file for test / higher sample (REQUIRED)
-g	Path to bam file for control / lower sample (REQUIRED)
-k	Path to kmer-dependent error file (baseline error rate for each k-mer combination)
-r	Path to reference genome fasta file (REQUIRED)
-b	Path to BED file associated with the reference genome (REQUIRED)
-o	Path to output file location
-d	Examines deletion signatures instead of U-to-C mismatches (useful for bisulfite-treated samples)
-a	Examines total mismatches at all positions instead of just those with a 'T' in the reference genome (enables mod-agnostic search)
-x	Read depth cutoff (positions that do not have sufficient coverage in both the test and control samples will be filtered out of the final output)
-y	Mismatch percentage delta cutoff (minimum acceptable difference in percentage mismatch between control and test samples)
-m	Sites with a p-value this low or lower will be outputted
```

Example Usage: 
```
Rscript ModDetect.R \
    -f <WT.bam> \
    -g <KD.bam> \
    -r <reference.fasta> \
    -o <output.csv> \
    -b <reference.bed> \
    -x <read coverage threshold> \
    -y <mismatch threshold>
```

Notes:

- The "R" folder must also be present in the same directory as the ModDetect.R script in order for the necessary packages to be installed. The script will take a bit longer when run for the first time to install these packages.
- The p-value is not an accurate test of statistical significance and should not be used as a selection factor in identifying sites.

## UNUAR Standard

### Information about standard
Seqeunce: GUUGAUUΨAACAUCGGUUGAUUΨAGCAUCGGUUGAUAΨAACAUCGGUUGAUAΨAGCAUCGGUUGAUGΨAACAUCGGUUGAUGΨAGCAUCGGUUGAUCΨAACAUCGGUUGAUCΨAGCAUCGAAAAAAAAAA

Pseudouridine positions: 8, 23, 38, 53, 68, 83, 98, 113

Following basecalling and demultiplexing (described above / in methods), information about mismatch rate was collected using the script ```bam_counts.R``` with the ```-all TRUE``` flag to capture all reference bases.

Given multiple bam files for the standards, this was folded into the wrapper script ```count_folder.sh```. The paths to the input bam folder, output directory, reference FASTA, and reference BED file are all hard-coded into the wrapper, which can be run as so:

```bash count_folder.sh```

Bam files are available on GEO at GSE314741.

The produces a single .csv file, UNUAR_counts.txt, which is used for generation of calibration curves.

### Calibration Curves 

```Rscript UNUARstandard/UNUAR_calibration_curves.R```

This script does the following:

1. Reads in the hard-coded file path to the UNUAR count sheet.
2. Calculates U-to-C mismatch rate as C/(U+C)
3. Creates plot S1C showing U-to-C mismatch rate across the whole sequence.
4. Generates calibration curves with the U-to-C mismatch rates with no pseudouridine subtracted.
5. Generates plot S1D showing calibration curve for each motif.
6. Saves formulas to correct UNUAR motif mismatch rate to absolute stoichiometry as derivation_formulas_subtract.csv

### Stoichiometry Correction

```Rscript Derivation_Subtract.R```

This script reads in hard-coded file for aggregated WT and KD mismatch rates does the following:

1. Averages mismatch rates across WT and KD replicates.
2. Calculate corrected stoichiometries using formulas derived above.
3. Makes boxplots of corrected and uncorrected values for WT and KD.
4. Makes heatmap of corrected and uncorrected values for WT and KD.

KO analysis was done separately using the same workflow.

## Borderline Site Identification

```Rscript borderline_dRNA.R```

This scripts read in the background "U" sites with adequate coverage (> 10 reads). It defines borderline sites as those NOT called with the KD, conforming to a UNUAR motif, and with a WT mismatch rate >= 15 in at least two of three biological replicates.

## Motif Enrichment Analysis

The motifs displayed in 1E were generated with the online STREME tool. This section speaks to the UNUAR motif enrichment in 1F.

```Rscript UNUAR_enrichment.R```

The script read in the hard-coded FASTAs with 5-mer motifs of PUS7-dependent and background sets. It filters for UNUAR motifs only. It then prints a composition summary and Fisher's exact test report of PUS7-dependent and background UNUAR sites. Finally, it generated the stacked bar chart in 1F.

## PUS7 KO

Additional analyses for the PUS7 KO.

```Rscript KO_analysis.R```

Generates plot S1I and S1J using PUS7 KO mismatch rates, as well as additional plots used in preliminary analyses.