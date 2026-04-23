# Figure 4: MPRA Sequence and Structure

This folder has code and files relevant to sequence and structure analysis, including the MPRA Nano-SHAPE-Amp pipeline.

## Sequence Analysis

```Rscript seqlogos.R```

Creates sequence logos for MPRA sites for background (4A), PUS7-dependent (4A), and quartile (4I) sets. 
Paths are hard-coded to input data within the script.

```Rscript UNUAR_enrichment.R```

The script reads in the hard-coded files with information about the MPRA background and PUS7-dependent sites. It filters for UNUAR motifs only. It then prints a composition summary and Fisher's exact test report of PUS7-dependent and background UNUAR sites, in vitro and in cellulo. Finally, it generated the stacked bar chart in 4B.

## RNAfold

To gather information of pairing probabilities and folding energies, Vienna RNA 2.5.1 is used. 
The following scripts are optimized to use on Sherlock, Stanford's HPC.

Use as follows:

```bash run_RNAfold.sh <fasta_file> <output_dir>```

This runs RNAfold on every sequence in the given fasta and outputs all files to the specified directory.

To concatenate these to a single csv file more suitable for downstream analysis, use the following:

```
Rscript ~/PUS7regulation2026/Figure4/RNAfold/extract_pairing_prob.R \
    -i <input directory> \
    -o <output csv> \
    -l <sequence length>
```

The above script extracts the pairing probabilites from the "ubox" lines in all dp.ps files output by RNAfold. This represents the total of all pairing probabilities for a given sequence, NOT the pairing probabilites of the MFE structures. 
The sequence length argument dictates how many columns will be created and should be equivalent to the length of the sequence.

```
Rscript ~/PUS7regulation2026/Figure4/RNAfold/fold_summary.R \
    -i <input directory> \
    -o <output csv>
```

The above script extracts information from all of the .fold files produced by RNAfold. 
This includes the sequence and dot-bracket notation for the MFE, MEA, and centroid structures, as well as the free energy (delta G) for those respective structure.

## Basecalling

``` 
bash run_basecaller_sup_pipeline.sh \
  -o <output_directory> \
  -s <script_directory> \
  -p <pod5_directory> \
  -b <number_of_batches> \
  -m <max_number_of_jobs> 
```

Runs two scripts in sequence: basecaller_sup_array.sbatch and demuxer_merge.sbatch

This defaults to running super-high-accuracy DNA basecalling with **dorado 1.1.0**, as set in the basecaller_sup_array.sbatch script.

To specifiy the GPU to basecall on, modify the SBATCH part of the basecaller_sup_array.sbatch. Default is as so:

```#SBATCH -C '[GPU_GEN:HPR|GPU_GEN:LOV|GPU_GEN:VLT]' ```

Which allows basecalling on any GPU in the Hopper, Lovelace, or Voltage generation. 

To change any parameters around basecalling (kit, basecalling model, modified bases, etc), update the dorado command in basecaller_sup_array.sbatch. The current command is:

```
dorado basecaller sup "$INPUT_DIR" \
    --kit-name EXP-PBC096 \
    --no-trim > "$BAM_OUTPUT_DIR/basecalled.bam"
``` 

Basecalling outputs a bam file, which is then demultiplexed by barcode into a single bam file per barcode in demuxer_merge.sbatch. Make changes to the demultiplexer with this command:

```
dorado demux --output-dir "$FINAL_RESULTS_DIR" \
    --no-trim \
    --no-classify \
    "$MERGED_BAM_FILE"
```

### Basecalling Job Fails

Sometimes, a GPU can go kaput halfway through your job. Never fear, you do not have to re-run the full script from above. Instead, you need to look at resume_failed_pipeline.sh.

Identify the tasks numbers that failed and update the user configuration settings accordingly. This will re-submit basecalling jobs only for the batches that intitially failed. It will then carry through merging both the basecalled batches that completed on the first attempt and those that did not. You will end up at the same spot: one bam file in a given directory.

## Mapping Reads: Nanopore Pipeline

This pipeline is optimized to run on Stanford's HPC Sherlock system. Make sure you adjust have the following installed and adjust package loading appropriately:

- samtools 1.16.1
- cutadapt 1.18 (requires Python 3.6)
- java 11
- bowtie2 2.3.4.1
- UMICollapse: https://github.com/Daniel-Liu-c0deb0t/UMICollapse 


Start from a bam file labeled with a barcode from Nanopore basecalling outputs,
ideally with super-high-accuracy basecalling.

The pipeline (```SHAPE_trim_map_dedup_mpra.sh```) goes through the following steps:
1. Convert bam to fastq, renaming from barcode in the process.
2. Trim Nanopore adapters from the fastq in sense and antisense directions. 
3. Reverse complement the antisense to unify orientation of reads.
4. Extract UMIs.
5. Trim the MPRA adapters in the sense direction.
6. Map reads with bowtie2.
7. Convert to sorted bam & index.
8. Deduplicate reads.
9. Convert to fastq for input into shapemapper

The following variables are hardcoded into process_single_file.sh:
- SENSE_ADAPTER_5PRIME="TTTCTGTTGGTGCTGATATTGCG"
- SENSE_ADAPTER_3PRIME="GAAGATAGAGCGACAGGCAAGT"
- ANTISENSE_ADAPTER_5PRIME="ACTTGCCTGTCGCTCTATCTTC"
- ANTISENSE_ADAPTER_3PRIME="CGCAATATCAGCACCAACAGAAA"
- POOL_SENSE_ADAPTER_5PRIME="GACGCTCTTCCGATCT"
- POOL_SENSE_ADAPTER_3PRIME="CACTCGGGCACCAAGGAC"
- UMI_PATTERN="NNNNNNNNNN"

To use, you need to submit it as a slurm array using ```submit_slurm_array.sh```:

Example command:
```
bash ~/PUS7regulation2026/Figure4/SHAPE_prep/submit_slurm_array.sh \
    --mail-user you@institution.edu \
    --script-path ~/PUS7regulation2026/Figure4/SHAPE_prep/SHAPE_trim_map_dedup_mpra.sh \
    --map-file <barcodes.txt> \
    --input-dir <basecalling output directory> \
    --output-dir <output dir> \
    --bowtie-index <bowtie index>
```

The barcodes file should be formatted as follows:

```
HEK293T_SHAPE_Rep1:barcode75 
HEK293T_DMSO_Rep1:barcode76 
```

This will submit an individual job for each bam file in the bam directory.
The final deduplicated fastq will be located in a directory inside the output directory.

## Mapping Reads: Illumina Pipeline

The Illumina data starts in a slightly different place, namely pair-ended reads, than the Nanopore data, so it requires a slightly different pipeline.

The same packages are still required:

- samtools 1.16.1
- cutadapt 1.18 (requires Python 3.6)
- java 11
- bowtie2 2.3.4.1
- UMICollapse: https://github.com/Daniel-Liu-c0deb0t/UMICollapse 

The pipeline (```process_single_file_illumina_downsample.sh```) goes through the following steps:
1. Starts with pair-ended .fastq.gz files.
2. Extract UMIs from the start of R2 and appends it to R1 headers.
3. Trims 3' adapters from R1.
4. Downsamples to ~20M reads to approximately match the read counts in the Nanopore sample.
5. Map reads with bowtie2.
6. Convert to sorted bam & index.
7. Deduplicate reads.
8. Convert to fastq for input into shapemapper

The following variables are hardcoded into process_single_file.sh:
- SENSE_ADAPTER_5PRIME="TTTCTGTTGGTGCTGATATTGCG"
- SENSE_ADAPTER_3PRIME="GAAGATAGAGCGACAGGCAAGT"
- ANTISENSE_ADAPTER_5PRIME="ACTTGCCTGTCGCTCTATCTTC"
- ANTISENSE_ADAPTER_3PRIME="CGCAATATCAGCACCAACAGAAA"
- POOL_SENSE_ADAPTER_5PRIME="GACGCTCTTCCGATCT"
- POOL_SENSE_ADAPTER_3PRIME="CACTCGGGCACCAAGGAC"
- UMI_PATTERN="NNNNNNNNNN"

To use, you need to submit it as a slurm array using ```submit_slurm_array_illumina.sh```:

Example command:
```
bash ~/PUS7regulation2026/Figure4/SHAPE_prep/submit_slurm_array_illumina.sh \
    --mail-user you@institution.edu \
    --script-path ~/PUS7regulation2026/Figure4/SHAPE_prep/process_single_file_illumina_downsample.sh" \
    --sample-map <map_file.tsv> \
    --output-dir <output dir> \
    --ref-index <bowtie index> 
    
```

The sample map file should be formatted as follows:

```
No_PUS_DMSO No_PUS_DMSO_S4_L001_R1_001.fastq.gz No_PUS_DMSO_S4_L001_R2_001.fastq.gz
No_PUS_SHAPE    No_PUS_SHAPE_S3_L001_R1_001.fastq.gz    No_PUS_SHAPE_S3_L001_R2_001.fastq.gz
```

This will submit an individual job for each pair of fastq files in the map file.
The final deduplicated fastq will be located in a directory inside the output directory.

## SHAPE-Mapper

This pipeline (```shapemapper_pipeline.sbatch```) goes through the following steps:

1. Runs shapemapper2-2.3 paired fastq files (DMSO and SHAPE)
2. Annotates RNAs with poor quality scores
3. Produces correlation plots of the replicates
4. Averages SHAPE reactivities across replicates
5. Performs SHAPE-informed RNA fold
6. Extracts pairing probabilities and structure feature summaries

This pipeline works the same for reads of Nanopore or Illumina origin, once they have been processed through the above pipelines.

One limitations of shapemapper is the number of sequences it can process at a given time. It will crash if your reference FASTA has more than ~100 sequences (represented by > header lines). 
To get around this, you should "chunk" your fasta into smaller lists. 
This script will run the chunks separately and then concatenate the results.

Shapemapper log output includes a list of "poor quality" RNAs, indicating an issue identified during processing.
This typically results from inadequate read coverage for a given sequence. 
Low coverage can skew reactivities, so these should be excluded from downstream analysis, which is what the annotation section does.

In the presence of multiple replicates, inverse variance weighting averaging is used to create one set of shape reactivites to go into RNAfold.

After running SHAPE-informed RNAfold, pairing probabilities and folding information is extracted using the same scripts as in the RNAfold folder.

This script is optimized to work on Sherlock, Stanford's HPC. The following packages are required for this to run:
- Vienna RNA 2.5.1
- shapemapper2-2.3 https://github.com/Weeks-UNC/shapemapper2

Example command to run a single condition with two replicates:

```
sbatch --array=1-1 \
    --mail-user=you@institution.edu \
    ~/PUS7regulation2026/Figure4/SHAPE_mapper/shapemapper_pipeline.sbatch \
        --sample-map <sample_map.tsv> \
        --output-dir <sample output dir> \
        --ref-fasta-dir <reference_fasta_chunks> \
        --ref-rnafold-fasta <reference_fasta> \
        --num-samples 1
```

The sample map should look as follows:
```
SampleID    Untreated   Treated
NoPus_Rep1	~/Rep1/final_deduplicated_fastq/NoPUS_DMSO.final.fastq	~/Rep1/final_deduplicated_fastq/NoPUS_SHAPE.final.fastq
NoPus_Rep2	~/Rep2/final_deduplicated_fastq/NoPUS_DMSO.final.fastq	~/Rep2/final_deduplicated_fastq/NoPUS_SHAPE.final.fastq
```

## Structure Visualization

These script automates the batch generation of 2D RNA MEA structure vector-based SVG diagrams using VARNA, without or without mapping SHAPE reactivity data to nucleotide colors (low, medium, high). It can also draw a custom bounding box and apply special formatting to highlight a specific nucleotide position of interest, such as the target uridine.

These scripts require python 3.9. 

Activate the virtual environment and use as follows:

```
source ~/PUS7regulation2026/Figure4/visualization/varna-env/bin/activate

python3 ~/PUS7regulation2026/Figure4/visualization/run_varna_highlight_SHAPE.py \
  --input-csv <fold summary csv> \
  --shape-dir <shape directory> \
  --highlight-csv <highlight csv> \
  --output-dir <output dir> 

python3 ~/PUS7regulation2026/Figure4/visualization/run_varna_highlight.py \
  --input-csv <fold summary csv> \
  --highlight-pos <int>> \
  --output-dir <output dir> 
```

The fold summary csv should contain the sequence and dot-bracket notation of the MEA structure, and which can be produced by ```fold_summary.R```.

The highlight csv should contain the site name and highlight position (1-indexed) (optional arguments to specify) for the residue that should be highlighted in pink for the SHAPE script. For the non-SHAPE script, the highlight is done at a single consistent position.

## Illumina and Nanopore Comparison

The source data is from the same SHAPE library made on the same day. Samples were split between Illumina and Nanopore sequencing.

To compare the SHAPE reactivies, first a correlation plot of the normalized reactivities was generated (S4C). 
Next, the Illumina reactivites were binarized as reacitve or non-reactive based on three thresholds: 0.3 (4D), 0.5, and 0.7 (S4D).
ROC curves are generated with an AUC value.

The files that go into this script should be the "profile_annotated.txt" files produced at the end of the shapemapper pipeline.
The script can be run as follows:

```
Rscript ~/PUS7regulation2026/Figure4/Illumina_Nanopore/comparison_Nano_Illum.R \
    --nanopore <profile_annotated.txt> \
    --illumina <profile_annotated.txt> \
    --out-prefix <output_dir/sample>
```

All plots are saved at the specified output location.

## Pairing Probability Metaplots & Analysis

Information about minimum free energy were extracted from RNAfold analysis as described above (section: RNAfold).

To create MFE plots, run ```Rscript MFE/boxplots.R```
File paths are hard-coded into the script.
This script generates boxplots comparing the minimum free energy of PUS7-dependent targets versus non-dependent targets.
Minimum free energy is from in vitro SHAPE-informed folding (4E) or RNAfold only (S4A) predictions.


Pairing probabilities were extracted from RNAfold analyses as described above (section: RNAfold).

To create metaplots, run ```Rscript pairing_prob/metaplots.R```.
File paths are hard-coded into the script.
The script provides the function to align pairing probabities around a center position, but this has already been done in the hard-coded files.
Next, the script generates multi-line metaplots for direct RNA (S4B, 4G, S4F), in vitro MPRA (S4B, 4G), in cellulo MPRA (S4B, 4G, S4F), 
in vitro and in cellulo quartiles (S4I) with RNA fold only (S4B), in vitro SHAPE (4G), and in cell SHAPE (S4F) pairing probabiilties.


To create target uridine pairing probabiity plots (4J), run ```Rscript pairing_prob/target_uridine.R```.
File paths are hard-coded into the script. 
It generates boxplots comparing the target uridine pairing probability of the in vitro SHAPE-informed structures across quartiles for both in vitro and in cellulo conditions.