#!/usr/bin/env Rscript

library(dplyr)
library(stringr)
library(readr)
library(ggplot2)
library(readxl)
library(scales)
library(tidyr)
library(purrr)


# ---- Directory ----
active_dir <- "~/PUS7regulation2026/Figure1/borderline"

# ---- Datasets ----
dRNA_background <- read.csv(paste0(active_dir, "/PUS7_10reads_0percent_finalized.csv"))

dRNA_sites <- read.csv("~/PUS7regulation2026/Figure1/UNUARstandard/PUS7_10reads_15percent_motif.csv")

# ---- Identify sites with 2 / 3 ----
dRNA_unuar_notcalled <- dRNA_background %>%
    filter(str_detect(motif, "T.TA[AG]")) %>%
    anti_join(dRNA_sites, by = c("chr" = "chr", "pos" = "pos"))

dRNA_borderline <- dRNA_unuar_notcalled %>%
    rowwise() %>%
    filter(sum(c(mm.perc.rep2, mm.perc.rep3, mm.perc.rep4) > 15, na.rm = TRUE) >= 2) %>%
    ungroup() %>%
    select(
    chr, pos, gene_name, gene_strand, motif, 
    # Renaming mm.perc to WT_mm
    WT_mm_rep2 = mm.perc.rep2, 
    WT_mm_rep3 = mm.perc.rep3, 
    WT_mm_rep4 = mm.perc.rep4,
    # Renaming ctrl.err to KD_mm
    KD_mm_rep2 = ctrl.err.rep2, 
    KD_mm_rep3 = ctrl.err.rep3, 
    KD_mm_rep4 = ctrl.err.rep4,
    overlap
  )

write.csv(dRNA_borderline, paste0(active_dir, "/dRNA_borderline_sites.csv"), row.names = FALSE)
