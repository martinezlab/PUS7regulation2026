#!/usr/bin/env Rscript

library(Biostrings)
library(ggseqlogo)
library(ggplot2)
library(dplyr)

create_sequence_logo <- function(sequences, title, subtitle = NULL, center_position) {
  if (length(center_position) == 1) {
    center_position <- rep(center_position, length(sequences))
  }
  # Define the center position and extract the 5mer
  center_position <- center_position
  motif_start <- center_position - 2
  motif_end <- center_position + 2
  extracted_motifs <- subseq(sequences, start = motif_start, end = motif_end)
  
  # Convert to character vector and replace T with U
  motif_chars <- chartr("T", "U", as.character(extracted_motifs))
  
  # Create position probability matrix
  ppm <- consensusMatrix(motif_chars, as.prob = TRUE)
  ppm <- ppm[c("A", "C", "G", "U"),]  # Keep only A, C, G, U rows and reorder
  
  # Create the plot
  p <- ggseqlogo(ppm, method = "prob") + 
    labs(title = title,
         subtitle = subtitle,
         y = "Probability") +
    scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, by = 0.2)) +
    theme_logo() +  # Keep the ggseqlogo default theme as base
    theme(
      axis.text.x = element_text(size = 9),
      axis.text.y = element_text(size = 9),
      axis.title.x = element_blank(),
      axis.title.y = element_text(size = 12),
      plot.title = element_text(size = 14, hjust = 0.5),
      plot.subtitle = element_text(size = 12, hjust = 0.5), 
      panel.background = element_blank(),
      panel.border = element_blank(),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      plot.margin = margin(10, 10, 10, 10, "pt")
    )
  # Save the plot as PDF
    filename <- paste0(title, "_seqlogo")
    path <- file.path(plot_dir, filename)
    ggsave(paste0(path, ".pdf"), p, width = 2.5, height = 1.6, units = "in")
    ggsave(paste0(path, ".png"), p, width = 2.5, height = 1.6, units = "in")
  # Return the plot object
  return(p)
}

invitro_full <- read.csv("invitroBID_MPRA_summary.csv", header = TRUE)

seq_ivPUS7 <- invitro_full %>% filter(PUS7_dependent = "Yes") %>% filter(!(is.na(psipos_motif)))
seq_NOTivPUS7 <- invitro_full %>% filter(PUS7_dependent = "No") %>% filter(!(is.na(psipos_motif)))
create_sequence_logo(sequences = seq_ivPUS7$psipos_motif, 
                        title = "invitro_PUS7", 
                        subtitle = paste("Based on", nrow(seq_ivPUS7), "sites"),
                        center_position = 3)
create_sequence_logo(sequences = seq_NOTivPUS7$psipos_motif, 
                        title = "invitro_background", 
                        subtitle = paste("Based on", nrow(seq_NOTivPUS7), "sites"),
                        center_position = 3)

incell_full <- read.csv("incellBID_WTsummary.csv", header = TRUE)

seq_icPUS7 <- incell_full %>% filter(PUS7_dependent = "Yes") %>% filter(!(is.na(psipos_motif)))
seq_NOTicPUS7 <- incell_full %>% filter(PUS7_dependent = "No") %>% filter(!(is.na(psipos_motif)))
create_sequence_logo(sequences = seq_icPUS7$psipos_motif, 
                        title = "incell_PUS7", 
                        subtitle = paste("Based on", nrow(seq_icPUS7), "sites"),
                        center_position = 3)
create_sequence_logo(sequences = seq_NOTicPUS7$psipos_motif, 
                        title = "incell_background", 
                        subtitle = paste("Based on", nrow(seq_NOTicPUS7), "sites"),
                        center_position = 3)

ivPUS7_Q1 <- read_tsv("invitro_quartiles/invitro_modified_sites_Q1_summary.tsv", col_names = TRUE)
ivPUS7_Q2 <- read_tsv("invitro_quartiles/invitro_modified_sites_Q2_summary.tsv", col_names = TRUE)
ivPUS7_Q3 <- read_tsv("invitro_quartiles/invitro_modified_sites_Q3_summary.tsv", col_names = TRUE)
ivPUS7_Q4 <- read_tsv("invitro_quartiles/invitro_modified_sites_Q4_summary.tsv", col_names = TRUE)

seq_ivPUS7_Q1 <- invitro_full %>% inner_join(ivPUS7_Q1, by = ("site" = "chr")) %>% filter(!(is.na(psipos_motif)))
seq_ivPUS7_Q2 <- invitro_full %>% inner_join(ivPUS7_Q2, by = ("site" = "chr")) %>% filter(!(is.na(psipos_motif)))
seq_ivPUS7_Q3 <- invitro_full %>% inner_join(ivPUS7_Q3, by = ("site" = "chr")) %>% filter(!(is.na(psipos_motif)))
seq_ivPUS7_Q4 <- invitro_full %>% inner_join(ivPUS7_Q4, by = ("site" = "chr")) %>% filter(!(is.na(psipos_motif)))

create_sequence_logo(sequences = seq_ivPUS7_Q1$psipos_motif, 
                        title = "invitro_PUS7_Q1", 
                        subtitle = paste("Based on", nrow(seq_ivPUS7_Q1), "sites"),
                        center_position = 3)
create_sequence_logo(sequences = seq_ivPUS7_Q2$psipos_motif, 
                        title = "invitro_PUS7_Q2", 
                        subtitle = paste("Based on", nrow(seq_ivPUS7_Q2), "sites"),
                        center_position = 3)
create_sequence_logo(sequences = seq_ivPUS7_Q3$psipos_motif, 
                        title = "invitro_PUS7_Q3", 
                        subtitle = paste("Based on", nrow(seq_ivPUS7_Q3), "sites"),
                        center_position = 3)
create_sequence_logo(sequences = seq_ivPUS7_Q4$psipos_motif, 
                        title = "invitro_PUS7_Q4", 
                        subtitle = paste("Based on", nrow(seq_ivPUS7_Q4), "sites"),
                        center_position = 3)

icPUS7_BothQuartiles <- read_tsv("incell_quartiles/PUS7_dep_union_by_Both_WT_quartiles.tsv", col_names = TRUE)

seq_icPUS7_Q1 <- incell_full %>% inner_join(icPUS7_BothQuartiles, by = ("chr" = "chr")) %>% filter(quartile == 1) %>% filter(!(is.na(psipos_motif)))
seq_icPUS7_Q2 <- incell_full %>% inner_join(icPUS7_BothQuartiles, by = ("chr" = "chr")) %>% filter(quartile == 2) %>% filter(!(is.na(psipos_motif)))
seq_icPUS7_Q3 <- incell_full %>% inner_join(icPUS7_BothQuartiles, by = ("chr" = "chr")) %>% filter(quartile == 3) %>% filter(!(is.na(psipos_motif)))
seq_icPUS7_Q4 <- incell_full %>% inner_join(icPUS7_BothQuartiles, by = ("chr" = "chr")) %>% filter(quartile == 4) %>% filter(!(is.na(psipos_motif)))
 
create_sequence_logo(sequences = seq_icPUS7_Q1$psipos_motif, 
                        title = "incell_PUS7_Q1", 
                        subtitle = paste("Based on", nrow(seq_icPUS7_Q1), "sites"),
                        center_position = 3)
create_sequence_logo(sequences = seq_icPUS7_Q2$psipos_motif, 
                        title = "incell_PUS7_Q2", 
                        subtitle = paste("Based on", nrow(seq_icPUS7_Q2), "sites"),
                        center_position = 3)
create_sequence_logo(sequences = seq_icPUS7_Q3$psipos_motif, 
                        title = "incell_PUS7_Q3", 
                        subtitle = paste("Based on", nrow(seq_icPUS7_Q3), "sites"),
                        center_position = 3)
create_sequence_logo(sequences = seq_icPUS7_Q4$psipos_motif, 
                        title = "incell_PUS7_Q4", 
                        subtitle = paste("Based on", nrow(seq_icPUS7_Q4), "sites"),
                        center_position = 3)