#!/usr/bin/env Rscript

library(dplyr)
library(stringr)
library(readr)
library(ggplot2)
library(tidyr)
library(purrr)
library(broom)
library(Biostrings)

# ---- Directory ----
active_dir <- "~/PUS7regulation2026/Figure1/motif"

# 5mer
dRNA_UNUAR_seq <- readDNAStringSet(paste0(active_dir, "/dRNA_UNUAR_motifs.fa"))
dRNA_background_seq <- readDNAStringSet(paste0(active_dir, "/dRNA_background.fasta"))

# --- PART 1: Convert to Dataframe ---
# Convert the DNAStringSet objects to standard dataframes
unuar_df <- tibble(
  sequence = as.character(dRNA_UNUAR_seq),
  set = "UNUAR"
)
background_df <- tibble(
  sequence = as.character(dRNA_background_seq),
  set = "Background"
)
background_unuar <- background_df %>%
  filter(str_detect(sequence, "T.TA[AG]"))

# Combine them into one master dataframe
combined_seqs <- bind_rows(unuar_df, background_unuar)


# --- PART 2: Absolute Counts & Relative Enrichment ---
composition_summary <- combined_seqs %>%
  group_by(set) %>%
  summarize(
    Total_Rows = n(),
    TGTAG = sum(sequence == "TGTAG", na.rm = TRUE),
    TCTAG = sum(sequence == "TCTAG", na.rm = TRUE),
    TNTAR = sum(str_detect(sequence, "^T[ACGT]TA[AG]$"), na.rm = TRUE),
    Combined_TG_TC = sum(sequence %in% c("TGTAG", "TCTAG"), na.rm = TRUE)
  ) %>%
  pivot_longer(
    cols = c(TGTAG, TCTAG, TNTAR, Combined_TG_TC),
    names_to = "Motif",
    values_to = "Absolute_Count"
  ) %>%
  # Calculate relative enrichment (percentage)
  mutate(Relative_Percent = (Absolute_Count / Total_Rows) * 100) %>%
  select(set, Motif, Absolute_Count, Total_Rows, Relative_Percent) %>%
  arrange(set, Motif)

cat("--- Absolute Counts & Relative Enrichment ---\n")
print(composition_summary)

# --- PART 3: Fisher's Exact Test (UNUAR vs Background) ---
fisher_results_unuar <- composition_summary %>%
  filter(Motif %in% c("TGTAG", "TCTAG", "Combined_TG_TC")) %>%
  filter(set %in% c("UNUAR", "Background")) %>%
  mutate(Non_Motif_Count = Total_Rows - Absolute_Count) %>%
  select(set, Motif, Absolute_Count, Non_Motif_Count) %>%
  pivot_wider(
    names_from = set,
    values_from = c(Absolute_Count, Non_Motif_Count)
  ) %>%
  mutate(
    # Run the Fisher's exact test for each motif
    stats = pmap(
      list(Absolute_Count_UNUAR, Non_Motif_Count_UNUAR, 
           Absolute_Count_Background, Non_Motif_Count_Background),
      ~ tidy(fisher.test(matrix(c(..1, ..2, ..3, ..4), nrow = 2)))
    )
  ) %>%
  unnest(stats) %>%
  select(
    Motif,
    Count_UNUAR = Absolute_Count_UNUAR,
    Count_Background = Absolute_Count_Background,
    p.value,
    Odds_Ratio = estimate
  )

cat("\n--- Fisher's Exact Test Results (UNUAR vs Background) ---\n")
print(fisher_results_unuar)

# ----- Stacked Bar Chart -----
theme_base_custom <- function(base_size = 16) {
  theme_minimal(base_size = base_size, base_family = "sans") %+replace%
    theme(
      plot.title = element_text(hjust = 0.5, size = rel(1.2), margin = margin(b = 5)),
      plot.subtitle = element_text(hjust = 0.5, size = rel(0.8), margin = margin(b = 5)),
      axis.title = element_text(size = rel(1.0)),
      axis.text = element_text(size = rel(0.8)),
      plot.margin = margin(5, 5, 5, 5)
    )
}
theme_bar <- function(base_size = 16) {
  theme_base_custom(base_size) +
    theme(legend.position = "none", panel.grid.major.x = element_blank(), panel.grid.minor = element_blank())
}

unuar_df <- tibble(
  sequence = as.character(dRNA_UNUAR_seq),
  set = "PUS7-dependent"
)
background_df <- tibble(
  sequence = as.character(dRNA_background_seq),
  set = "Background"
)
background_unuar <- background_df %>%
  filter(str_detect(sequence, "T.TA[AG]"))

combined_seqs <- bind_rows(unuar_df, background_unuar)

plot_data <- combined_seqs %>%
  mutate(motif = case_when(
    str_detect(sequence, "TGTAG") ~ "UGUAG",
    str_detect(sequence, "TCTAG") ~ "UCUAG",
    TRUE ~ "Other UNUAR"
  )) %>%
  group_by(set, motif) %>%
  summarise(Count = n(), .groups = "drop") %>%
  group_by(set) %>%
  mutate(Percentage = (Count / sum(Count)) * 100)

plot_data$motif <- factor(plot_data$motif, levels = c("UGUAG", "UCUAG", "Other UNUAR"))
plot_data$set <- factor(plot_data$set, levels = c("Background", "PUS7-dependent"))

ggplot(plot_data, aes(x = set, y = Percentage, fill = motif)) +
  geom_bar(stat = "identity", width = 0.5, color = "black") +
  scale_fill_manual(values = c("Other UNUAR" = "gray70", 
                               "UCUAG" = "#d598b7",
                               "UGUAG" = "#a0cdb9")) + 
  theme_bar(base_size = 10) +
  labs(
    title = "direct RNA Motifs",
    x = NULL, 
    y = "Percentage of Sites (%)",
    fill = "Motif"
  ) + 
  theme(legend.position = "right", axis.text.x = element_text(angle = 45, hjust = 1))

plot_path <- paste0(active_dir, "/1F_bar_dRNA_UNUAR_motif_composition")
ggsave(paste0(plot_path, ".png"), width = 3, height = 3)
ggsave(paste0(plot_path, ".pdf"), width = 3, height = 3)
