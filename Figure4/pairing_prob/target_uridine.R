#!/usr/bin/env Rscript

library(dplyr)
library(readr)
library(tidyr)
library(ggplot2)
library(stringr)

# ---- Plot Themes ----
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
theme_boxplot <- function(base_size = 16) {
  theme_base_custom(base_size) +
    theme(legend.position = "none", panel.grid.major.x = element_blank(), panel.grid.minor = element_blank())
}
incell <- "#e25098"
invitro <- "#a2a2fc"

plot_dir <- "~/PUS7regulation2026/Figure4/plots" # Adjust to your actual target folder
# ---- Datasets ----
ivSHAPE_pairprob_pool1 <- read.csv("IV_NoPus_pairprob_aligned.csv")

icPUS7_BothQuartiles <- read_tsv("~/PUS7regulation2026/Figure4/sequence/incell_quartiles/PUS7_dep_union_by_Both_WT_quartiles.tsv", col_names = TRUE)

ivPUS7_Q1 <- read_tsv("~/PUS7regulation2026/Figure4/sequence/invitro_quartiles/invitro_modified_sites_Q1_summary.tsv", col_names = TRUE)
ivPUS7_Q2 <- read_tsv("~/PUS7regulation2026/Figure4/sequence/invitro_quartiles/invitro_modified_sites_Q2_summary.tsv", col_names = TRUE)
ivPUS7_Q3 <- read_tsv("~/PUS7regulation2026/Figure4/sequence/invitro_quartiles/invitro_modified_sites_Q3_summary.tsv", col_names = TRUE)
ivPUS7_Q4 <- read_tsv("~/PUS7regulation2026/Figure4/sequence/invitro_quartiles/invitro_modified_sites_Q4_summary.tsv", col_names = TRUE)

# ---- Execution ----
ivPUS7_SHAPE_pairprob <- ivSHAPE_pairprob_pool1 %>%
  mutate(quartile = case_when(
    rna_name %in% ivPUS7_Q1$chr ~ "1",
    rna_name %in% ivPUS7_Q2$chr ~ "2",
    rna_name %in% ivPUS7_Q3$chr ~ "3",
    rna_name %in% ivPUS7_Q4$chr ~ "4"
  )) %>%
  filter(!(is.na(quartile)))

p <- ggplot(ivPUS7_SHAPE_pairprob, aes(x = quartile, y = pos_0)) +
  geom_boxplot(fill = invitro) +
  labs(
    x = "Quartile",
    y = "Target Uridine Pairing Probability",
    title = "In Vitro"
  ) +
  theme_boxplot(base_size = 12)
path <- file.path(plot_dir, "4_boxplot_psipairprob_invitro")
ggsave(paste0(path, ".pdf"), p, width = 1.8, height = 2.25, units = "in")
p_big <- p + theme_boxplot(base_size = 24)
ggsave(paste0(path, ".png"), p_big, width = 6, height = 6, units = "in")

# with WT both data
icPUS7_SHAPE_pairprob <- ivSHAPE_pairprob_pool1 %>%
  left_join(icPUS7_BothQuartiles, by = c("rna_name" = "chr")) %>%
  filter(!(is.na(quartile))) %>% 
  mutate(quartile = as.character(quartile))

p <- ggplot(icPUS7_SHAPE_pairprob, aes(x = quartile, y = pos_0)) +
  geom_boxplot(fill = incell) +
  labs(
    x = "Quartile",
    y = "Target Uridine Pairing Probability",
    title = "In Cellulo"
  ) +
  theme_boxplot(base_size = 12)
path <- file.path(plot_dir, "4_boxplot_psipairprob_incell_both")
ggsave(paste0(path, ".pdf"), p, width = 1.8, height = 2.25, units = "in")
p_big <- p + theme_boxplot(base_size = 24)
ggsave(paste0(path, ".png"), p_big, width = 6, height = 6, units = "in")