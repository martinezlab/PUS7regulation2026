#!/usr/bin/env Rscript

library(dplyr)
library(stringr)
library(readr)
library(ggplot2)
library(tidyr)
library(purrr)
library(broom)

# ---- Directory ----
active_dir <- "~/PUS7regulation2026/Figure4/sequence"

pool1 <- read.csv(paste0(active_dir, "/pool1_background.csv"))

# ---- in cell UNUAR only ----

# Identify modified sites
incell_full <- read.csv("incellBID_WTsummary.csv", header = TRUE)
incell_modified <- incell_full %>% 
    filter(str_detect(psipos_motif, "^T[ACGT]TA[AG]$"), PUS7_dependent == "Yes") %>%
    select(site, psipos_motif) %>%
    mutate(Comparison_Group = "Modified")

# Identify background distribution
unmodified_data <- anti_join(pool1, incell_modified, by = join_by("clean_name" == "site")) %>%
    filter(str_detect(psipos_5mer, "^T[ACGT]TA[AG]$")) %>%
    select(site = clean_name, psipos_motif = psipos_5mer) %>%
    mutate(Comparison_Group = "Unmodified")

# Absolute Counts
summary_counts <- bind_rows(unmodified_data, incell_modified) %>%
  group_by(Comparison_Group) %>%
  summarize(
    Total_TNTAR_Rows = n(),
    TGTAG = sum(psipos_motif == "TGTAG", na.rm = TRUE),
    TCTAG = sum(psipos_motif == "TCTAG", na.rm = TRUE)
  )

# Create final report
final_report <- summary_counts %>%
  pivot_longer(
    cols = c(TGTAG, TCTAG), 
    names_to = "Target_Motif", 
    values_to = "Motif_Count"
  ) %>%
  mutate(Non_Motif_Count = Total_TNTAR_Rows - Motif_Count) %>%
  pivot_wider(
    names_from = Comparison_Group, 
    values_from = c(Motif_Count, Non_Motif_Count, Total_TNTAR_Rows),
    values_fill = 0
  ) %>%
  mutate(
    Percent_Unmodified = (Motif_Count_Unmodified / Total_TNTAR_Rows_Unmodified) * 100,
    Percent_Modified = (Motif_Count_Modified / Total_TNTAR_Rows_Modified) * 100,
    stats = pmap(
      list(Motif_Count_Modified, Non_Motif_Count_Modified, Motif_Count_Unmodified, Non_Motif_Count_Unmodified),
      ~ tidy(fisher.test(matrix(c(..1, ..2, ..3, ..4), nrow = 2)))
    )
  ) %>%
  unnest(stats) %>%
  select(Target_Motif, Percent_Unmodified, Percent_Modified, p.value, Odds_Ratio = estimate)

cat("in cellulo UNUAR enrichment")
print(final_report)

# ---- in vitro UNUAR only ----

# Identify modified sites
invitro_full <- read.csv("invitroBID_MPRA_summary.csv", header = TRUE)
invitro_modified <- invitro_full %>% 
    filter(str_detect(psipos_motif, "^T[ACGT]TA[AG]$"), PUS7_dependent == "Yes") %>%
    select(site, psipos_motif) %>%
    mutate(Comparison_Group = "Modified")

# Identify background distribution
unmodified_data <- anti_join(pool1, invitro_modified, by = join_by("clean_name" == "site")) %>%
    filter(str_detect(psipos_5mer, "^T[ACGT]TA[AG]$")) %>%
    select(site = clean_name, psipos_motif = psipos_5mer) %>%
    mutate(Comparison_Group = "Unmodified")

# Absolute Counts
summary_counts <- bind_rows(unmodified_data, invitro_modified) %>%
  group_by(Comparison_Group) %>%
  summarize(
    Total_TNTAR_Rows = n(),
    TGTAG = sum(psipos_motif == "TGTAG", na.rm = TRUE),
    TCTAG = sum(psipos_motif == "TCTAG", na.rm = TRUE)
  )

# Create final report
final_report <- summary_counts %>%
  pivot_longer(
    cols = c(TGTAG, TCTAG), 
    names_to = "Target_Motif", 
    values_to = "Motif_Count"
  ) %>%
  mutate(Non_Motif_Count = Total_TNTAR_Rows - Motif_Count) %>%
  pivot_wider(
    names_from = Comparison_Group, 
    values_from = c(Motif_Count, Non_Motif_Count, Total_TNTAR_Rows),
    values_fill = 0
  ) %>%
  mutate(
    Percent_Unmodified = (Motif_Count_Unmodified / Total_TNTAR_Rows_Unmodified) * 100,
    Percent_Modified = (Motif_Count_Modified / Total_TNTAR_Rows_Modified) * 100,
    stats = pmap(
      list(Motif_Count_Modified, Non_Motif_Count_Modified, Motif_Count_Unmodified, Non_Motif_Count_Unmodified),
      ~ tidy(fisher.test(matrix(c(..1, ..2, ..3, ..4), nrow = 2)))
    )
  ) %>%
  unnest(stats) %>%
  select(Target_Motif, Percent_Unmodified, Percent_Modified, p.value, Odds_Ratio = estimate)

cat("in vitro UNUAR enrichment")
print(final_report)


# ---- Stacked Bar Chart ----

pool1_unuar <- pool1 %>%
  filter(str_detect(psipos_5mer, "T.TA[AG]"))

pool1_bg <- pool1_unuar %>%
  select(sequence = psipos_5mer) %>%
  mutate(set = "Background")

pool1_invitro <- pool1_unuar %>%
  filter(clean_name %in% invitro_modified$site) %>%
  select(sequence = psipos_5mer) %>%
  mutate(set = "in vitro PUS7")

pool1_incell <- pool1_unuar %>%
  filter(clean_name %in% incell_modified$site) %>%
  select(sequence = psipos_5mer) %>%
  mutate(set = "in cell PUS7")

pool1_plot_data <- bind_rows(pool1_bg, pool1_invitro, pool1_incell) %>%
  mutate(motif = case_when(
    str_detect(sequence, "TGTAG") ~ "UGUAG",
    str_detect(sequence, "TCTAG") ~ "UCUAG",
    TRUE ~ "Other UNUAR"
  )) %>%
  group_by(set, motif) %>%
  summarise(Count = n(), .groups = "drop") %>%
  group_by(set) %>%
  mutate(Percentage = (Count / sum(Count)) * 100)

pool1_plot_data$motif <- factor(pool1_plot_data$motif, levels = c("UGUAG", "UCUAG", "Other UNUAR"))
pool1_plot_data$set <- factor(pool1_plot_data$set, levels = c("Background", "in vitro PUS7", "in cell PUS7"))

ggplot(pool1_plot_data, aes(x = set, y = Percentage, fill = motif)) +
  geom_bar(stat = "identity", width = 0.5, color = "black") +
  scale_fill_manual(values = c("Other UNUAR" = "gray70", 
                               "UCUAG" = "#d598b7",
                               "UGUAG" = "#a0cdb9")) + 
  theme_bar(base_size = 10) +
  labs(
    title = "Pool 1 Motifs",
    x = NULL, 
    y = "Percentage of Sites (%)",
    fill = "Motif"
  ) + 
  theme(legend.position = "right", axis.text.x = element_text(angle = 45, hjust = 1))

plot_path <- paste0(active_dir, "/4B_bar_MPRA_UNUAR_motif_composition")
ggsave(paste0(plot_path, ".png"), width = 4, height = 3)
ggsave(paste0(plot_path, ".pdf"), width = 4, height = 3)