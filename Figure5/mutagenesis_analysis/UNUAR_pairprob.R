#!/usr/bin/env Rscript

library(dplyr)
library(tidyr)
library(stringr)
library(ggplot2)

# ---- Plot Themes ----
input_color <- "#eee8aa"
invitro_color <- "#a2a2fc"
unmodified <- "#d3d3d3"

theme_base_custom <- function(base_size = 16) {
  theme_minimal(base_size = base_size, base_family = "sans") %+replace%
    theme(
      plot.title = element_text(hjust = 0.5, size = rel(1.2), margin = margin(b = 5)),
      plot.subtitle = element_text(hjust = 0.5, size = rel(0.8), margin = margin(b = 5)),
      axis.title = element_text(size = rel(1.0)),
      axis.text = element_text(size = rel(1.0)),
      plot.margin = margin(5, 5, 5, 5)
    )
}
theme_boxplot <- function(base_size = 16) {
  theme_base_custom(base_size) +
    theme(legend.position = "none", panel.grid.major.x = element_blank(), panel.grid.minor = element_blank())
}
theme_line <- function(base_size = 16) {
  theme_base_custom(base_size) +
    theme(legend.position = "right", panel.grid.minor = element_blank())
}

# ---- Data ----

active_dir <- "~/PUS7regulation2026/Figure5/mutagenesis_analysis"
pool2_full <- read.csv(file.path(active_dir, "mutagenesis_data.csv"))

unuar_mutations <- pool2_full %>%
  filter(
    heatmap_27 == TRUE,
    sequence_type %in% c("UNUAR", "original_pool1")
    )

unuar_pairprob <- read.csv(file.path(active_dir, "UNUAR_pairprob.csv"))

# ---- Pairing Probabilities - Line Graphs ----
selected_unuar_mutations <- unuar_mutations %>%
  mutate(psipos_5mer = str_replace_all(psipos_5mer, "T", "U")) %>%
  filter(parent_name %in% c("ZNF74_chr22_20405560", "STIM1_chr11_4092507", "IGF2BP1_chr17_49043468"))

POSITIONS_OF_INTEREST <- c(63, 64, 65, 66, 67)
prob_col_names <- paste0("X", POSITIONS_OF_INTEREST)
manual_order <- c("UUUAA", "UAUAA", "UCUAA", "UGUAA", "UAUAG", "UUUAG", "UGUAG", "UCUAG")

plot_data <- selected_unuar_mutations %>%
  left_join(
    select(unuar_pairprob, rna_name, all_of(prob_col_names)),
    by = c("chr" = "rna_name")
  ) %>%
  mutate(
    prob_5mer_avg = rowMeans(select(., all_of(prob_col_names)), na.rm = TRUE),
    psipos_5mer_ordered = factor(psipos_5mer, levels = manual_order)) %>%
  pivot_longer(
    cols = c(delta_delrate, X65, X66, prob_5mer_avg),
    names_to = "metric",
    values_to = "value"
  ) %>%
  mutate(metric = case_when(
    metric == "delta_delrate"   ~ "Delta Deletion Fraction",
    metric == "X65"             ~ "Pairing Prob (Target Uridine)",
    metric == "X66"             ~ "Pairing Prob (Pos +1)",
    metric == "prob_5mer_avg"   ~ "Pairing Prob (Target Motif)",
    TRUE                        ~ metric
  ))

all_metric_names <- unique(plot_data$metric)
black_line_metric <- "Delta Deletion Fraction"
other_metrics <- all_metric_names[all_metric_names != black_line_metric]
manual_colors <- c(
  setNames("black", black_line_metric),
  setNames(scales::hue_pal()(length(other_metrics)), other_metrics)
)

p_unuar_line_combined <- ggplot(plot_data, 
  aes(x = psipos_5mer_ordered, y = value, group = metric, color = metric)) +
  geom_line(linewidth = 1) +
  geom_point(size = 2) +
  scale_color_manual(values = manual_colors) +
  facet_wrap(~parent_name, ncol = 5) +
  labs(
    x = "UNUAR Motif",
    y = "Value",
    color = "Metric"
  ) +
  theme_line(base_size = 12) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "bottom"
  )

plot_path <- file.path(active_dir, "line_unuar_pairprob")
ggsave(paste0(plot_path, ".pdf"), p_unuar_line_combined, width = 6, height = 3, units = "in")
ggsave(paste0(plot_path, ".png"), p_unuar_line_combined, width = 6, height = 4, units = "in", dpi = 600)

# ---- Pairing Probabilities - Boxplots ----

delrate_comparison <- unuar_mutations %>%
  filter(psipos_5mer %in% c("UGUAG", "UCUAG")) %>%
  select(parent_name, psipos_5mer, delta_delrate) %>%
  pivot_wider(
    names_from = psipos_5mer,
    values_from = delta_delrate
  ) %>%
  na.omit() %>%
  mutate(
    delrate_diff = UGUAG - UCUAG,
    group = case_when(
      delrate_diff > 0.05  ~ "UGUAG > UCUAG",
      delrate_diff < -0.05 ~ "UGUAG < UCUAG",
      TRUE                 ~ "Equal"
    )
  ) %>%
  select(parent_name, group)

prob_plot_data_split <- unuar_mutations %>%
  filter(psipos_5mer %in% c("UGUAG", "UCUAG")) %>%
  inner_join(delrate_comparison, by = "parent_name") %>%
  inner_join(
    select(unuar_pairprob, rna_name, all_of(prob_col_names)),
    by = c("chr" = "rna_name")
  ) %>%
  rowwise() %>%
  mutate(XAverage = mean(c_across(all_of(prob_col_names)))) %>%
  ungroup() %>%
  pivot_longer(
    cols = c(all_of(prob_col_names), XAverage),
    names_to = "metric",
    values_to = "pair_prob"
  ) %>%
  filter(metric %in% c("X65", "X66", "XAverage")) %>%
  mutate(
    metric = case_when(
      metric == "X65"      ~ "Target Uridine",
      metric == "X66"      ~ "Pos +1",
      metric == "XAverage" ~ "Target Motif",
      TRUE                 ~ metric 
    ),
    group = factor(group, levels = c("UGUAG > UCUAG", "Equal", "UGUAG < UCUAG")),
    metric = factor(metric, levels = c(
      "Target Uridine", "Pos +1", "Target Motif"
    ))
  )

p_unuar_prob_boxplot <- ggplot(prob_plot_data_split, aes(x = group, y = pair_prob, fill = psipos_5mer)) +
  geom_boxplot(position = position_dodge(width = 0.85), width = 0.75) +
  facet_wrap(~metric, ncol = 3) +
  scale_fill_manual(
    name = "Mutation Motif",
    values = c("UGUAG" = "#a0cdb9", "UCUAG" = "#d598b7")
  ) +
  labs(
    x = "Deletion Fraction Comparison",
    y = "Pairing Probability"
  ) +
  theme_boxplot(base_size = 12) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "bottom"
  )
plot_path <- file.path(active_dir, "boxplot_UNUAR_pairprob")
ggsave(paste0(plot_path, ".pdf"), p_unuar_prob_boxplot, width = 6, height = 3.5, units = "in")
ggsave(paste0(plot_path, ".png"), p_unuar_prob_boxplot, width = 6, height = 3.5, units = "in", dpi = 600)
