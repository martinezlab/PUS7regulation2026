#!/usr/bin/env Rscript

library(dplyr)
library(ggplot2)
library(tidyr)
library(forcats)

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

# ---- Data ----

active_dir <- "~/PUS7regulation2026/Figure5/mutagenesis_analysis"
pool2_full <- read.csv(file.path(active_dir, "mutagenesis_data.csv"))

unuar_mutations <- pool2_full %>%
  filter(
    heatmap_27 == TRUE,
    sequence_type %in% c("UNUAR", "original_pool1")
    )

# --- Summary Stats ---

# Delta Deletion Rate
delrate_summary <- unuar_mutations %>%
  group_by(psipos_5mer) %>%
  summarise(
    mean_delta_delrate = mean(delta_delrate, na.rm = TRUE),
    min_delta_delrate = min(delta_delrate, na.rm = TRUE),
    max_delta_delrate = max(delta_delrate, na.rm = TRUE),
    .groups = 'drop'
  )

# Modified or Not
category_counts <- unuar_mutations %>%
  count(psipos_5mer, category) %>%
  pivot_wider(
    names_from = category, 
    values_from = n, 
    values_fill = 0
  )

# Join the two summary tables together.
full_summary_table <- delrate_summary %>%
  left_join(category_counts, by = "psipos_5mer")


# --- Boxplot ---

p_unuar_boxplot <- ggplot(unuar_mutations, aes(x = fct_reorder(psipos_5mer, delta_delrate), y = delta_delrate)) +
  geom_boxplot(fill = invitro_color) +
  theme_boxplot() +
  labs(
    title = "UNUAR Mutations",
    x = "UNUAR Motif",
    y = "Delta Deletion Fraction"
  ) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1)
  )

plot_path <- file.path(active_dir, "UNUAR_boxplot")
ggsave(paste0(plot_path, ".pdf"), p_unuar_boxplot, width = 6, height = 4, units = "in")
ggsave(paste0(plot_path, ".png"), p_unuar_boxplot, width = 6, height = 4, units = "in", dpi = 600)
