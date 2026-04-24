#!/usr/bin/env Rscript

library(dplyr)
library(ggplot2)

# ---- Plot Themes ----
disruption <- "#fce87e"
compensation <- "#cfa2fc"
invitro <- "#a2a2fc"
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
theme_bar <- function(base_size = 16) {
  theme_base_custom(base_size) +
    theme(legend.position = "none", panel.grid.major.x = element_blank(), panel.grid.minor = element_blank())
}

create_structure_barplot <- function(sequence_name, input_data) {

  plot_data <- input_data %>%
    filter(parent_name == sequence_name) %>%
    filter(sequence_type %in% c("original_pool1", "total_disruption", "total_compensation")) %>%
    mutate(
      labels = case_when(
        sequence_type == "original_pool1" ~ "Parental",
        sequence_type == "total_disruption" ~ "Disruption",
        sequence_type == "total_compensation" ~ "Compensation",
        TRUE ~ sequence_type
      ),
      labels = factor(labels, levels = c("Parental", "Disruption", "Compensation"))
    )
  plot_title <- sub("_", ":", sub("_", " ", sequence_name))
  gene_name <- strsplit(sequence_name, "_")[[1]][1]
  plot_base_path <- file.path(active_dir, paste0("barplot_structure_", gene_name))

  p <- ggplot(plot_data, aes(x = labels, y = delta_delrate, fill = labels)) +
    stat_summary(fun = mean, geom = "bar", aes(fill = labels)) +
    stat_summary(fun.data = mean_se, geom = "errorbar", width = 0.2) +
    geom_jitter(width = 0.2, alpha = 0.5) +
    scale_y_continuous(expand = expansion(mult = c(0, 0.05))) +
    scale_fill_manual(values = c("Parental" = invitro, "Disruption" = disruption, "Compensation" = compensation)) +
    labs(
      title = plot_title,
      x = "Mutation Type",
      y = "Deletion Fraction"
    ) +
    theme_bar(base_size = 10) +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      legend.position = "none"
    )
  ggsave(paste0(plot_base_path, ".pdf"), p, width = 2, height = 2.5, units = "in")
  ggsave(paste0(plot_base_path, ".png"), p, width = 2, height = 2.5, units = "in", dpi = 600)
}

# ---- Data ----

active_dir <- "~/PUS7regulation2026/Figure5/mutagenesis_analysis"
pool2_full <- read.csv(file.path(active_dir, "mutagenesis_data.csv"))

# ---- Structure Boxplots ----
create_structure_barplot(
  sequence_name = "RHBDD2_chr7_75888787",
  input_data = pool2_full
)
create_structure_barplot(
  sequence_name = "PPP1R12B_chr1_202587327",
  input_data = pool2_full
)
create_structure_barplot(
  sequence_name = "IGF2BP1_chr17_49043468",
  input_data = pool2_full
)
