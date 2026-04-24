#!/usr/bin/env Rscript

library(dplyr)
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

# ---- Data ----

active_dir <- "~/PUS7regulation2026/Figure5/mutagenesis_analysis"
pool2_full <- read.csv(file.path(active_dir, "mutagenesis_data.csv"))

# ---- Accessibility - Psi ----

access0_make <- pool2_full %>%
  filter(
    sequence_type %in% c("makepair", "original_pool1")
  ) %>%
  group_by(parent_name) %>%
  filter(n_distinct(sequence_type) == 2) %>%
  ungroup() %>%
  mutate(sequence_type = factor(sequence_type, levels = c("original_pool1", "makepair")))

access0_break <- pool2_full %>%
  filter(
    sequence_type %in% c("breakpair", "original_pool1")
  ) %>%
  group_by(parent_name) %>%
  filter(n_distinct(sequence_type) == 2) %>%
  ungroup() %>%
  mutate(sequence_type = factor(sequence_type, levels = c("original_pool1", "breakpair")))

# --- Boxplot ---

p_access0_make_boxplot <- ggplot(access0_make, aes(x = sequence_type, y = delta_delrate)) +
  geom_boxplot(fill = invitro_color) +
  theme_boxplot() +
  labs(
    title = "Make Pair at Target Uridine",
    subtitle = paste0("Based on ", nrow(access0_make)/2, " sites"),
    x = "Mutation Type",
    y = "Delta Deletion Fraction"
  ) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1)
  )

plot_path <- file.path(active_dir, "access_psi_make_boxplot")
ggsave(paste0(plot_path, ".pdf"), p_access0_make_boxplot, width = 4, height = 5, units = "in")
ggsave(paste0(plot_path, ".png"), p_access0_make_boxplot, width = 4, height = 5, units = "in", dpi = 600)

p_access0_break_boxplot <- ggplot(access0_break, aes(x = sequence_type, y = delta_delrate)) +
  geom_boxplot(fill = invitro_color) +
  theme_boxplot() +
  labs(
    title = "Break Pair at Target Uridine",
    subtitle = paste0("Based on ", nrow(access0_break)/2, " sites"),
    x = "Mutation Type",
    y = "Delta Deletion Fraction"
  ) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1)
  )

plot_path <- file.path(active_dir, "access_psi_break_boxplot")
ggsave(paste0(plot_path, ".pdf"), p_access0_break_boxplot, width = 4, height = 5, units = "in")
ggsave(paste0(plot_path, ".png"), p_access0_break_boxplot, width = 4, height = 5, units = "in", dpi = 600)


# ---- Accessibility - Psi and +1 ----

access01_make <- pool2_full %>%
  filter(
    sequence_type %in% c("makepair01", "original_pool1")
  ) %>%
  group_by(parent_name) %>%
  filter(n_distinct(sequence_type) == 2) %>%
  ungroup() %>%
  mutate(sequence_type = factor(sequence_type, levels = c("original_pool1", "makepair01")))

access01_break <- pool2_full %>%
  filter(
    sequence_type %in% c("breakpair01", "original_pool1")
  ) %>%
  group_by(parent_name) %>%
  filter(n_distinct(sequence_type) == 2) %>%
  ungroup() %>%
  mutate(sequence_type = factor(sequence_type, levels = c("original_pool1", "breakpair01")))

# --- Boxplot ---

p_access01_make_boxplot <- ggplot(access01_make, aes(x = sequence_type, y = delta_delrate)) +
  geom_boxplot(fill = invitro_color) +
  theme_boxplot() +
  labs(
    title = "Make Pair at Target Uridine and +1",
    subtitle = paste0("Based on ", nrow(access01_make)/2, " sites"),
    x = "Mutation Type",
    y = "Delta Deletion Fraction"
  ) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1)
  )

plot_path <- file.path(active_dir, "access_psi01_make_boxplot")
ggsave(paste0(plot_path, ".pdf"), p_access01_make_boxplot, width = 4, height = 5, units = "in")
ggsave(paste0(plot_path, ".png"), p_access01_make_boxplot, width = 4, height = 5, units = "in", dpi = 600)

p_access01_break_boxplot <- ggplot(access01_break, aes(x = sequence_type, y = delta_delrate)) +
  geom_boxplot(fill = invitro_color) +
  theme_boxplot() +
  labs(
    title = "Break Pair at Target Uridine and +1",
    subtitle = paste0("Based on ", nrow(access01_break)/2, " sites"),
    x = "Mutation Type",
    y = "Delta Deletion Fraction"
  ) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1)
  )

plot_path <- file.path(active_dir, "access_psi01_break_boxplot")
ggsave(paste0(plot_path, ".pdf"), p_access01_break_boxplot, width = 4, height = 5, units = "in")
ggsave(paste0(plot_path, ".png"), p_access01_break_boxplot, width = 4, height = 5, units = "in", dpi = 600)


# ---- Accessibility - 5mer ----

access5mer_make <- pool2_full %>%
  filter(
    sequence_type %in% c("makepair5mer", "original_pool1")
  ) %>%
  group_by(parent_name) %>%
  filter(n_distinct(sequence_type) == 2) %>%
  ungroup() %>%
  mutate(sequence_type = factor(sequence_type, levels = c("original_pool1", "makepair5mer")))

access5mer_break <- pool2_full %>%
  filter(
    sequence_type %in% c("breakpair5mer", "original_pool1")
  ) %>%
  group_by(parent_name) %>%
  filter(n_distinct(sequence_type) == 2) %>%
  ungroup() %>%
  mutate(sequence_type = factor(sequence_type, levels = c("original_pool1", "breakpair5mer")))

# --- Boxplot ---

p_access5mer_make_boxplot <- ggplot(access5mer_make, aes(x = sequence_type, y = delta_delrate)) +
  geom_boxplot(fill = invitro_color) +
  theme_boxplot() +
  labs(
    title = "Make Pair at 5mer",
    subtitle = paste0("Based on ", nrow(access5mer_make)/2, " sites"),
    x = "Mutation Type",
    y = "Delta Deletion Fraction"
  ) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1)
  )

plot_path <- file.path(active_dir, "access_5mer_make_boxplot")
ggsave(paste0(plot_path, ".pdf"), p_access5mer_make_boxplot, width = 4, height = 5, units = "in")
ggsave(paste0(plot_path, ".png"), p_access5mer_make_boxplot, width = 4, height = 5, units = "in", dpi = 600)

p_access5mer_break_boxplot <- ggplot(access5mer_break, aes(x = sequence_type, y = delta_delrate)) +
  geom_boxplot(fill = invitro_color) +
  theme_boxplot() +
  labs(
    title = "Break Pair at 5mer",
    subtitle = paste0("Based on ", nrow(access5mer_break)/2, " sites"),
    x = "Mutation Type",
    y = "Delta Deletion Fraction"
  ) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1)
  )

plot_path <- file.path(active_dir, "access_5mer_break_boxplot")
ggsave(paste0(plot_path, ".pdf"), p_access5mer_break_boxplot, width = 4, height = 5, units = "in")
ggsave(paste0(plot_path, ".png"), p_access5mer_break_boxplot, width = 4, height = 5, units = "in", dpi = 600)


# ---- Accessibility - Unpaired2 ----

access_unpaired2_make <- pool2_full %>%
  filter(
    sequence_type %in% c("makepair_unpaired2", "original_pool1")
  ) %>%
  group_by(parent_name) %>%
  filter(n_distinct(sequence_type) == 2) %>%
  ungroup() %>%
  mutate(sequence_type = factor(sequence_type, levels = c("original_pool1", "makepair_unpaired2")))

access_unpaired2_break <- pool2_full %>%
  filter(
    sequence_type %in% c("breakpair_unpaired2", "original_pool1")
  ) %>%
  group_by(parent_name) %>%
  filter(n_distinct(sequence_type) == 2) %>%
  ungroup() %>%
  mutate(sequence_type = factor(sequence_type, levels = c("original_pool1", "breakpair_unpaired2")))

# --- Boxplot ---

p_access_unpaired2_make_boxplot <- ggplot(access_unpaired2_make, aes(x = sequence_type, y = delta_delrate)) +
  geom_boxplot(fill = invitro_color) +
  theme_boxplot() +
  labs(
    title = "Make Pair at Unpaired2 (+6,7,8)",
    subtitle = paste0("Based on ", nrow(access_unpaired2_make)/2, " sites"),
    x = "Mutation Type",
    y = "Delta Deletion Fraction"
  ) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1)
  )

plot_path <- file.path(active_dir, "access_unpaired2_make_boxplot")
ggsave(paste0(plot_path, ".pdf"), p_access_unpaired2_make_boxplot, width = 4, height = 5, units = "in")
ggsave(paste0(plot_path, ".png"), p_access_unpaired2_make_boxplot, width = 4, height = 5, units = "in", dpi = 600)

p_access_unpaired2_break_boxplot <- ggplot(access_unpaired2_break, aes(x = sequence_type, y = delta_delrate)) +
  geom_boxplot(fill = invitro_color) +
  theme_boxplot() +
  labs(
    title = "Break Pair at Unpaired2 (+6,7,8)",
    subtitle = paste0("Based on ", nrow(access_unpaired2_break)/2, " sites"),
    x = "Mutation Type",
    y = "Delta Deletion Fraction"
  ) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1)
  )

plot_path <- file.path(active_dir, "access_unpaired2_break_boxplot")
ggsave(paste0(plot_path, ".pdf"), p_access_unpaired2_break_boxplot, width = 4, height = 5, units = "in")
ggsave(paste0(plot_path, ".png"), p_access_unpaired2_break_boxplot, width = 4, height = 5, units = "in", dpi = 600)
