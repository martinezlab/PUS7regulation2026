#!/usr/bin/env Rscript

library(dplyr)
library(tidyr)
library(ggplot2)
library(scales)

# ---- Directory ----
active_dir <- "~/PUS7regulation2026/Figure1/UNUARstandard"

# ---- Datasets ----

formulas <- read.csv(paste0(active_dir,"/derivation_formulas_subtract.csv"))

dRNA_sites <- read.csv(paste0(active_dir,"/PUS7_10reads_15percent_motif.csv"))

dRNA_processed <- dRNA_sites %>%
  rowwise() %>%
  mutate(
    avg_mm_WT = mean(c(mm.perc.rep2, mm.perc.rep3, mm.perc.rep4), na.rm = TRUE),
    avg_mm_KD = mean(c(ctrl.err.rep2, ctrl.err.rep3, ctrl.err.rep4), na.rm = TRUE),
    motif = gsub("T", "U", motif)
    ) %>%
  ungroup() %>% 
  select(chr, pos, gene_name, motif, gene_strand, avg_mm_WT, avg_mm_KD) %>%
  mutate(chr_pos_motif = paste0(chr, "_", pos, "_", motif))

# ---- Apply Formulas ----

dRNA_adjusted <- dRNA_processed %>%
  left_join(formulas, by = "motif") %>%
  mutate(
    adj_mm_WT = (avg_mm_WT - intercept) / slope,
    adj_mm_KD = (avg_mm_KD - intercept) / slope
  ) %>%
  select(chr, pos, chr_pos_motif, gene_name, gene_strand, motif, avg_mm_WT, avg_mm_KD, adj_mm_WT, adj_mm_KD)

write.csv(x = dRNA_adjusted, 
          file = paste0(active_dir, "/UNUAR/dRNA_UtoC_adjusted_subtract.csv"), row.names = FALSE)

# ---- Plot ----

dRNA_WT <- "#daa520"
dRNA_KD <- "#edd290"
dRNA_adj_WT <- "#bc8f8f"
dRNA_adj_KD <- "#dec7c7"

theme_boxplot <- function(base_size = 18, base_family = "sans") {
  theme_minimal(base_size = base_size) %+replace%
    theme(
      # Text elements
      text = element_text(family = base_family),
      plot.title = element_text(
        hjust = 0.5, 
        size = 22,
        margin = margin(b = 5)
      ),
      plot.subtitle = element_text(
        hjust = 0.5, 
        size = 12,
        margin = margin(b = 5)
      ),
      axis.title = element_text(size = 12),
      axis.text = element_text(size = 12),
      axis.text.y = element_text(margin = margin(r = 0), hjust = 1),
      
      # Remove legend
      legend.position = "none",
      
      # Grid lines
      panel.grid.major.x = element_blank(),
      panel.grid.minor = element_blank(),
      
      # Margins
      plot.margin = margin(t = 0, r = 0, b = 0, l = 0)
    )
}

# --- Original Values ---
boxplot_data <- dRNA_adjusted %>%
  select(avg_mm_WT, avg_mm_KD) %>%
  pivot_longer(
    cols = c(avg_mm_WT, avg_mm_KD),
    names_to = "category",
    values_to = "mismatch_rate"
  ) %>%
  mutate(category = factor(category, levels = c("avg_mm_WT", "avg_mm_KD")))

color_map <- c(
  "avg_mm_WT" = dRNA_WT,
  "avg_mm_KD" = dRNA_KD
)

# label outliers
y_limit <- 100
arrow_data <- boxplot_data %>%
  group_by(category) %>%
  summarise(has_outliers = any(mismatch_rate > y_limit)) %>%
  filter(has_outliers == TRUE)

p <- ggplot(boxplot_data, aes(x = category, y = mismatch_rate, fill = category)) +
  geom_boxplot() +
  scale_fill_manual(values = color_map) +
  scale_x_discrete(labels = c("avg_mm_WT" = "WT", "avg_mm_KD" = "KD")) +
  geom_point(
    data = arrow_data,
    aes(x = category, y = y_limit), # Place the arrow at the top of the plot
    shape = 24, # An upward-pointing triangle
    size = 2,
    fill = dRNA_WT
  ) +
  coord_cartesian(ylim = c(0, y_limit)) +
  labs(
    title = "Comparison of Mismatch Rates",
    subtitle = "Original Rate",
    x = NULL,
    y = "U-to-C Mismatch (%)"
  ) +
  theme_boxplot() +
  theme(legend.position = "none")

ggsave("plots/original_dRNA.png", p, width = 4, height = 4, unit = "in")
ggsave("plots/original_dRNA.pdf", p, width = 1.8, height = 2, unit = "in")

# --- Corrected Values ---

boxplot_data <- dRNA_adjusted %>%
  select(adj_mm_WT, adj_mm_KD) %>%
  pivot_longer(
    cols = c(adj_mm_WT, adj_mm_KD),
    names_to = "category",
    values_to = "mismatch_rate"
  ) %>%
  mutate(category = factor(category, levels = c("adj_mm_WT", "adj_mm_KD")))

color_map <- c(
  "adj_mm_WT" = dRNA_adj_WT,
  "adj_mm_KD" = dRNA_adj_KD
)

# label outliers
y_limit <- 100
arrow_data <- boxplot_data %>%
  group_by(category) %>%
  summarise(has_outliers = any(mismatch_rate > y_limit)) %>%
  filter(has_outliers == TRUE)

p <- ggplot(boxplot_data, aes(x = category, y = mismatch_rate, fill = category)) +
  geom_boxplot() +
  scale_fill_manual(values = color_map) +
  scale_x_discrete(labels = c("adj_mm_WT" = "WT", "adj_mm_KD" = "KD")) +
  geom_point(
    data = arrow_data,
    aes(x = category, y = y_limit), 
    shape = 24, 
    size = 2,
    fill = dRNA_adj_WT
  ) +
  coord_cartesian(ylim = c(0, y_limit)) +
  labs(
    title = "Comparison of Mismatch Rates",
    subtitle = "Corrected Rate with Subtraction",
    x = NULL,
    y = "Psi Stoichiometry"
  ) +
  theme_boxplot() +
  theme(legend.position = "none")

ggsave("plots/corrected_dRNA_subtract.png", p, width = 4, height = 4, unit = "in")
ggsave("plots/corrected_dRNA_subtract.pdf", p, width = 1.8, height = 2, unit = "in")

# --- Heatmap ---

theme_heat <- function(base_size = 18, base_family = "sans") {
  theme_minimal(base_size = base_size) %+replace%
    theme(
      # Text elements
      text = element_text(family = base_family),
      plot.title = element_text(
        hjust = 0.5, 
        size = base_size,
        margin = margin(b = 5)
      ),
      plot.subtitle = element_text(
        hjust = 0.5, 
        size = 12,
        margin = margin(b = 5)
      ),
      axis.title = element_text(size = 12),
      axis.text = element_text(size = 12),
      axis.text.y = element_text(margin = margin(r = 0), hjust = 1),
      
      legend.position = "right",
      legend.text = element_text(size = 10), 
      legend.title = element_text(size = 12),
      
      # Grid lines
      panel.grid.major.y = element_blank(),
      panel.grid.minor = element_blank(),
      
      panel.grid = element_blank(),
      
      # Margins
      plot.margin = margin(t = 0, r = 0, b = 0, l = 0)
    )
}

ordered_levels <- dRNA_adjusted %>%
  arrange(adj_mm_WT) %>%
  pull(gene_pos_motif) 

heatmap_data <- dRNA_adjusted %>%
  pivot_longer(
    cols = c(avg_mm_WT, avg_mm_KD, adj_mm_WT, adj_mm_KD),
    names_to = "category",
    values_to = "mismatch_rate"
  ) %>%
  mutate(category = factor(category, levels = c("avg_mm_WT", "avg_mm_KD", "adj_mm_WT", "adj_mm_KD")))

p_heatmap <- ggplot(heatmap_data %>%
                      mutate(gene_pos_motif = factor(gene_pos_motif, levels = ordered_levels)),
                    aes(x = category, y = gene_pos_motif, fill = mismatch_rate)) +
  geom_tile() +
  scale_fill_gradient(
    low = "white",
    high = dRNA_adj_WT,
    name = "(%)",
    limits = c(0, 100),
    oob = scales::squish # Forces values >100 to the max color
  ) +
  scale_x_discrete(labels = c("avg_mm_WT" = "WT", "avg_mm_KD" = "KD", "adj_mm_WT" = "WT", "adj_mm_KD" = "KD")) +
  labs(
    title = "Comparison of Mismatch Rates",
    x = NULL,
    y = "Sites"
  ) +
  theme_heat() +
  theme(axis.text.y = element_text(size = 8)) # Adjust size based on number of sites

print(p_heatmap)

# Save the heatmap with dimensions that suit a potentially tall plot
ggsave("plots/heatmap_dRNA.png", p_heatmap, width = 4, height = 8, unit = "in")
ggsave("plots/heatmap_dRNA.pdf", p_heatmap, width = 5.5, height = 6, unit = "in")
