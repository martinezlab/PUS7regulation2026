#!/usr/bin/env Rscript

library(dplyr)
library(readr)
library(ggplot2)
library(tidyr)

# ---- Directory ----
active_dir <- "~/PUS7regulation2026/Figure6/read_counts"

# ---- Plot themes ----
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
theme_scatter <- function(base_size = 16, base_family = "sans") {
  theme_base_custom(base_size) +
  theme(legend.position = "none",
        panel.grid.minor = element_blank(),
        panel.border = element_rect(color = "black", fill = NA, size = 1),
        plot.margin = margin(t = 5, r = 5, b = 5, l = 5))
}

# ---- Colors ----
incell <- "#e25098"
cellline <- "#00a86b"

# ---- Datasets ----

incell_counts <- read_tsv(file.path(active_dir, "incell_read_counts.tsv"))

incell_all  <- read_tsv("~/PUS7regulation2026/Figure3/plots/BIDdetect_data_incell_delpos.txt")
incell_PUS7 <- read_tsv("~/PUS7regulation2026/Figure3/plots/PUS7_dep_union_significant_summary.tsv")
incell_sum  <- read_tsv("~/PUS7regulation2026/Figure3/plots/WT_mod_Both_significance.tsv")
CTspec_PUS7 <- read_tsv(file.path(active_dir, "celltype_spec_PUS7_union_293T_high_summary.tsv"))

# ---- Normalize Read Counts ----

incell_norm <- incell_all %>%
  mutate(Filename = paste0(celltype, "_", vector, "_", rep, "_", treat, ".bam")) %>%
  left_join(incell_counts, by = "Filename") %>%
  mutate(norm_reads = totalReads / Read_Count) %>%
  select(celltype, vector, rep, treat, chr, delrate, norm_reads)

hepg2_wt_reads <- incell_norm %>%
  filter(vector %in% c("P102", "P4", "WT"), celltype == "HepG2") %>%
  group_by(chr) %>%
  summarize(avg_reads = mean(norm_reads))
  
hepg2_wt_delrate <- incell_all %>% 
  dplyr::select(chr, celltype, vector, treat, rep, delrate) %>%
  filter(celltype == "HepG2", vector %in% c("P102", "P4", "WT")) %>%
  group_by(chr, treat) %>%
  summarize(avg_delrate = mean(delrate)) %>%
  pivot_wider(names_from = treat, values_from = avg_delrate) %>%
  mutate(delta_delrate = BS - input)

hepg2_wt <- hepg2_wt_delrate %>%
  left_join(hepg2_wt_reads, by = "chr") %>%
  mutate(highlight_group = case_when(
    chr %in% CTspec_PUS7$chr ~ "PUS7-dependent, cell-type specific",
    chr %in% incell_PUS7$chr ~ "PUS7-dependent",
    TRUE ~ "Not PUS7-dependent"
  ))

hek293t_wt_reads <- incell_norm %>%
  filter(vector %in% c("P102", "P4", "WT"), celltype == "293T") %>%
  group_by(chr) %>%
  summarize(avg_reads = mean(norm_reads))
  
hek293t_wt_delrate <- incell_all %>% 
  dplyr::select(chr, celltype, vector, treat, rep, delrate) %>%
  filter(celltype == "293T", vector %in% c("P102", "P4", "WT")) %>%
  group_by(chr, treat) %>%
  summarize(avg_delrate = mean(delrate)) %>%
  pivot_wider(names_from = treat, values_from = avg_delrate) %>%
  mutate(delta_delrate = BS - input)

hek293t_wt <- hek293t_wt_delrate %>%
  left_join(hek293t_wt_reads, by = "chr") %>%
  mutate(highlight_group = case_when(
    chr %in% CTspec_PUS7$chr ~ "PUS7-dependent, cell-type specific",
    chr %in% incell_PUS7$chr ~ "PUS7-dependent",
    TRUE ~ "Not PUS7-dependent"
  ))

cell_line_reads <- hek293t_wt %>%
  inner_join(hepg2_wt, by = c("chr", "highlight_group"), suffix = c("_293t", "_hepg2")) %>%
  dplyr::select(chr, avg_reads_293t, avg_reads_hepg2, delta_delrate_293t, delta_delrate_hepg2, highlight_group)


# --- Plot Norm Read Counts ---

# Calculate Pearson and Spearman correlations
pearson_cor <- cor(cell_line_reads$avg_reads_293t, cell_line_reads$avg_reads_hepg2, method = "pearson", use = "complete.obs")

cor_label <- paste0("Pearson r = ", round(pearson_cor, 3))

p <- ggplot(cell_line_reads, aes(x = avg_reads_hepg2, y = avg_reads_293t, color = highlight_group)) +
  geom_point(size = 0.5, alpha = 0.6) +
  annotate("text", x = Inf, y = Inf, label = cor_label, hjust = 1.1, vjust = 1.5, color = "black", size = 2) +
  labs(
    title = "In Cellulo",
    x = "HepG2 WT Normalized Read Counts",
    y = "HEK293T WT Normalized Read Counts"
  ) +
  coord_cartesian(xlim = c(0, 0.02), ylim = c(0, 0.02)) +
  theme_scatter(base_size = 8) +
  theme(legend.position = "bottom") +
  scale_color_manual(
    values = c("Not PUS7-dependent" = "#707070", 
               "PUS7-dependent" = incell,
               "PUS7-dependent, cell-type specific" = cellline),
    name = NULL,
    guide = "legend"
  )

plot_path <- paste0(active_dir, "/scatter_CTnormreads")
ggsave(paste0(plot_path, ".png"), p, width = 2, height = 2.5, unit = "in")
ggsave(paste0(plot_path, ".pdf"), p, width = 2, height = 2.5, unit = "in")

# ---- Plot Norm Counts v Deletion Fraction ----

# -- HepG2 --
# Calculate Pearson and Spearman correlations
pearson_cor <- cor(cell_line_reads$delta_delrate_hepg2, cell_line_reads$avg_reads_hepg2, method = "pearson", use = "complete.obs")

cor_label <- paste0("Pearson r = ", round(pearson_cor, 3))

p <- ggplot(cell_line_reads, aes(x = avg_reads_hepg2, y = delta_delrate_hepg2, color = highlight_group)) +
  geom_point(size = 0.5, alpha = 0.6) +
  annotate("text", x = Inf, y = Inf, label = cor_label, hjust = 1.1, vjust = 1.5, color = "black", size = 2) +
  labs(
    title = "HepG2 WT",
    x = "Normalized Read Counts",
    y = "Delta Deletion Fraction"
  ) +
  theme_scatter(base_size = 8) +
  theme(legend.position = "bottom") +
  coord_cartesian(xlim = c(0, 0.015), ylim = c(0, 1)) +
  scale_color_manual(
    values = c("Not PUS7-dependent" = "#707070", 
               "PUS7-dependent" = incell,
               "PUS7-dependent, cell-type specific" = cellline),
    name = NULL,
    guide = "legend"
  )

plot_path <- paste0(active_dir, "/scatter_HepG2normreadvdelrate")
ggsave(paste0(plot_path, ".png"), p, width = 2, height = 2.5, unit = "in")
ggsave(paste0(plot_path, ".pdf"), p, width = 2, height = 2.5, unit = "in")

# -- HEK293T --
# Calculate Pearson and Spearman correlations
pearson_cor <- cor(cell_line_reads$delta_delrate_293t, cell_line_reads$avg_reads_293t, method = "pearson", use = "complete.obs")

cor_label <- paste0("Pearson r = ", round(pearson_cor, 3))

p <- ggplot(cell_line_reads, aes(x = avg_reads_293t, y = delta_delrate_293t, color = highlight_group)) +
  geom_point(size = 0.5, alpha = 0.6) +
  annotate("text", x = Inf, y = Inf, label = cor_label, hjust = 1.1, vjust = 1.5, color = "black", size = 2) +
  labs(
    title = "HEK293T WT",
    x = "Normalized Read Counts",
    y = "Delta Deletion Fraction"
  ) +
  theme_scatter(base_size = 8) +
  theme(legend.position = "bottom") +
  coord_cartesian(xlim = c(0, 0.015), ylim = c(0, 1)) +
  scale_color_manual(
    values = c("Not PUS7-dependent" = "#707070", 
               "PUS7-dependent" = incell,
               "PUS7-dependent, cell-type specific" = cellline),
    name = NULL,
    guide = "legend"
  )

plot_path <- paste0(active_dir, "/scatter_HEK293Tnormreadvdelrate")
ggsave(paste0(plot_path, ".png"), p, width = 2, height = 2.5, unit = "in")
ggsave(paste0(plot_path, ".pdf"), p, width = 2, height = 2.5, unit = "in")