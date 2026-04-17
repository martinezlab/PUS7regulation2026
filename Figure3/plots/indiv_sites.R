#!/usr/bin/env Rscript
library(dplyr)
library(stringr)
library(readr)
library(ggplot2)
library(scales)
library(tidyr)
library(purrr)

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
theme_heat <- function(base_size = 16) {
  theme_base_custom(base_size) +
    theme(legend.position = "right", panel.grid = element_blank())
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
invitro <- "#a2a2fc"
both <- "#c279ca"
input <- "#eee8aa"

KD <- "#ed93c0"
OE <- "#a21b5d"

# ---- Data ----
invitro_all <- read_tsv("invitro_standardized_input.tsv")
incell_all <- read_tsv("BIDdetect_data_incell_delpos.txt")

# ---- Plots ----

create_deletion_plot_facet <- function(data, sites, filename, experiment_type,
                                       plot_width = 7, plot_height = 2.5, cellline = c("293T", "HepG2")) {
  # Process all sites
  all_data <- map_df(sites, ~process_site(data, ., experiment_type, cellline))
  # Find the maximum delrate value
  max_delrate <- max(all_data$delrate, na.rm = TRUE)
  # Set up color scales and x-axis limits based on experiment type
  if (experiment_type == "incell") {
    color_scale <- scale_fill_manual(values = c(
      "in" = input,
      "KD" = KD,
      "WT" = incell,
      "OE" = OE
    ))
    x_limits <- c("in", "KD", "WT", "OE")
  } else if (experiment_type == "invitro") {
    color_scale <- scale_fill_manual(values = c(
      "in" = input,
      "BS" = invitro
    ))
    x_limits <- c("in", "BS")
  } else if (experiment_type == "endo") {
    color_scale <- scale_fill_manual(values = c(
      "in" = input,
      "KD" = KD,
      "WT" = incell
    ))
    x_limits <- c("in", "KD", "WT")
  } else if (experiment_type == "all") {
    color_scale <- scale_fill_manual(values = c(
      "in" = input, 
      "KD" = KD, 
      "WT" = incell, 
      "OE" = OE,
      "IV" = invitro
    ))
    x_limits <- c("in", "KD", "WT", "OE", "IV")
  }
  # Create faceted plot
  facet_plot <- ggplot(all_data, aes(x = expression, y = delrate)) +
    stat_summary(fun = mean, geom = "bar", aes(fill = expression)) +
    stat_summary(fun.data = mean_se, geom = "errorbar", width = 0.2) +
    geom_jitter(width = 0.2, alpha = 0.5) +
    labs(
      x = "",
      y = "Deletion Rate"
    ) +
    theme_bar(base_size = 12) +
    color_scale +
    scale_x_discrete(limits = x_limits) +
    scale_y_continuous(limits = c(0, max_delrate), expand = expansion(mult = c(0, 0.1))) +
    facet_wrap(~ gene_info + chr_info, scales = "free_y", ncol = 5) +
    theme(strip.text = element_text(size = 12))
  ggsave(filename, facet_plot, width = plot_width, height = plot_height, units = "in")
  return(facet_plot)
}


# In vitro specific sites
mpra_sites<- c("HDAC6_chrX_48823353", "RHBDD2_chr7_75888787", "STIM1_chr11_4092507")
pdf_path <- file.path(plot_dir, "3_barplot_iv_mpra_sites.pdf")
create_deletion_plot_facet(invitro_all, mpra_sites, pdf_path, "invitro", 
                            plot_width = 5, plot_height = 2.25)
png_path <- file.path(plot_dir, "3_barplot_iv_mpra_sites.png")
create_deletion_plot_facet(invitro_all, mpra_sites, png_path, "invitro", 
                            plot_width = 5, plot_height = 2.25)

# In cell specific sites, 293T only
pdf_path <- file.path(plot_dir, "3_barplot_ic_mpra_sites_293T.pdf")
create_deletion_plot_facet(incell_all, mpra_sites, pdf_path, "incell", 
                            plot_width = 5, plot_height = 2.25, cellline = "293T")
png_path <- file.path(plot_dir, "3_barplot_ic_mpra_sites_293T.png")
create_deletion_plot_facet(incell_all, mpra_sites, png_path, "incell", 
                            plot_width = 5, plot_height = 2.25, cellline = "293T")

# In cell KD only sites, 293T only
selected_sites_KDonly <- c("ACTL6A_chr3_179587991", "FKBP1B_chr2_24063441", "EML4_chr2_42329888")
pdf_path <- file.path(plot_dir, "S3_barplot_KDonly_sites.pdf")
create_deletion_plot_facet(all_data, selected_sites_KDonly, pdf_path, "all", 
                            plot_width = 5, plot_height = 2, cellline = "293T")
png_path <- file.path(plot_dir, "S3_barplot_KDonly_sites.png")
create_deletion_plot_facet(all_data, selected_sites_KDonly, png_path, "all", 
                            plot_width = 5, plot_height = 2, cellline = "293T")

# In cell OE only sites, 293T only
selected_sites_OEonly <- c("TMEM109_chr11_60923252", "TGFBR2_chr3_30612922", "CCT7_chr2_73252776")
pdf_path <- file.path(plot_dir, "S3_barplot_OEonly_sites.pdf")
create_deletion_plot_facet(all_data, selected_sites_OEonly, pdf_path, "all", 
                            plot_width = 5, plot_height = 2, cellline = "293T")
png_path <- file.path(plot_dir, "S3_barplot_OEonly_sites.png")
create_deletion_plot_facet(all_data, selected_sites_OEonly, png_path, "all", 
                            plot_width = 5, plot_height = 2, cellline = "293T")