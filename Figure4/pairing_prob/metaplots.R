#!/usr/bin/env Rscript

library(dplyr)
library(readr)
library(tidyr)
library(ggplot2)
library(stringr)
library(RColorBrewer)

incell <- "#e25098"
invitro <- "#a2a2fc"
dRNA <- "#daa520"

# ---- Functions ----
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
theme_line <- function(base_size = 16) {
  theme_base_custom(base_size) +
    theme(legend.position = "right", panel.grid.minor = element_blank())
}

align_probabilities <- function(prob_df, info_df, psipos_col_name) {
  # Identify the raw probability data columns (e.g., X1, X2, ...)
  prob_cols <- grep("^X[0-9]+$", names(prob_df), value = TRUE)
  data_merged <- prob_df %>%
    left_join(info_df, by = c("rna_name" = "clean_name"))
  # Apply the alignment logic to each row individually
  aligned_matrix <- t(apply(data_merged, 1, function(row) {
    # Extract the center position for this row
    psipos <- as.numeric(row[psipos_col_name])
    if (is.na(psipos)) {
      return(rep(NA_real_, 121)) # Return a vector of NAs of the correct length
    }
    # Define the 101-position window and create a blank result vector
    start_pos <- psipos - 60
    end_pos <- psipos + 60
    target_indices <- start_pos:end_pos
    result_window <- rep(NA_real_, 121)
    # Identify which positions in the window are valid (i.e., not off the edge of the data)
    valid_indices_mask <- target_indices >= 1 & target_indices <= length(prob_cols)
    # Extract the valid probability values from the original data
    prob_values <- as.numeric(row[prob_cols])
    values_to_place <- prob_values[target_indices[valid_indices_mask]]
    # Place these values into the correct spots in the result window
    result_window[valid_indices_mask] <- values_to_place
    return(result_window)
  }))
  # Convert the result matrix to a clean, named dataframe
  aligned_df <- as.data.frame(aligned_matrix)
  positions <- -60:60
  position_labels <- ifelse(positions < 0, paste0("m", abs(positions)), as.character(positions))
  names(aligned_df) <- paste0("pos_", position_labels)
  aligned_df <- bind_cols(rna_name = prob_df$rna_name, aligned_df)
  return(aligned_df)
}

create_multi_line_plot <- function(data_list,
                                   title,
                                   filename,
                                   subtitle = NULL,
                                   window_size = 20,
                                   error_metric = "se",
                                   palette_name = "Set1"
                                   color_map = NULL,
                                   max_legend_chars_per_row = 45,
                                   legend_label_width = 30,
                                   plot_width = 2, plot_height = 3) {
  # --- 1. Set up Colors and Legend Layout ---
  source_names <- names(data_list)
  plot_colors <- color_map
  
  if (is.null(plot_colors)) {
    n_colors <- length(source_names)
    max_brewer_colors <- RColorBrewer::brewer.pal.info[palette_name, "maxcolors"] # nolint: line_length_linter.
    if (n_colors > max_brewer_colors) {
      get_palette <- colorRampPalette(brewer.pal(max_brewer_colors, palette_name)) # nolint: line_length_linter.
      colors_to_use <- get_palette(n_colors)
    } else {
      colors_to_use <- brewer.pal(n_colors, palette_name)
    }
    plot_colors <- setNames(colors_to_use, source_names)
  }
  
  # Automatically calculate the number of legend rows based on total label length
  overhead_per_item <- 5 
  total_legend_chars <- sum(nchar(source_names)) + (length(source_names) * overhead_per_item)
  legend_rows <- ceiling(total_legend_chars / max_legend_chars_per_row)
  
  # --- 2. Summarize Data ---
  summary_fn <- if (error_metric == "sd") {
    list(mean = ~mean(.x, na.rm = TRUE), error = ~sd(.x, na.rm = TRUE))
  } else {
    list(mean = ~mean(.x, na.rm = TRUE), error = ~sd(.x, na.rm = TRUE) / sqrt(sum(!is.na(.x))))
  }
  
  summary_values <- bind_rows(data_list, .id = "source") %>%
    mutate(source = factor(source, levels = source_names)) %>%
    group_by(source) %>%
    summarise(across(starts_with("pos_"), .fns = summary_fn), .groups = "drop") %>%
    pivot_longer(
      cols = -source,
      names_to = c("Position", ".value"),
      names_pattern = "pos_(.*)_(mean|error)"
    ) %>%
    mutate(Position = as.numeric(sub("m", "-", Position))) %>%
    filter(between(Position, -window_size, window_size)) %>%
    arrange(source, Position)
  
  # --- 3. Create the Plot ---
  wrapped_labels <- str_wrap(source_names, width = legend_label_width)
  
  plot <- ggplot(summary_values, aes(x = Position, y = mean, color = source, fill = source)) +
    geom_ribbon(aes(ymin = mean - error, ymax = mean + error), alpha = 0.25, linetype = 0) +
    geom_line(linewidth = 1) +
    geom_point(size = 0.5) +
    geom_point(data = . %>% filter(Position == 0), size = 1.5, color = "#c279ca", show.legend = FALSE) +
    scale_color_manual(values = plot_colors, labels = wrapped_labels) +
    scale_fill_manual(values = plot_colors, labels = wrapped_labels) +
    labs(
      title = title,
      subtitle = subtitle,
      x = "Position",
      y = "Average Pairing Probability",
      color = "Data Source",
      fill = "Data Source"
    ) +
    guides(color = guide_legend(nrow = legend_rows), fill = guide_legend(nrow = legend_rows)) +
    scale_x_continuous(
      breaks = seq(-window_size, window_size, by = ifelse(window_size > 20, 10, 5)),
      limits = c(-window_size, window_size),
      expand = c(0, 0)
    ) +
    scale_y_continuous(breaks = seq(0, 1, 0.2), expand = c(0, 0)) +
    coord_cartesian(ylim = c(0, 1)) +
    theme_line(base_size = 10) +
    theme(
      aspect.ratio = 1,
      panel.border = element_rect(color = "black", fill = NA, linewidth = 1),
      legend.position = "bottom",
      legend.title = element_blank(),
      plot.margin = margin(10, 10, 10, 10)
    )
  
  # --- 4. Saving ---
  print(paste0(plot_dir, filename, ".pdf"))
  ggsave(filename = paste0(plot_dir, filename, ".pdf"), plot = plot, width = plot_width, height = plot_height)
  ggsave(filename = paste0(plot_dir, filename, ".png"), plot = plot, width = plot_width, height = plot_height, dpi = 600)
  
  print(plot)
  return(plot)
}

# ---- Datasets ----

ivSHAPE_pairprob_pool1 <- read.csv("IV_NoPus_pairprob_aligned.csv")
rnafold_pairprob_pool1 <- read.csv("pool1_noadapters_rnafold_pairprob_aligned.csv")
icSHAPE_pairprob_pool1 <- read.csv("HEK293T_pairprob_aligned.csv")

iv_randomU_background <- read.csv("preprint_IV_pairprob_background.csv")
ic_randomU_background <- read.csv("IC_pairprob_background.csv")

dRNA_RNAfold <- read.csv("dRNA_rnafold_pairprob_aligned.csv")
dRNA_background_RNAfold <- read.csv("dRNA_background_downsamp_rnafold_pairprob_aligned.csv")

dRNA_sites <- read.csv("~/PUS7regulation2026/Figure1/UNUARstandard/PUS7_10reads_15percent_motif.csv")
pool1 <- read.csv("~/PUS7regulation2026/Figure4/sequence/pool1_background.csv")
dRNA_unuar <- dRNA_sites %>%
  left_join(pool1 %>%
              dplyr::select(clean_name, chr, psi),
            by = c("chr" = "chr", "pos" = "psi"))

invitro_full <- read.csv("~/PUS7regulation2026/Figure4/sequence/invitroBID_MPRA_summary.csv", header = TRUE)
invitro_PUS7 <- invitro_full %>% filter(PUS7_dependent == "Yes") %>% filter(!(is.na(psipos_motif)))

incell_full <- read.csv("~/PUS7regulation2026/Figure4/sequence/incellBID_WTsummary.csv", header = TRUE)
incell_PUS7 <- incell_full %>% filter(PUS7_dependent == "Yes") %>% filter(!(is.na(psipos_motif)))

icPUS7_BothQuartiles <- read_tsv("~/PUS7regulation2026/Figure4/sequence/incell_quartiles/PUS7_dep_union_by_Both_WT_quartiles.tsv", col_names = TRUE)
ivPUS7_Q1 <- read_tsv("~/PUS7regulation2026/Figure4/sequence/invitro_quartiles/invitro_modified_sites_Q1_summary.tsv", col_names = TRUE)
ivPUS7_Q2 <- read_tsv("~/PUS7regulation2026/Figure4/sequence/invitro_quartiles/invitro_modified_sites_Q2_summary.tsv", col_names = TRUE)
ivPUS7_Q3 <- read_tsv("~/PUS7regulation2026/Figure4/sequence/invitro_quartiles/invitro_modified_sites_Q3_summary.tsv", col_names = TRUE)
ivPUS7_Q4 <- read_tsv("~/PUS7regulation2026/Figure4/sequence/invitro_quartiles/invitro_modified_sites_Q4_summary.tsv", col_names = TRUE)

plot_dir <- "~/PUS7regulation2026/Figure4/plots"

# ---- Plots - in vitro SHAPE ----

# direct RNA sites with in vitro data
dRNA_data <- list(
    "Background Sites" = iv_randomU_background,
    "Direct RNA Sites" = ivSHAPE_pairprob_pool1 %>%
        filter(rna_name %in% dRNA_unuar$clean_name))
dRNA_colors <- c(
    "Direct RNA Sites" = dRNA,
    "Background Sites" = "darkgray")
create_multi_line_plot(data_list = dRNA_data, title = "Direct RNA", filename = "/4_metaplot_dRNA_SHAPE",
    window_size = 20, color_map = dRNA_colors)

# in vitro PUS7 sites with in vitro data
iv_data <- list(
    "Background Sites" = iv_randomU_background,
    "In Vitro PUS7 Sites" = ivSHAPE_pairprob_pool1 %>%
        filter(rna_name %in% invitro_PUS7$chr))
iv_colors <- c(
    "In Vitro PUS7 Sites" = invitro,
    "Background Sites" = "darkgray")
create_multi_line_plot(data_list = iv_data, title = "In Vitro", filename = "/4_metaplot_ivPUS7_SHAPE",
    window_size = 20, color_map = iv_colors)

# in cell PUS7 sites with in vitro data
ic_data <- list(
    "Background Sites" = iv_randomU_background,
    "In Cell PUS7 Sites" = ivSHAPE_pairprob_pool1 %>%
        filter(rna_name %in% incell_PUS7$chr))
ic_colors <- c(
    "In Cell PUS7 Sites" = incell,
    "Background Sites" = "darkgray")
create_multi_line_plot(data_list = ic_data, title = "In Cell", filename = "/4_metaplot_icPUS7_SHAPE",
    window_size = 20, color_map = ic_colors)

# in vitro MPRA sites by quartile

IV_pairprob_Q1 <- ivSHAPE_pairprob_pool1 %>%
  filter(rna_name %in% ivPUS7_Q1$chr)
IV_pairprob_Q2 <- ivSHAPE_pairprob_pool1 %>%
  filter(rna_name %in% ivPUS7_Q2$chr)
IV_pairprob_Q3 <- ivSHAPE_pairprob_pool1 %>%
  filter(rna_name %in% ivPUS7_Q3$chr)
IV_pairprob_Q4 <- ivSHAPE_pairprob_pool1 %>%
  filter(rna_name %in% ivPUS7_Q4$chr)

ivQ_colors <- c(
    "Q1 In Vitro PUS7, In Vitro SHAPE" = "#4a4a85",
    "Q2 In Vitro PUS7, In Vitro SHAPE" = "#7b9dbe",
    "Q3 In Vitro PUS7, In Vitro SHAPE" = "#b5cadf",
    "Q4 In Vitro PUS7, In Vitro SHAPE" = invitro
)
invitro_data <- list(
    "Q2 In Vitro PUS7, In Vitro SHAPE" = IV_pairprob_Q2,
    "Q3 In Vitro PUS7, In Vitro SHAPE" = IV_pairprob_Q3,
    "Q1 In Vitro PUS7, In Vitro SHAPE" = IV_pairprob_Q1,
    "Q4 In Vitro PUS7, In Vitro SHAPE" = IV_pairprob_Q4
)
create_multi_line_plot(
    data_list = invitro_data,
    title = "In Vitro Quartiles SHAPE",
    filename = "/S4_metaplot_iv_quartiles",
    window_size = 20,
    color_map = ivQ_colors,
    error_metric = "se",
    plot_width = 3.5,
    plot_height = 4.5
)

# in cell MPRA sites by quartile
IC_pairprob_Q1 <- ivSHAPE_pairprob_pool1 %>%
  inner_join(icPUS7_BothQuartiles, by = c("rna_name" = "chr")) %>% filter(quartile == 1)
IC_pairprob_Q2 <- ivSHAPE_pairprob_pool1 %>%
  inner_join(icPUS7_BothQuartiles, by = c("rna_name" = "chr")) %>% filter(quartile == 2)
IC_pairprob_Q3 <- ivSHAPE_pairprob_pool1 %>%
  inner_join(icPUS7_BothQuartiles, by = c("rna_name" = "chr")) %>% filter(quartile == 3)
IC_pairprob_Q4 <- ivSHAPE_pairprob_pool1 %>%
  inner_join(icPUS7_BothQuartiles, by = c("rna_name" = "chr")) %>% filter(quartile == 4)

icQ_colors <- c(
    "Q1 In Cell PUS7, In Vitro SHAPE" = "#991e5c",
    "Q2 In Cell PUS7, In Vitro SHAPE" = "#997e7e",
    "Q3 In Cell PUS7, In Vitro SHAPE" = "#f5a9a9",
    "Q4 In Cell PUS7, In Vitro SHAPE" = incell
)
incell_data <- list(
    "Q2 In Cell PUS7, In Vitro SHAPE" = IC_pairprob_Q2,
    "Q3 In Cell PUS7, In Vitro SHAPE" = IC_pairprob_Q3,
    "Q1 In Cell PUS7, In Vitro SHAPE" = IC_pairprob_Q1,
    "Q4 In Cell PUS7, In Vitro SHAPE" = IC_pairprob_Q4
)
create_multi_line_plot(
    data_list = incell_data,
    title = "In Cellulo Quartiles SHAPE",
    filename = "/S4_metaplot_ic_quartiles",
    window_size = 20,
    color_map = icQ_colors,
    error_metric = "se",
    plot_width = 3.5,
    plot_height = 4.5
)

# ---- Plots - RNAfold Only ----

# direct RNA sites
dRNA_data <- list(
    "Background Sites" = dRNA_background_RNAfold,
    "Direct RNA PUS7 Sites" = dRNA_RNAfold)
dRNA_colors <- c(
    "Direct RNA PUS7 Sites" = dRNA,
    "Background Sites" = "darkgray")
create_multi_line_plot(data_list = dRNA_data, title = "Direct RNA", filename = "/S4_metaplot_dRNA_RNAfold",
    window_size = 20, color_map = dRNA_colors)
# in vitro PUS7 sites
iv_data <- list(
    "Background Sites" = rnafold_pairprob_pool1 %>%
        filter(!(rna_name %in% invitro_PUS7$chr)),
    "In Vitro PUS7 Sites" = rnafold_pairprob_pool1 %>%
        filter(rna_name %in% invitro_PUS7$chr))
iv_colors <- c(
    "In Vitro PUS7 Sites" = invitro,
    "Background Sites" = "darkgray")
create_multi_line_plot(data_list = iv_data, title = "In Vitro", filename = "/S4_metaplot_ivPUS7_RNAfold",
    window_size = 20, color_map = iv_colors)

# in cell PUS7 sites 
ic_data <- list(
    "Background Sites" = rnafold_pairprob_pool1 %>%
        filter(!(rna_name %in% incell_PUS7$chr)),
    "In Cell PUS7 Sites" = rnafold_pairprob_pool1 %>%
        filter(rna_name %in% incell_PUS7$chr))
ic_colors <- c(
    "In Cell PUS7 Sites" = incell,
    "Background Sites" = "darkgray")
create_multi_line_plot(data_list = ic_data, title = "In Cell", filename = "/S4_metaplot_icPUS7_RNAfold",
    window_size = 20, color_map = ic_colors)

# ---- Plots - in cell SHAPE ----

# direct RNA sites with in cell data
dRNA_data <- list(
    "Background Sites" = ic_randomU_background,
    "Direct RNA Sites" = icSHAPE_pairprob_pool1 %>%
        filter(rna_name %in% dRNA_unuar$clean_name))
dRNA_colors <- c(
    "Direct RNA Sites" = dRNA,
    "Background Sites" = "darkgray")
create_multi_line_plot(data_list = dRNA_data, title = "Direct RNA", filename = "/S4_metaplot_dRNA_icSHAPE",
    window_size = 20, color_map = dRNA_colors)

# in cell PUS7 sites with in cell data
ic_data <- list(
    "Background Sites" = ic_randomU_background,
    "In Cell PUS7 Sites" = icSHAPE_pairprob_pool1 %>%
        filter(rna_name %in% incell_PUS7$chr))
ic_colors <- c(
    "In Cell PUS7 Sites" = incell,
    "Background Sites" = "darkgray")
create_multi_line_plot(data_list = ic_data, title = "In Cell", filename = "/S4_metaplot_icPUS7_icSHAPE",
    window_size = 20, color_map = ic_colors)