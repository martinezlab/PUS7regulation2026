#!/usr/bin/env Rscript

library(dplyr)
library(readr)
library(tidyr)
library(ggplot2)
library(ggrepel)
library(scales)
library(tidyr)

# ---- Directory ----
active_dir <- "~/PUS7regulation2026/Figure1/KO"

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

# ---- Load data ----
borderline <- read.csv(file.path(active_dir, "borderline/dRNA_borderline_sites.csv"), header = TRUE)
PUS7dep <- read.csv(file.path(active_dir, "UNUARstandard/PUS7_10reads_15percent_motif.csv"), header = TRUE)
KOdata <- read_tsv(file.path(active_dir, "KO/HepG2_PUS7KO_target_counts.csv"))

# ---- Clean & math on dataframes ----
KOdata <- KOdata %>%
  dplyr::mutate(KO_mm_rate = (C.count / (T.count + C.count)) * 100) %>%
  select(!(delrate))

borderline_clean <- borderline %>%
  mutate(
    WT_mm_avg = rowMeans(cbind(WT_mm_rep2, WT_mm_rep3, WT_mm_rep4), na.rm = TRUE),
    KD_mm_avg = rowMeans(cbind(KD_mm_rep2, KD_mm_rep3, KD_mm_rep4), na.rm = TRUE)
  )

PUS7dep_clean <- PUS7dep %>%
  rename(
    WT_mm_rep2 = mm.perc.rep2, WT_mm_rep3 = mm.perc.rep3, WT_mm_rep4 = mm.perc.rep4,
    KD_mm_rep2 = ctrl.err.rep2, KD_mm_rep3 = ctrl.err.rep3, KD_mm_rep4 = ctrl.err.rep4) %>%
  mutate(
    WT_mm_avg = rowMeans(cbind(WT_mm_rep2, WT_mm_rep3, WT_mm_rep4), na.rm = TRUE),
    KD_mm_avg = rowMeans(cbind(KD_mm_rep2, KD_mm_rep3, KD_mm_rep4), na.rm = TRUE))

# Join with KO Data and call new sites
borderline_joined <- inner_join(borderline_clean, KOdata, by = c("chr", "pos")) %>%
  mutate(
    diff_rep2 = WT_mm_rep2 - KO_mm_rate,
    diff_rep3 = WT_mm_rep3 - KO_mm_rate,
    diff_rep4 = WT_mm_rep4 - KO_mm_rate
  ) %>%
  select(chr, pos, motif,
  WT_mm_rep2, WT_mm_rep3, WT_mm_rep4, 
  KO_mm_rate, 
  diff_rep2, diff_rep3, diff_rep4,
  KD_mm_rep2, KD_mm_rep3, KD_mm_rep4, 
  WT_mm_avg, KD_mm_avg) %>%
  mutate(
    pool2 = paste(chr, pos, sep = "_") %in% pool2_keys,
    KO_called = {
      # Define the 4 logical conditions
      # We use coalesce() or is.na check to turn NA into FALSE for the sum
      cond1 <- coalesce(diff_rep2 >= 15, FALSE)
      cond2 <- coalesce(diff_rep3 >= 15, FALSE)
      cond3 <- coalesce(diff_rep4 >= 15,  FALSE)
      cond4 <- pool2
      
      # Sum the logicals (TRUE = 1, FALSE = 0) 
      # and check if the total is 2 or more
      (cond1 + cond2 + cond3 + cond4) >= 2
    }
  ) %>%
  arrange(desc(KO_called), desc(pool2))
write.csv(borderline_joined, file.path(active_dir, "borderline_KO_comparison.csv"), row.names = FALSE)

PUS7dep_joined <- inner_join(PUS7dep_clean, KOdata, by = c("chr", "pos")) %>%
  mutate(
    diff_rep2 = WT_mm_rep2 - KO_mm_rate,
    diff_rep3 = WT_mm_rep3 - KO_mm_rate,
    diff_rep4 = WT_mm_rep4 - KO_mm_rate
  ) %>%
  select(chr, pos, motif,
  WT_mm_rep2, WT_mm_rep3, WT_mm_rep4, 
  KO_mm_rate, 
  diff_rep2, diff_rep3, diff_rep4,
  KD_mm_rep2, KD_mm_rep3, KD_mm_rep4, 
  WT_mm_avg, KD_mm_avg) %>%
  mutate(
    pool2 = paste(chr, pos, sep = "_") %in% pool2_keys,
    KO_called = {
      pool2 = paste(chr, pos, sep = "_") %in% pool2_keys
      cond1 <- coalesce(diff_rep2 >= 15, FALSE)
      cond2 <- coalesce(diff_rep3 >= 15, FALSE)
      cond3 <- coalesce(diff_rep4 >= 15,  FALSE)
      (cond1 + cond2 + cond3 + pool2) >= 2
    }
  ) %>%
  arrange(desc(KO_called))
write.csv(PUS7dep_joined, file.path(active_dir, "UNUAR_KO_comparison.csv"), row.names = FALSE)


# ---- Plot everything with uncorrected values ----
generate_ko_scatter_plots <- function(df, plot_title, prefix, active_dir, kd_color = "#e25098", ko_color = "#ed93c0") {

  if (!dir.exists(active_dir)) {
    dir.create(active_dir, recursive = TRUE)
  }
  
  x_vars <- c("WT_mm_rep2", "WT_mm_rep3", "WT_mm_rep4", 
              "KD_mm_rep2", "KD_mm_rep3", "KD_mm_rep4", 
              "WT_mm_avg", "KD_mm_avg")
  
  # Helper function to generate clean x-axis labels
  get_clean_xlab <- function(var_name) {
    if (var_name == "WT_mm_avg") return("WT Average Mismatch %")
    if (var_name == "KD_mm_avg") return("KD Average Mismatch %")
    if (grepl("rep", var_name)) {
      parts <- strsplit(var_name, "_")[[1]]
      return(paste0(parts[1], " Mismatch % Rep", gsub("rep", "", parts[3])))
    }
    return(var_name)
  }
  
  # 3. Create the grouping columns for the 4 different versions
  plot_df <- df %>%
    mutate(
      UNUAR_match = grepl("^T[ACGT]TA[AG]$", motif),
      v1_group = ifelse(KO_called, "KO Called", "Other"),
      v2_group = ifelse(UNUAR_match, "UNUAR Motif", "Other"),
      v3_group = motif,
      v4_group = ifelse(KO_called, motif, "Other")
    )
  
  # Set up a consistent color palette for motifs to use in V3 and V4
  all_motifs <- unique(plot_df$motif)
  motif_colors <- hue_pal()(length(all_motifs))
  names(motif_colors) <- all_motifs
  # Add "Other" to the motif palette for V4
  v4_colors <- motif_colors
  v4_colors["Other"] <- "#797979"
  
  # 4. Loop through each x-variable to generate and save plots
  for (x_var in x_vars) {
    
    xlab <- get_clean_xlab(x_var)
    
    # Helper function to define the base plot aesthetics
    create_base_plot <- function(data_sorted) {
      ggplot(data_sorted, aes(x = .data[[x_var]], y = KO_mm_rate)) +
        scale_x_continuous(limits = c(0, 100)) +
        scale_y_continuous(limits = c(0, 100)) +
        labs(x = xlab, y = "KO Mismatch %", title = plot_title) +
        theme_scatter(base_size = 10) +
        theme(
          aspect.ratio = 1,
          legend.position = "right",
          legend.title = element_blank()
        )
    }
    # Helper function to handle saving
    save_plot <- function(p, version_suffix) {
      base_filename <- file.path(active_dir, paste0(prefix, "_", x_var, "_vs_KO_", version_suffix))
      ggsave(paste0(base_filename, ".png"), plot = p, width = 3.5, height = 3, units = "in")
      ggsave(paste0(base_filename, ".pdf"), plot = p, width = 2.5, height = 2, units = "in")
    }
    
    # --- Version 1: KO Called ---
    df_v1 <- plot_df %>% arrange(v1_group != "Other")
    p1 <- create_base_plot(df_v1) +
      geom_point(aes(color = v1_group), size = 1, alpha = 0.9) +
      scale_color_manual(values = c("Other" = "#797979", "KO Called" = "#daa250"))
    save_plot(p1, "v1_KO_called")
    
    # --- Version 2: UNUAR Motif ---
    df_v2 <- plot_df %>% arrange(v2_group != "Other")
    p2 <- create_base_plot(df_v2) +
      geom_point(aes(color = v2_group), size = 1, alpha = 0.9) +
      scale_color_manual(values = c("Other" = "#797979", "UNUAR Motif" = "#c279ca"))
    save_plot(p2, "v2_UNUAR_motif")
    
    # --- Version 3: All Motifs ---
    # Randomize plotting order to avoid one motif burying another
    set.seed(42)
    df_v3 <- plot_df %>% sample_frac(1)
    p3 <- create_base_plot(df_v3) +
      geom_point(aes(color = v3_group), size = 1, alpha = 0.9) +
      scale_color_manual(values = motif_colors)
    save_plot(p3, "v3_All_motifs")
    
    # --- Version 4: KO Called colored by Motif ---
    df_v4 <- plot_df %>% arrange(v4_group != "Other")
    p4 <- create_base_plot(df_v4) +
      geom_point(aes(color = v4_group), size = 1, alpha = 0.9) +
      scale_color_manual(values = v4_colors)
    save_plot(p4, "v4_KOCalled_by_motif")
  }
  # --- 2. The Connected Plot (WT_mm_avg vs KD_mm_avg and KO_mm_rate) ---
  df_connected_orig <- plot_df %>%
    mutate(chr_pos = paste(chr, pos, sep = "_")) %>%
    pivot_longer(
      cols = c(KD_mm_avg, KO_mm_rate),
      names_to = "Condition",
      values_to = "Mismatch_Rate"
    ) %>%
    mutate(Condition = ifelse(Condition == "KD_mm_avg", "KD", "KO"))
  
  p_conn_orig <- ggplot(df_connected_orig, aes(x = WT_mm_avg, y = Mismatch_Rate)) +
    geom_point(aes(color = Condition), size = 1, alpha = 0.9) +
    scale_x_continuous(limits = c(0, 100)) +
    scale_y_continuous(limits = c(0, 100)) +
    scale_color_manual(values = c("KD" = kd_color, "KO" = ko_color)) +
    labs(x = "WT Average Mismatch %", y = "PUS7 Depletion Mismatch %", title = plot_title) +
    theme_scatter(base_size = 10) +
    theme(
      aspect.ratio = 1,
      legend.position = "right",
      legend.title = element_blank()
    )
  
  # Final manual save for the connected plot
  base_filename_conn <- file.path(active_dir, paste0(prefix, "_WTavg_vs_KDKO"))
  ggsave(paste0(base_filename_conn, ".png"), plot = p_conn_orig, width = 3.5, height = 3, units = "in")
  ggsave(paste0(base_filename_conn, ".pdf"), plot = p_conn_orig, width = 2.5, height = 2, units = "in")
}

# PUS7-dependent sites
generate_ko_scatter_plots(
  df = PUS7dep_joined, 
  plot_title = "UNUAR Sites", 
  prefix = "UNUAR", 
  active_dir = paste0(active_dir, "/PUS7dep_plots")
)

# Borderline sites
generate_ko_scatter_plots(
  df = borderline_joined, 
  plot_title = "Borderline Sites", 
  prefix = "borderline", 
  active_dir = paste0(active_dir, "/borderline_plots")
)