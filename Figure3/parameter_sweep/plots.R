#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(stringr)
  library(ggplot2)
  library(tidyr)
})

# ---- Directory ----
active_dir <- "~/Figure3/parameter_sweep"

input_raw_file <- file.path(active_dir, "BIDdetect_data_invitro_delpos.txt")
sweep_dir <- file.path(active_dir, "invitro")
prefix <- "invitro"

# UNUAR motif definition as TNTAR
target_regex <- "^T.TA[AG]$"

# ---- 1. Calculate Background Enrichment ----
raw_data <- read_tsv(input_raw_file, show_col_types = FALSE)

unique_sites <- raw_data %>%
  dplyr::select(chr, pos, psipos_5mer) %>%
  dplyr::distinct() %>%
  tidyr::drop_na(psipos_5mer)

total_background_sites <- nrow(unique_sites)
background_unuar_count <- sum(stringr::str_detect(unique_sites$psipos_5mer, target_regex))
background_freq <- background_unuar_count / total_background_sites

cat(sprintf("Background UNUAR frequency: %.2f%% (%d / %d sites)\n", 
            background_freq * 100, background_unuar_count, total_background_sites))

# ---- 2. Parse Sweep Iterations ----
# Find all modified_sites.tsv files from the sweep
# These are numerous and not uploaded, but should be created by invitro_thresholds.R
mod_files <- list.files(sweep_dir, pattern = "_modified_sites\\.tsv$", full.names = TRUE)

results_list <- list()

for (file in mod_files) {
  # Extract iteration parameters from the filename
  filename <- basename(file)
  iter_id <- str_remove(filename, paste0(prefix, "_"))
  iter_id <- str_remove(iter_id, "_modified_sites\\.tsv$")
  
  # Read the modified sites
  mod_data <- read_tsv(file, show_col_types = FALSE) %>%
    tidyr::drop_na(psipos_5mer)
  
  total_mod <- nrow(mod_data)
  
  if (total_mod == 0) {
    next # Skip if no modified sites
  }
  
  kmers <- mod_data$psipos_5mer
  
  # Strict Conformity & Enrichment Calculation
  strict_match_count <- sum(stringr::str_detect(kmers, target_regex))
  
  # The conformity score is now strictly the frequency of perfect matches
  conformity_score <- strict_match_count / total_mod 
  fold_enrichment <- conformity_score / background_freq
  
  # Parse individual parameters for coloring/faceting
  cov_val <- as.numeric(str_replace(str_extract(iter_id, "Cov_?\\d+"), "Cov_?", ""))
  
  pval_raw <- str_replace(str_extract(iter_id, "pval_?\\d+_\\d+"), "pval_?", "")
  pval_val <- as.numeric(str_replace(pval_raw, "_", "."))
  
  sesoi_raw <- str_replace(str_extract(iter_id, "SESOI_?\\d+_\\d+"), "SESOI_?", "")
  sesoi_val <- as.numeric(str_replace(sesoi_raw, "_", "."))
  
  results_list[[iter_id]] <- tibble(
    Iteration = iter_id,
    Coverage = as.factor(cov_val),
    P_Value = as.factor(pval_val),
    SESOI = as.factor(sesoi_val),
    Total_Sites = total_mod,
    Conformity_Score = conformity_score,
    Fold_Enrichment = fold_enrichment
  )
}

plot_data <- bind_rows(results_list)

# ---- 3. Plotting ----
theme_custom <- function(base_size = 14) {
  theme_minimal(base_size = base_size, base_family = "sans") %+replace%
    theme(
      panel.border = element_rect(color = "black", fill = NA, size = 1),
      plot.title = element_text(hjust = 0.5, size = rel(1.2), margin = margin(b = 5)),
      plot.subtitle = element_text(hjust = 0.5, size = rel(0.8), margin = margin(b = 5)),
      axis.title = element_text(size = rel(1.0)),
      axis.text = element_text(size = rel(0.8)),
      plot.margin = margin(5, 5, 5, 5),
      legend.position = "right"
    )
}

plot_data$Coverage <- factor(plot_data$Coverage, levels = c(20, 50, 100))

p_pareto <- ggplot(plot_data, aes(x = Total_Sites, y = Conformity_Score, color = SESOI, shape = P_Value)) +
  geom_point(size = 3, alpha = 0.8) +
  scale_color_viridis_d(option = "turbo", direction = -1) +
  facet_wrap(~ Coverage, labeller = label_both) + # Splits plot into panels by Coverage
  labs(
    title = "In Vitro",
    # subtitle = paste("Background Strict Match Frequency:", round(background_freq * 100, 2), "%"),
    x = "Total Modified Sites",
    y = "Fraction UNUAR",
    shape = "Max p-value",
    color = "Min Delta\nDeletion\nFraction"
  ) +
  theme_custom(base_size = 12) +
  theme(strip.text = element_text(size = 12))

plot_path_pareto <- file.path(active_dir, paste0("S3_invitro_ParetoFront"))
ggsave(paste0(plot_path_pareto, ".pdf"), p_pareto, width = 7.5, height = 2.75, units = "in")
ggsave(paste0(plot_path_pareto, ".png"), p_pareto, width = 7.5, height = 2.75, units = "in", dpi = 300)

# ---- Metaplots ----

ivSHAPE_pairprob_pool1 <- read.csv(file.path(active_dir, "R/Pool1_SHAPE/IV_NoPus_pairprob_aligned.csv"))
pool1 <- read.csv(file.path(active_dir, "R/Pool1_Sequences/pool1_cleaned_20251107.csv"))

ivSHAPE_pairprob_pool1 <- ivSHAPE_pairprob_pool1 %>%
  dplyr::rename(chr = rna_name)

valid_chrs <- pool1 %>%
  tidyr::drop_na(psipos_5mer) %>%
  dplyr::pull(clean_name) %>%
  unique()

# 3. Define the specific iterations to extract
target_conditions <- tibble::tribble(
  ~cov, ~pval, ~sesoi,
  20,   0.1,   0.01,
  20,   0.05,  0.05,
  50,   0.05,  0.01,
  50,   0.05,  0.05,
  100,  0.05,  0.05,
  100,  0.001, 0.2
)

# Extract all unique levels across the whole sweep so the color and shape mapping 
# remains completely identical to your Pareto front, even when plotting one line.
all_sesoi_levels <- as.character(sort(unique(target_conditions$sesoi)))
all_pval_levels <- as.character(sort(unique(target_conditions$pval)))

plot_data_list <- list()
facet_levels <- c()

# 4. Extract, Process, Plot, and Save Individually
for (i in 1:nrow(target_conditions)) {
  c_cov <- target_conditions$cov[i]
  c_pval <- target_conditions$pval[i]
  c_sesoi <- target_conditions$sesoi[i]
  
  # Format decimals to match the standard filename structure (e.g., 0.05 -> 0_05)
  pval_str <- stringr::str_replace(as.character(c_pval), "\\.", "_")
  sesoi_str <- stringr::str_replace(as.character(c_sesoi), "\\.", "_")
  
  facet_label <- sprintf("Coverage > %d\np < %s\nDelFrac > %s", c_cov, c_pval, c_sesoi)
  facet_levels <- c(facet_levels, facet_label)

  # Reconstruct the expected base filename
  base_name <- sprintf("%s_Cov_%d_pval_%s_SESOI_%s", prefix, c_cov, pval_str, sesoi_str)
  input_file <- file.path(sweep_dir, paste0(base_name, "_modified_sites.tsv"))
  
  if (file.exists(input_file)) {
    mod_data <- read_tsv(input_file, show_col_types = FALSE)
    
    # Filter for valid chromosomes, join, and summarize the single condition
    plot_df <- mod_data %>%
      dplyr::filter(chr %in% valid_chrs) %>%
      dplyr::inner_join(ivSHAPE_pairprob_pool1, by = "chr") %>%
      dplyr::summarise(
        dplyr::across(
          dplyr::starts_with("pos_"), 
          list(
            mean = ~mean(.x, na.rm = TRUE), 
            se = ~sd(.x, na.rm = TRUE) / sqrt(sum(!is.na(.x)))
          )
        )
      ) %>%
      tidyr::pivot_longer(
        cols = dplyr::starts_with("pos_"),
        names_to = c("Position", ".value"),
        names_pattern = "pos_(.*)_(mean|se)"
      ) %>%
      dplyr::mutate(
        Position = as.numeric(sub("m", "-", Position)),
        # Lock in the factor levels so the single line grabs the correct color/shape
        SESOI = factor(as.character(c_sesoi), levels = all_sesoi_levels),
        P_Value = factor(as.character(c_pval), levels = all_pval_levels),
        Facet_Label = facet_label
      ) %>%
      dplyr::arrange(Position)
    
    plot_data_list[[i]] <- plot_df

    # Build the individual plot
    p <- ggplot(plot_df, aes(x = Position, y = mean, color = SESOI, fill = SESOI)) +
      geom_ribbon(aes(ymin = mean - se, ymax = mean + se), alpha = 0.2, linetype = 0) +
      geom_line(aes(linetype = P_Value), linewidth = 1) +
      geom_point(aes(shape = P_Value), size = 0.5, data = . %>% filter(Position %% 5 == 0)) + 
      # drop = FALSE prevents ggplot from discarding unused factor levels in the legend/palette
      scale_color_viridis_d(option = "turbo", direction = -1, drop = FALSE) +
      scale_fill_viridis_d(option = "turbo", direction = -1, drop = FALSE) +
      scale_shape_discrete(drop = FALSE) +
      # scale_linetype_discrete(drop = FALSE) +
      labs(
        # title = "In Vitro Folding Probability",
        title = sprintf("Coverage: %d | Max p-value: %s | SESOI: %s", c_cov, c_pval, c_sesoi),
        x = "Position",
        y = "Average Pairing Probability",
        color = "Min Delta\nDeletion\nFraction (SESOI)",
        fill  = "Min Delta\nDeletion\nFraction (SESOI)",
        shape = "Max p-value",
        linetype = "Max p-value"
      ) +
      coord_cartesian(ylim = c(0, 1)) +
      scale_y_continuous(breaks = seq(0, 1, 0.2), expand = c(0, 0)) +
      theme_custom(base_size = 8) +
      theme(
        # legend.key.width = unit(1.5, "cm"),
        legend.position = "none",
        panel.grid.minor = element_blank()
      )
    
    # Save directly to sweep_dir with matching nomenclature
    out_pdf <- file.path(sweep_dir, paste0(base_name, "_metaplot.pdf"))
    out_png <- file.path(sweep_dir, paste0(base_name, "_metaplot.png"))
    
    ggsave(out_pdf, plot = p, width = 2, height = 2)
    ggsave(out_png, plot = p, width = 2, height = 2)
    
  } else {
    warning(sprintf("File not found: %s", input_file))
  }
}

# 5. Combine Data and Make Faceted Plot
combined_plot_df <- dplyr::bind_rows(plot_data_list) %>%
  dplyr::mutate(Facet_Label = factor(Facet_Label, levels = facet_levels))

p_faceted <- ggplot(combined_plot_df, aes(x = Position, y = mean)) +
  geom_ribbon(aes(ymin = mean - se, ymax = mean + se), alpha = 0.2, linetype = 0) +
  geom_line(linewidth = 1) + #aes(linetype = P_Value), 
  geom_point(size = 0.5, data = . %>% filter(Position %% 5 == 0)) + #aes(shape = P_Value), 
  facet_wrap(~ Facet_Label, nrow = 1) + 
  labs(
    # title = "In Vitro Folding Probability",
    # subtitle = "Position -20 to +20 across select sweep conditions",
    x = "Position",
    y = "Average Pairing Probability"
  ) +
  scale_x_continuous(breaks = seq(-20, 20, by = 10), limits = c(-20, 20), expand = c(0, 0)) +
  coord_cartesian(ylim = c(0, 1)) +
  scale_y_continuous(breaks = seq(0, 1, 0.2), expand = c(0, 0)) +
  theme_custom(base_size = 8) +
  theme(
    # legend.key.width = unit(1.5, "cm"),
    legend.position = "none",
    panel.grid.minor = element_blank(),
    strip.text = element_text(size = 8), # Slightly smaller facet text to fit
    panel.spacing = unit(0.5, "lines") # Give plots breathing room
  )

# 6. Save the Faceted Plot
out_pdf_faceted <- file.path(active_dir, "S3_ivSHAPE_Metaplots_Faceted.pdf")
out_png_faceted <- file.path(active_dir, "S3_ivSHAPE_Metaplots_Faceted.png")

# Width is set wide to accommodate 6 plots in a single row nicely
ggsave(out_pdf_faceted, plot = p_faceted, width = 7.5, height = 1.75)
ggsave(out_png_faceted, plot = p_faceted, width = 7.5, height = 1.75)
