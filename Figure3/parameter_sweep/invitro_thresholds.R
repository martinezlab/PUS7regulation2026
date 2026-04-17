#!/usr/bin/env Rscript

# ------------------------------------------------------------------------------
# Script: invitro_thresholds.R
# Description: Parameter sweep for Nano-BID-Amp modification analysis.
# ------------------------------------------------------------------------------

# 1. Check R Version
if (getRversion() < "4.3.0") {
  stop("ERROR: R version 4.3.0 or higher is required.")
}

# ---- Install and Load Libraries ----
required_packages <- c(
  "dplyr", "readr", "car", "tidyr", "purrr",
  "ggplot2", "ggtext", "blme", "future", "furrr", "ggseqlogo", "Biostrings"
)

missing_packages <- required_packages[!sapply(required_packages, requireNamespace, quietly = TRUE)]

if (length(missing_packages) > 0) {
  cat("Installing missing packages:", paste(missing_packages, collapse = ", "), "\n")
  if ("Biostrings" %in% missing_packages) {
    if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
    BiocManager::install("Biostrings")
    missing_packages <- missing_packages[missing_packages != "Biostrings"]
  }
  if (length(missing_packages) > 0) install.packages(missing_packages, repos = "http://cran.us.r-project.org")
}

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(car)
  library(tidyr)
  library(purrr)
  library(ggplot2)
  library(ggtext)
  library(blme)
  library(future)
  library(furrr)
  library(ggseqlogo)
  library(Biostrings)
})

# ==============================================================================
# ---- USER CONFIGURATION ----
# ==============================================================================

# Directories and Files
active_dir <- "~/Figure3/parameter_sweep"
input_file <- file.path(active_dir, "BIDdetect_data_invitro_delpos.txt")
output_dir <- file.path(active_dir, "invitro")
prefix     <- "invitro"

# Parallelization
num_cores <- 8

# Thresholds to Sweep
cov_thresholds   <- c(20, 50, 100)
sesoi_thresholds <- c(0.01, 0.02, 0.05, 0.1, 0.2)
pval_thresholds  <- c(0.001, 0.01, 0.05, 0.1)

# Baseline Definition (Must exist in the arrays above)
base_cov   <- 20
base_pval  <- 0.05
base_sesoi <- 0.05

# ==============================================================================

# ---- Directory Setup ----
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

log_file <- file.path(output_dir, paste0(prefix, "_analysis_sweep_log.txt"))
sink(log_file, split = TRUE) 

cat("--- Parameter Sweep Started ---\n")
cat("Output Directory:", output_dir, "\n\n")

# ---- Functions ----

create_sequence_logo <- function(kmers, title, subtitle = NULL) {
  valid_kmers <- na.omit(kmers)
  if (length(valid_kmers) == 0) return(NULL)
  
  motif_chars <- chartr("T", "U", as.character(valid_kmers))
  
  ppm <- Biostrings::consensusMatrix(motif_chars, as.prob = TRUE)
  
  for (nuc in c("A", "C", "G", "U")) {
    if (!(nuc %in% rownames(ppm))) ppm <- rbind(ppm, setNames(rep(0, ncol(ppm)), nuc))
  }
  ppm <- ppm[c("A", "C", "G", "U"), , drop = FALSE]
  
  p <- ggseqlogo(ppm, method = "prob") + 
    labs(title = title,
         subtitle = subtitle,
         y = "Probability") +
    scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, by = 0.2)) +
    theme_logo() +
    theme(
      axis.text.x = element_text(size = 9),
      axis.text.y = element_text(size = 9),
      axis.title.x = element_blank(),
      axis.title.y = element_text(size = 10),
      plot.title = element_text(size = 8, hjust = 0.5),
      plot.subtitle = element_text(size = 10, hjust = 0.5), 
      panel.background = element_blank(),
      panel.border = element_blank(),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      plot.margin = margin(10, 10, 10, 10, "pt")
    )

  safe_title <- gsub("\\.", "_", title)
  filename <- paste0(safe_title, "_seqlogo")
  path <- file.path(output_dir, filename)
  ggsave(paste0(path, ".pdf"), p, width = 2.75, height = 1.6, units = "in")
  ggsave(paste0(path, ".png"), p, width = 2.75, height = 1.6, units = "in")
  
  return(p)
}

# ---- Data Loading and Validation ----
if (!file.exists(input_file)) {
  stop(paste("Error: Input file not found at", input_file))
}
raw_data <- read_tsv(input_file, show_col_types = FALSE)

# Standardize treatment
raw_data <- raw_data %>%
  mutate(treat = case_when(
    vector == "noPUS" ~ "input",
    treat %in% c("in", "input") ~ "input",
    treat %in% c("BID", "BS") ~ "BS",
    TRUE ~ treat
  ))

# Set up parallel backend
plan(multisession, workers = num_cores)

sweep_summary_list <- list()
site_tracking_list <- list()

# ---- Outer Loop: Coverage ----
for (cov in cov_thresholds) {
  cat("\n=========================================\n")
  cat(sprintf("Processing Coverage Threshold: >= %d\n", cov))
  
  # Filter individual replicates based on coverage
  data_cov <- raw_data %>% dplyr::filter(totalReads >= cov)
  
  # Master summary for this coverage
  master_cov_summary <- data_cov %>%
    dplyr::group_by(chr, pos) %>%
    dplyr::summarise(
      avg_delrate_BS = mean(delrate[treat == "BS"], na.rm = TRUE),
      avg_delrate_input = mean(delrate[treat == "input"], na.rm = TRUE),
      psipos_5mer = na.omit(psipos_5mer)[1], 
      .groups = 'drop'
    ) %>%
    dplyr::mutate(delta_delrate = avg_delrate_BS - avg_delrate_input)
  
  # Prep data for modeling (need enough replicates)
  data_for_testing <- data_cov %>%
    dplyr::group_by(chr, pos) %>%
    dplyr::filter(n() >= 4, length(unique(treat)) >= 2, length(unique(rep)) >= 2) %>%
    dplyr::ungroup()
  
  cat(n_distinct(data_for_testing$chr, data_for_testing$pos), "sites ready for modeling at this coverage.\n")
  
  # Run modeling
  model_results <- data_for_testing %>%
    dplyr::group_by(chr, pos) %>%
    nest() %>%
    dplyr::mutate(stats = future_map(data, ~ {
      nested <- .x
      nested$treat <- factor(nested$treat, levels = c("input", "BS"))
      
      out <- tibble(
        p_raw = NA_real_,
        log_odds_ratio = NA_real_
      )
      
      glm_model <- tryCatch({
        bglmer(delrate ~ (1 | rep) + treat, data = nested,
               family = "binomial", weights = totalReads, fixef.prior = normal)
      }, error = function(e) NULL)
      
      if (!is.null(glm_model)) {
        # Significance (p-value)
        p_anova <- tryCatch({ Anova(glm_model, type = "III") }, error = function(e) NULL)
        if (!is.null(p_anova) && "treat" %in% rownames(p_anova)) {
          out$p_raw <- p_anova["treat", "Pr(>Chisq)"]
        }
        # Effect Size (Log-Odds Ratio)
        out$log_odds_ratio <- tryCatch({ fixef(glm_model)["treatBS"] }, error = function(e) NA_real_)
      }
      return(out)
    })) %>%
    dplyr::select(-data) %>%
    tidyr::unnest(stats) %>%
    ungroup() %>%
    mutate(p_adj = p.adjust(p_raw, method = "BH"))
  
  # Join stats back to master summary
  cov_final <- master_cov_summary %>%
    left_join(model_results, by = c("chr", "pos"))
  
  # ---- Inner Loop: SESOI and p-value ----
  for (pval in pval_thresholds) {
    for (sesoi in sesoi_thresholds) {
      
    # Determine Modified Sites
    mod_sites <- cov_final %>%
        dplyr::filter(!is.na(p_adj) & p_adj < pval & delta_delrate > sesoi) %>%
        dplyr::arrange(dplyr::desc(delta_delrate))
      
      # Build ID strings
      iter_id <- paste("Cov", cov, "pval", pval, "SESOI", sesoi, sep="_")
      iter_name <- gsub("\\.", "_", iter_id) # Safe name for files
      full_prefix <- paste0(prefix, "_", iter_name)
      
      # Store site lists for comparison
      site_keys <- paste(mod_sites$chr, mod_sites$pos, sep=":")
      site_tracking_list[[iter_id]] <- site_keys
      
      # Generate Sequence Logo
      if (nrow(mod_sites) > 0) {
        create_sequence_logo(
          kmers = mod_sites$psipos_5mer, 
          title = full_prefix, 
          subtitle = sprintf("n = %d modified sites", nrow(mod_sites))
        )
      }
      
      # Save Modified Sites TSV
      mod_out_path <- file.path(output_dir, paste0(full_prefix, "_modified_sites.tsv"))
      write_tsv(mod_sites, mod_out_path)
      
    }
  }
}

plan(sequential)

# ---- Baseline Comparison & Summary Generation ----
cat("\n--- Running Baseline Comparisons ---\n")

baseline_id <- paste("Cov", base_cov, "pval", base_pval, "SESOI", base_sesoi, sep="_")
if (!(baseline_id %in% names(site_tracking_list))) {
  warning("Baseline configuration was not generated in the sweep. Summary comparison will fail.")
} else {
  baseline_sites <- site_tracking_list[[baseline_id]]
  
  # Compute Gained/Lost metrics for all iterations
  summary_rows <- lapply(names(site_tracking_list), function(iter) {
    current_sites <- site_tracking_list[[iter]]
    
    # Parse parameters from ID
    params <- strsplit(iter, "_")[[1]]
    
    tibble(
      Iteration_ID = iter,
      Coverage_Min = as.numeric(params[2]),
      P_Value_Max = as.numeric(params[4]),
      SESOI_Min = as.numeric(params[6]),
      Total_Modified = length(current_sites),
      Retained_vs_Baseline = length(intersect(current_sites, baseline_sites)),
      Gained_vs_Baseline = length(setdiff(current_sites, baseline_sites)),
      Lost_vs_Baseline = length(setdiff(baseline_sites, current_sites))
    )
  }) %>% bind_rows()
  
  # Save Summary Table
  summary_path <- file.path(output_dir, paste0(prefix, "_parameter_sweep_summary.tsv"))
  write_tsv(summary_rows, summary_path)
  
  cat("Parameter Sweep complete. Comparison summary saved to:\n", summary_path, "\n")
  print(summary_rows)
}

sink() # Stop logging