#!/usr/bin/env Rscript

if (getRversion() < "4.3.0") {
  stop("ERROR: R version 4.3.0 or higher is required.")
}
required_packages <- c("dplyr", "readr", "car", "tidyr", "purrr", "ggplot2", "blme", "forcats")
missing_packages <- required_packages[!sapply(required_packages, requireNamespace, quietly = TRUE)]
if (length(missing_packages) > 0) {
  cat("Installing missing packages:", paste(missing_packages, collapse = ", "), "\n")
  install.packages(missing_packages, repos = "http://cran.us.r-project.org")
}
suppressPackageStartupMessages({lapply(required_packages, library, character.only = TRUE)})


# ---- Plot Themes ----
invitro <- "#a2a2fc"
pupup <- "#fce87e"
unuar <- "#e25098"
distal <- "#d3d3d3"

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

# ---- Data ----
fis1_data <- read.csv("~/PUS7regulation2026/Figure5/ASO/ASO_FIS1all_data.csv")

# ---- FIS1 Plot ----

# plot just the FIS1 site
plot2_data <- fis1_data %>% filter(treat == "BS") %>%
  mutate(chr_formatted = gsub("^([^_]+)_([^_]+)_(.*)$", "\\1\n\\2:\\3", chr)) %>% 
  filter(chr == "FIS1_chr7_101245032")

p <- ggplot(data = plot2_data, mapping = aes(x = method, fill = method)) +
  stat_summary(
    aes(y = delrate), 
    fun = \(x) mean(x, na.rm = TRUE),
    geom = "bar", 
    width = 0.8
  ) +
  geom_jitter(aes(y = delrate), width = 0.2, alpha = 0.5) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.1))) +
  scale_fill_manual(values = c(
    "No ASO" = invitro, 
    "Pupup" = pupup,
    "Control" = distal, 
    "UNUAR" = unuar
  )) +
  labs(
    x = "",
    y = "Deletion Fraction"
  ) +
  theme_bar(base_size = 11) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "none",
    strip.text = element_text(size = rel(0.8))
  ) +
  facet_wrap(~ chr_formatted, ncol = 5)
plot_path <- file.path(active_dir, "barplot_ASO_FIS1only")
ggsave(paste0(plot_path, ".pdf"), p, width = 1.5, height = 2, units = "in")
ggsave(paste0(plot_path, ".png"), p, width = 1.5, height = 2, units = "in", dpi = 600)


# plot all sites but FIS1
plot_data <- fis1_data %>% filter(treat == "BS") %>%
  mutate(chr_formatted = gsub("^([^_]+)_([^_]+)_(.*)$", "\\1\n\\2:\\3", chr)) %>% 
  filter(!(chr == "FIS1_chr7_101245032"))

p <- ggplot(data = plot_data, mapping = aes(x = method, fill = method)) +
  stat_summary(
    aes(y = delrate), 
    fun = \(x) mean(x, na.rm = TRUE),
    geom = "bar", 
    width = 0.8
  ) +
  geom_jitter(aes(y = delrate), width = 0.2, alpha = 0.5) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.1))) +
  scale_fill_manual(values = c(
    "No ASO" = invitro, 
    "Pupup" = pupup,
    "Control" = distal, 
    "UNUAR" = unuar
  )) +
  labs(
    x = "",
    y = "Deletion Fraction"
  ) +
  theme_bar(base_size = 11) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "none",
    strip.text = element_text(size = rel(0.8))
  ) +
  facet_wrap(~ chr_formatted, ncol = 5)
plot_path <- file.path(active_dir, "barplot_ASO_FIS1others")
ggsave(paste0(plot_path, ".pdf"), p, width = 4.7, height = 2, units = "in")
ggsave(paste0(plot_path, ".png"), p, width = 4.7, height = 2, units = "in", dpi = 600)

# ---- Statistics ----

active_dir <- "~/PUS7regulation2026/Figure5/mutagenesis_analysis"
output_dir      <- file.path(active_dir, "FIS_stats")
raw_dir <- file.path(output_dir, "data_raw")
sum_dir <- file.path(output_dir, "data_summary")

dir.create(raw_dir, showWarnings=FALSE, recursive=TRUE)
dir.create(sum_dir, showWarnings=FALSE, recursive=TRUE)

# Execution Options
sesoi           <- 0.05    # Single threshold for TOST and effect size filtering

# Log File
log_file <- file.path(output_dir, "stats_log.txt")
sink(log_file, split = TRUE)

# ---- Functions ----
# Function: Factor Dependency  Analysis
perform_anova_test <- function(data, factor_column, random_effects_formula = "") {
  # Compares two models to test for a significant interaction between 'treat' and another factor.
  
  formula_factor_str <- paste("delrate ~ (1|rep)", random_effects_formula, "+", factor_column, "+ L1andBS.ind + L2andBS.ind")
  formula_nofactor_str <- paste("delrate ~ (1|rep)", random_effects_formula, "+", factor_column, "+ BS.ind")
  
  model_factor_formula <- as.formula(formula_factor_str)
  model_nofactor_formula <- as.formula(formula_nofactor_str)

  results <- data %>%
    group_by(chr, pos) %>%
    nest() %>%
    mutate(result = map(data, ~ {
      nested_data <- .x

      if (nrow(nested_data) < 2 || length(unique(nested_data$treat)) < 2 || 
          length(unique(nested_data$rep)) < 2 || length(unique(nested_data[[factor_column]])) < 2) {
        return(tibble(term = NA, p.value = NA, note = "Not enough data"))
      }

      glm_model_factor <- tryCatch(bglmer(model_factor_formula, data = nested_data, 
                                         family="binomial", weights=totalReads, fixef.prior=normal), 
                                   error = function(e) NULL)
      glm_model_nofactor <- tryCatch(bglmer(model_nofactor_formula, data = nested_data, 
                                           family="binomial", weights=totalReads, fixef.prior=normal), 
                                     error = function(e) NULL)
      
      if (is.null(glm_model_factor) || is.null(glm_model_nofactor)) {
        return(tibble(term = NA, p.value = NA, note = "Model fitting failed"))
      }
      
      p_anova <- tryCatch(anova(glm_model_factor, glm_model_nofactor, test = "ChiSq"), error = function(e) NULL)

      if (!is.null(p_anova)) {
        anova_results <- as_tibble(p_anova, rownames = "term")
        return(select(anova_results, term, p.value = `Pr(>Chisq)`) %>% mutate(note = ""))
      } else {
        return(tibble(term = NA, p.value = NA, note = "ANOVA failed"))
      }
    }, .options = furrr_options(seed = TRUE))) %>%
    unnest(result) %>%
    filter(term == "glm_model_factor") # Keep only the row with the p-value
  
  results %>% mutate(p.value.BH = p.adjust(p.value, method = "BH"))
}

# Function: Run the main factor dependency analysis (above)
run_factor_dependency_analysis <- function(data, factor_column, levels, random_effects_formula, sesoi, direction = "positive", return_all_sites = TRUE) {
  # Encapsulates the entire workflow for a single factor dependency test.
  # Returns a data frame of significant sites that also pass the effect size filter.
  
  baseline_level <- levels[1]
  experimental_level <- levels[2]
  
  prepared_data <- data %>%
    mutate(!!sym(factor_column) := factor(!!sym(factor_column), levels = levels)) %>%
    mutate(BS.ind = ifelse(treat == "BS", 1, 0),
           L1andBS.ind = ifelse(treat == "BS" & !!sym(factor_column) == baseline_level, 1, 0),
           L2andBS.ind = ifelse(treat == "BS" & !!sym(factor_column) == experimental_level, 1, 0))

  sig_results <- perform_anova_test(prepared_data, factor_column, random_effects_formula)

  delta_delrates <- prepared_data %>%
    group_by(chr, pos, .data[[factor_column]], treat) %>%
    summarise(mean_delrate = mean(delrate, na.rm = TRUE), .groups = "drop") %>%
    pivot_wider(names_from = treat, values_from = mean_delrate) %>%
    mutate(delta_delrate = BS - input) %>%
    select(chr, pos, .data[[factor_column]], delta_delrate) %>%
    pivot_wider(names_from = all_of(factor_column), 
                values_from = delta_delrate,
                names_prefix = "delta_delrate_") %>%
    mutate(dd_delrate = .data[[paste0("delta_delrate_", experimental_level)]] - .data[[paste0("delta_delrate_", baseline_level)]])

  # Join all statistical results with all effect sizes
  all_sites_results <- sig_results %>%
    left_join(delta_delrates, by = "chr")
    
  # Conditionally return ALL sites or filter them
  if (return_all_sites) {
    # If TRUE, return everything. No filtering is applied.
    return(all_sites_results)
  } else {
    # If FALSE (the default), perform the original filtering on BOTH p-value and effect size.
    final_sites <- all_sites_results %>%
      filter(!is.na(p.value.BH), p.value.BH < 0.05) %>%
      filter(
        case_when(
          !is.na(dd_delrate) & direction == "positive" ~ dd_delrate > sesoi,
          !is.na(dd_delrate) & direction == "both"     ~ abs(dd_delrate) > sesoi,
          TRUE                                         ~ FALSE
        )
      )
    return(final_sites)
  
  return(final_sites)
}}

# Function: Save summary and raw data for a list of sites
save_site_subset <- function(sites_df, input_data, file_prefix) {
  # Standardizes the saving of a summary table and its corresponding raw data.
  
  if (nrow(sites_df) > 0) {
    # Save the summary table (the list of significant sites)
    summary_path <- file.path(sum_dir, paste0(file_prefix, "_summary.tsv"))
    write_tsv(sites_df, summary_path, progress = FALSE)
    
    # Filter the main raw data to get only the rows for these sites
    raw_data_subset <- input_data %>%
      semi_join(sites_df, by = c("chr"))
    
    # Save the raw data subset
    raw_path <- file.path(raw_dir, paste0(file_prefix, "_raw_data.tsv"))
    write_tsv(raw_data_subset, raw_path, progress = FALSE)
    
    cat("   Saved", nrow(sites_df), "sites to summary and raw data files with prefix:", file_prefix, "\n")
  } else {
    cat("   Found 0 sites to save for prefix:", file_prefix, "\n")
  }
}

# Function: Save a list of sites and print summary statistics for a given column
summarize_and_save_sites <- function(sites_df, input_data, file_prefix, summary_col, summary_title) {
  # 1. Save the summary and raw data files using the existing helper function
  save_site_subset(sites_df, input_data, file_prefix)
  # 2. Check if there are any sites to summarize
  if (nrow(sites_df) > 0 && summary_col %in% colnames(sites_df)) {
    cat("\n--- ", summary_title, " ---\n")
    # 3. Calculate summary statistics for the specified column.
    summary_stats <- sites_df %>%
      ungroup() %>%
      summarise(
        n_sites = n(),
        mean_value = mean(.data[[summary_col]], na.rm = TRUE),
        sd_value = sd(.data[[summary_col]], na.rm = TRUE),
        min_value = min(.data[[summary_col]], na.rm = TRUE),
        max_value = max(.data[[summary_col]], na.rm = TRUE),
        .groups = 'drop'
      )
    print(summary_stats, width = Inf)
  }
}

# ---- Analysis ----

ASO_Pupup <- fis1_data %>%
  filter(
    method %in% c("No ASO", "Pupup")
  )

final_sites <- run_factor_dependency_analysis(
    data = ASO_Pupup,
    factor_column = "method",
    levels = c("No ASO", "Pupup"),
    random_effects_formula = "",
    sesoi = sesoi,
    direction = "both"
    )

summarize_and_save_sites(
    sites_df = final_sites,
    input_data = ASO_Pupup,
    file_prefix = "ASO_Pupup",
    summary_col = "dd_delrate",
    summary_title = "Summary for Pupup ASO"
    )

ASO_UNUAR <- fis1_data %>%
  filter(
    method %in% c("No ASO", "UNUAR")
  )

final_sites <- run_factor_dependency_analysis(
    data = ASO_UNUAR,
    factor_column = "method",
    levels = c("No ASO", "UNUAR"),
    random_effects_formula = "",
    sesoi = sesoi,
    direction = "both"
    )

summarize_and_save_sites(
    sites_df = final_sites,
    input_data = ASO_UNUAR,
    file_prefix = "ASO_UNUAR",
    summary_col = "dd_delrate",
    summary_title = "Summary for UNUAR ASO"
    )

ASO_Control <- fis1_data %>%
  filter(
    method %in% c("No ASO", "Control")
  )

final_sites <- run_factor_dependency_analysis(
    data = ASO_Control,
    factor_column = "method",
    levels = c("No ASO", "Control"),
    random_effects_formula = "",
    sesoi = sesoi,
    direction = "both"
    )

summarize_and_save_sites(
    sites_df = final_sites,
    input_data = ASO_Control,
    file_prefix = "ASO_Control",
    summary_col = "dd_delrate",
    summary_title = "Summary for Control ASO"
    )

sink()