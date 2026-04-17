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

# ---- Data ----
invitro_sum <- read_tsv("invitro_modification_significance.tsv")
invitro_PUS7 <- read_tsv("invitro_modified_sites_summary.tsv")

incell_PUS7 <- read_tsv("PUS7_dep_union_significant_summary.tsv")
incell_sum <- read_tsv("WT_mod_Both_significance.tsv")

PUS7_union <- full_join(incell_PUS7, invitro_PUS7, by = "chr")

# ---- Make Plots ----
# scatter plot
mpra_comp <- incell_sum %>%
    full_join(invitro_sum, by = "chr", suffix = c("_incell", "_invitro")) %>%
    select(chr, delta_delrate_incell, category_incell, 
              delta_delrate_invitro, category_invitro) %>%
    mutate(significance = case_when(
        category_incell == "Modified" & category_invitro == "Modified" 
            & chr %in% incell_PUS7$chr ~ "Intersection",
        category_incell == "Modified" & category_invitro == "Unmodified" 
            & chr %in% incell_PUS7$chr ~ "In Cellulo Only",
        category_incell == "Unmodified" & category_invitro == "Modified" 
            ~ "In Vitro Only",
        TRUE ~ "Inconclusive"
    )) %>%
    filter(chr %in% PUS7_union$chr)

p <- ggplot(mpra_comp, aes(x = delta_delrate_incell, y = delta_delrate_invitro, color = significance)) +
  geom_point() +
  labs(
    title = "MPRA PUS7 Comparison",
    x = "In Cell WT Delta Deletion Fraction",
    y = "In Vitro Delta Deletion Fraction"
  ) +
  coord_fixed(ratio = 1, xlim = c(0, 1), ylim = c(0, 1)) +
  theme_scatter(base_size = 12) +
  theme(legend.position = "right") +
  scale_color_manual(
    values = c("Intersection" = both, 
               "In Vitro Only" = invitro, 
               "In Cellulo Only" = incell),
    name = NULL, # Hides the legend title
    guide = "legend" # Ensures a legend is shown
  )
path <- file.path(plot_dir, "3_scatter_mpra_WT")
ggsave(paste0(path, ".pdf"), p, width = 4, height = 4, units = "in")
p_big <- p + theme_boxplot(base_size = 24)
ggsave(paste0(path, ".png"), p_big, width = 6, height = 6, units = "in")

# make the OE comparison plot
incell_oe_sum <- incell_all %>%
  group_by(chr) %>%
  summarize(
    avg_delrate_OE = mean(delrate[treat == "BS" & vector == "P3"], na.rm = TRUE),
    avg_delrate_input = mean(delrate[treat == "input"], na.rm = TRUE),
    delta_delrate_oe = avg_delrate_OE - avg_delrate_input,
  )
mpra_oe_comp <- incell_sum %>%
    full_join(invitro_sum, by = "chr", suffix = c("_incell", "_invitro")) %>%
    select(chr, delta_delrate_incell, category_incell, 
              delta_delrate_invitro, category_invitro) %>%
    mutate(significance = case_when(
        category_incell == "Modified" & category_invitro == "Modified" 
            & chr %in% incell_PUS7$chr ~ "Intersection",
        category_incell == "Modified" & category_invitro == "Unmodified" 
            & chr %in% incell_PUS7$chr ~ "In Cellulo Only",
        category_incell == "Unmodified" & category_invitro == "Modified" 
            ~ "In Vitro Only",
        TRUE ~ "Inconclusive"
    )) %>%
    filter(chr %in% PUS7_union$chr) %>%
    left_join(incell_oe_sum, by = "chr")
p <- ggplot(mpra_oe_comp, aes(x = delta_delrate_oe, y = delta_delrate_invitro, color = significance)) +
  geom_point() +
  labs(
    title = "MPRA PUS7 OE Comparison",
    x = "In Cell OE Delta Deletion Fraction",
    y = "In Vitro Delta Deletion Fraction"
  ) +
  coord_fixed(ratio = 1, xlim = c(0, 1), ylim = c(0, 1)) +
  theme_scatter(base_size = 12) +
  theme(legend.position = "right") +
  scale_color_manual(
    values = c("Intersection" = both, 
               "In Vitro Only" = invitro, 
               "In Cellulo Only" = incell),
    name = NULL, # Hides the legend title
    guide = "legend" # Ensures a legend is shown
  )
path <- file.path(plot_dir, "3_scatter_mpra_OE")
ggsave(paste0(path, ".pdf"), p, width = 4, height = 4, units = "in")
p_big <- p + theme_boxplot(base_size = 24)
ggsave(paste0(path, ".png"), p_big, width = 6, height = 6, units = "in")

# heat map
PUS7_heat_data <- incell_PUS7_avgs %>%
    inner_join(invitro_PUS7, by = "chr") %>%
    mutate(avg_delrate_IV = avg_delrate_BS, avg_delrate_input = avg_delrate_input.x) %>%
    select(chr, avg_delrate_input, avg_delrate_KD, avg_delrate_WT, avg_delrate_OE, avg_delrate_IV)
# determine site list
# sites in PUS7 union
# where KD < WT < OE, separeted by 0.1 --> 33 sites
PUS7_avgs_heat <- PUS7_heat_data %>%
  filter(chr %in% PUS7_union$chr) %>%
  filter(avg_delrate_KD < avg_delrate_WT & 
           avg_delrate_WT < avg_delrate_OE &
           (avg_delrate_WT - avg_delrate_KD) >= 0.1 &
           (avg_delrate_OE - avg_delrate_WT) >= 0.1)
PUS7_avgs_heat_long <- PUS7_avgs_heat %>%
  pivot_longer(
    cols = starts_with("avg_delrate"),
    names_to = "condition",
    values_to = "delrate"
  ) %>%
  mutate(condition = factor(condition, 
                            levels = c("avg_delrate_input", 
                                       "avg_delrate_KD", 
                                       "avg_delrate_WT", 
                                       "avg_delrate_OE",
                                       "avg_delrate_IV")))
# Create ordered factor for chr based on WT delrate values
ordered_chr <- PUS7_avgs_heat %>%
  arrange(avg_delrate_WT) %>%
  pull(chr)
p <- ggplot(PUS7_avgs_heat_long %>%
         mutate(chr = factor(chr, levels = ordered_chr)), 
       aes(x = condition, y = chr, fill = delrate)) +
  geom_tile() +
  scale_fill_gradient(low = "white",
                      high = OE,
                      na.value = "lightgray",
                      name = "Deletion Rate") +
  theme_heat(base_size = 12) +
  theme(axis.text.y = element_text(size = 7)) +
  labs(x = "", y = "Sites", title = "Deletion Rates by Site") +
  scale_x_discrete(labels = c(
    "avg_delrate_input" = "in",
    "avg_delrate_KD" = "KD",
    "avg_delrate_WT" = "WT",
    "avg_delrate_OE" = "OE",
    "avg_delrate_IV" = "IV"
  ))
path <- file.path(plot_dir, "3_heatmap_selectedsites")
ggsave(paste0(path, ".pdf"), p, width = 5, height = 4, units = "in")
p_big <- p + theme_boxplot(base_size = 24)
ggsave(paste0(path, ".png"), p_big, width = 6, height = 6, units = "in")

