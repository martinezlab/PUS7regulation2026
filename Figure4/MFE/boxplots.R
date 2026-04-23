#!/usr/bin/env Rscript

library(dplyr)
library(ggplot2)

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
theme_boxplot <- function(base_size = 16) {
  theme_base_custom(base_size) +
    theme(legend.position = "none", panel.grid.major.x = element_blank(), panel.grid.minor = element_blank())
}

both <- "#c279ca"
plot_dir <- "~/PUS7regulation2026/Figure4/plots"

# ---- Datasets ----
rnafold_deltaG_pool1 <- read.csv("pool1_noadapters_rnafold_deltag.csv")
ivSHAPE_deltaG_pool1 <- read.csv("IV_NoPus_rnafold_summary.csv")

invitro_full <- read.csv("~/PUS7regulation2026/Figure4/sequence/invitroBID_MPRA_summary.csv", header = TRUE)
invitro_PUS7 <- invitro_full %>% filter(PUS7_dependent == "Yes") %>% filter(!(is.na(psipos_motif)))
incell_full <- read.csv("~/PUS7regulation2026/Figure4/sequence/incellBID_WTsummary.csv", header = TRUE)
incell_PUS7 <- incell_full %>% filter(PUS7_dependent == "Yes") %>% filter(!(is.na(psipos_motif)))
PUS7_intersection <- inner_join(incell_PUS7, invitro_PUS7, by = "site")

# --- RNAfold Only ---
RNAfold_MFE <- rnafold_deltaG_pool1 %>%
  mutate(labels = case_when(
    sequence_id %in% PUS7_intersection$site ~ "PUS7-dependent",
    !sequence_id %in% PUS7_intersection$site ~ "Not Dependent"
  )) %>%
  mutate(labels = factor(labels, levels = c("Not Dependent", "PUS7-dependent")))

p <- ggplot(RNAfold_MFE, aes(x = labels, y = mfe_energy)) +
  geom_boxplot(fill = both) +
  labs(
    x = "",
    y = "MFE Structure Energy",
    title = "Intersection"
  ) +
  theme_boxplot(base_size = 12)
path <- file.path(plot_dir, "boxplot_rnafoldMFE_intersection")
ggsave(paste0(path, ".pdf"), p, width = 2, height = 2.25, units = "in")
p_big <- p + theme_boxplot(base_size = 24)
ggsave(paste0(path, ".png"), p_big, width = 6, height = 6, units = "in")


# --- in vitro SHAPE informed ---
ivSHAPE_MFE <- ivSHAPE_deltaG_pool1 %>%
  mutate(labels = case_when(
    sequence_id %in% PUS7_intersection$site ~ "PUS7-dependent",
    !sequence_id %in% PUS7_intersection$site ~ "Not Dependent"
  )) %>%
  mutate(labels = factor(labels, levels = c("Not Dependent", "PUS7-dependent")))

p <- ggplot(ivSHAPE_MFE, aes(x = labels, y = mfe_energy)) +
  geom_boxplot(fill = both) +
  labs(
    x = "",
    y = "MFE Structure Energy",
    title = "Intersection"
  ) +
  theme_boxplot(base_size = 12)
path <- file.path(plot_dir, "boxplot_MFE_intersection")
ggsave(paste0(path, ".pdf"), p, width = 1.8, height = 2.25, units = "in")
p_big <- p + theme_boxplot(base_size = 24)
ggsave(paste0(path, ".png"), p_big, width = 6, height = 6, units = "in")
