#!/usr/bin/env Rscript

library(dplyr)
library(tidyr)
library(stringr)
library(readr)
library(readxl)
library(writexl)
library(ggplot2)
library(ggpubr)
library(ggrepel)

# ---- Directory ----
active_dir <- "~/PUS7regulation2026/Figure2/endo_sites"

# ---- Datasets ----
# --- endo sites ---
endo_data_1 <- read.csv(paste0(active_dir, "/endo1_counts.csv"))
endo_1 <- endo_data_1 %>%
    filter(pos == "287" | pos == "392", vector %in% c("pLKO", "P36")) %>%
    group_by(chr, pos, celltype, treat) %>%
    summarise(mean_delrate = mean(delrate, na.rm = TRUE), .groups = "drop") %>%
    pivot_wider(names_from = treat, values_from = mean_delrate) %>%
    mutate(delta_delrate = BS - input) %>%
    select(chr, celltype, delta_delrate) %>%
    pivot_wider(names_from = celltype, 
                values_from = delta_delrate,
                names_prefix = "endo_")

endo_data_2 <- read.csv(paste0(active_dir, "/endo2_counts.csv"))
endo_2 <- endo_data_2 %>%
  filter(pos == "145", vector %in% c("pLKO", "P36")) %>%
    group_by(chr, pos, celltype, treat) %>%
    summarise(mean_delrate = mean(delrate, na.rm = TRUE), .groups = "drop") %>%
    pivot_wider(names_from = treat, values_from = mean_delrate) %>%
    mutate(delta_delrate = BS - input) %>%
    select(chr, celltype, delta_delrate) %>%
    pivot_wider(names_from = celltype, 
                values_from = delta_delrate,
                names_prefix = "endo_")

endo_data_3 <- read.csv(paste0(active_dir, "/R/endoBID/20250606/endo3_counts.csv"))
endo_3 <- endo_data_3 %>%
  filter(pos == "479", vector %in% c("pLKO", "P36")) %>%
    group_by(chr, pos, celltype, treat) %>%
    summarise(mean_delrate = mean(delrate, na.rm = TRUE), .groups = "drop") %>%
    pivot_wider(names_from = treat, values_from = mean_delrate) %>%
    mutate(delta_delrate = BS - input) %>%
    select(chr, celltype, delta_delrate) %>%
    pivot_wider(names_from = celltype, 
                values_from = delta_delrate,
                names_prefix = "endo_")

endo_data_4 <- read_tsv(paste0(active_dir, "/endo4_counts.txt"))
endo_4 <- endo_data_4 %>%
    filter(pos == "261" & chr == "STIM1_chr11_4092507", vector %in% c("pLKO", "P36", "WT")) %>%
    group_by(chr, pos, celltype, treat) %>%
    summarise(mean_delrate = mean(delrate, na.rm = TRUE), .groups = "drop") %>%
    pivot_wider(names_from = treat, values_from = mean_delrate) %>%
    mutate(delta_delrate = BS - input) %>%
    select(chr, celltype, delta_delrate) %>%
    pivot_wider(names_from = celltype, 
                values_from = delta_delrate,
                names_prefix = "endo_")

endo_sites <- bind_rows(endo_1, endo_2, endo_3, endo_4)

comp_endo <- endo_sites %>%
    arrange(desc(endo_293T)) %>%
    rename(endo_HEK293T_WT = endo_293T, endo_HepG2_WT = endo_HepG2) %>%
    mutate(chr_pos = str_remove(id, "^[^_]+_"))

# ---- BID sites ----
bid_hela <- read_excel(paste0(active_dir, "/BID_Dai_2023.xlsx"), sheet = "Supplementary Table 5",
  skip = 2)
bid_293t <- read_excel(paste0(active_dir, "/BID_Dai_2023.xlsx"), sheet = "Supplementary Table 6",
  skip = 2)
bid_a549 <- read_excel(paste0(active_dir, "/BID_Dai_2023.xlsx"), sheet = "Supplementary Table 7",
  skip = 2)

bid_all <- bid_hela %>%
  select(chr, pos, BID_HeLa_WT = Deletion_Ave) %>%
  full_join(bid_293t %>% select(chr, pos, BID_HEK293T_WT = Deletion_Ave), by = c("chr", "pos")) %>%
  full_join(bid_a549 %>% select(chr, pos, BID_A549_WT = Deletion_Ave), by = c("chr", "pos")) %>% 
  mutate(chr_pos = paste0(chr, "_", pos))

comp_endo_bid <- comp_endo %>%
  left_join(bid_all, join_by(chr_pos == chr_pos))

# ---- PRAISE sites ----
praise <- read_excel(paste0(active_dir, "/PRAISE PUS-dependent sites Zhang2023.xlsx"), sheet = "Supplementary Dataset 2", skip = 1)

praise_sites <- praise %>%
  # separate(chr_site, into = c("chr", "pos", "extra"), 
  extract(chr_site, into = c("chr", "pos"), 
    regex = "^([^_]+).*?([0-9]+)",
    remove = FALSE, convert = TRUE) %>%
  rowwise() %>%
  mutate(PRAISE_HEK293T_WT = mean(c(rep1_deletion_ratio, rep2_deletion_ratio), na.rm = TRUE)) %>%
  ungroup() %>%
  select(chr, pos, PRAISE_HEK293T_WT) %>% 
  mutate(chr_pos = paste0(chr, "_", pos))

comp_endo_bid_praise <- comp_endo_bid %>%
  left_join(praise_sites, join_by(chr == chr, pos == pos))

# ---- Scatter Plot ----
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

p <- ggplot(comp_endo_bid_praise, 
            aes(x = endo_HEK293T_WT, y = PRAISE_HEK293T_WT)) +
  geom_point(size = 1) +
  stat_cor(method = "pearson", aes(label = ..r.label..), label.x = 0.05, label.y = 0.95, size = 3) +
  geom_text_repel(aes(label = site), size = 2, segment.size = 0.2) +
  labs(
    title = "HEK293T Deletion Fractions",
    x = "Endogenous WT Delta Deletion Fraction\n(BS - Input)",
    y = "PRAISE WT Deletion Fraction\n(BS)"
  ) +
  coord_fixed(ratio = 1, xlim = c(0, 1), ylim = c(0, 1)) +
  theme_scatter(base_size = 8)

path <- file.path(out_dir, "S2_scatter_293T_endo_praise_R")
ggsave(paste0(path, ".pdf"), p, width = 2.5, height = 2.5, units = "in")
ggsave(paste0(path, ".png"), p, width = 2.5, height = 2.5, units = "in")

p <- ggplot(comp_endo_bid_praise, 
            aes(x = endo_HEK293T_WT, y = BID_HEK293T_WT)) +
  geom_point(size = 1) +
  stat_cor(method = "pearson", aes(label = ..r.label..), label.x = 0.05, label.y = 0.95, size = 3) +
  geom_text_repel(aes(label = site), size = 2, segment.size = 0.2) +
  labs(
    title = "HEK293T Deletion Fractions",
    x = "Endogenous WT Delta Deletion Fraction\n(BS - Input)",
    y = "BID WT Deletion Fraction\n(BS)"
  ) +
  coord_fixed(ratio = 1, xlim = c(0, 1), ylim = c(0, 1)) +
  theme_scatter(base_size = 8)

path <- file.path(out_dir, "S2_scatter_293T_endo_bid_R")
ggsave(paste0(path, ".pdf"), p, width = 2.5, height = 2.5, units = "in")
ggsave(paste0(path, ".png"), p, width = 2.5, height = 2.5, units = "in")