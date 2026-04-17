#!/usr/bin/env Rscript

library(dplyr)
library(tidyverse)
library(ggplot2)
library(scales)
library(RColorBrewer)
library(broom)

# ---- Directory ----
active_dir <- "~/PUS7regulation2026/Figure1/UNUARstandard"

# ---- Data Loading & Calculations ----

data <- read.table(paste0(active_dir,"/UNUAR_counts.txt", header = TRUE))

data <- data %>%
  mutate(
    mm_UtoC = case_when(
      ref != "T" ~ 0,
      (T.count + C.count) == 0 ~ 0,
      TRUE ~ C.count / (T.count + C.count)
  )
)

data <- data %>%
  tidyr::extract(
    sample,
    into = "UNpsiAR_conc",
    regex = "_([0-9.]+)",
    remove = FALSE,
    convert = TRUE
  )

data <- data %>%
  group_by(pos) %>%
  mutate(
    mm_UtoC_perc = mm_UtoC * 100
  ) %>%
  ungroup()

# ---- Plot the U-to-C Mismatch Rate ----

all_samples <- unique(data$sample)
other_samples <- all_samples[all_samples != "UNpsiAR_0"]

other_colors <- brewer.pal(n = length(other_samples), name = "RdYlGn")

color_map <- setNames(
  c("black", other_colors),
  c("UNpsiAR_0", other_samples)
)

p <- ggplot(data, aes(x = pos, y = mm_UtoC, group = sample, color = sample)) +
  geom_line(linewidth = 1) +
  scale_color_manual(values = color_map) +
  scale_y_continuous(labels = scales::percent_format(accuracy = 0.1)) +
  scale_x_continuous(
    breaks = c(8, 23, 38, 53, 68, 83, 98, 113),
    minor_breaks = seq(0, max(data$pos), by = 5)
  ) +
  labs(
    title = "U-to-C Mismatch Rate in UNUAR",
    subtitle = "Calculated as C Reads / (U + C Reads) for T positions",
    x = "Position",
    y = "Mismatch Rate",
    color = "Sample ID"
  ) +
  theme_minimal(base_family = "sans") +
  theme(
    legend.position = "bottom",
    plot.title = element_text(size = 18, hjust = 0.5),
    plot.subtitle = element_text(size = 12, hjust = 0.5),
    axis.title = element_text(size = 12, hjust = 0.5),
    axis.text = element_text(size = 12),
    panel.grid.major.x = element_line(color = "darkgray"),
    legend.key.spacing.x = unit(0.1, 'cm'),
    legend.key.spacing.y = unit(0.1, 'cm'),
    legend.key.width = unit(0.25, 'cm'),
    legend.title = element_blank()
  )

ggsave("plots/UtoC_mismatch.png", p, height = 5, width = 8, units = "in")
ggsave("plots/UtoC_mismatch.pdf", p, height = 4, width = 5, units = "in")

# ---- Calibration Curves - With Background Subtraction ----

motif_lookup <- tibble(
  pos =   c(8,   23,    38,    53,    68,    83,    98,    113),
  motif = c("UUUAA", "UUUAG", "UAUAA", "UAUAG", "UGUAA", "UGUAG", "UCUAA", "UCUAG")
)

plot_data_norm <- data %>%
  filter(pos %in% motif_lookup$pos) %>%
  left_join(motif_lookup, by = "pos") %>%
  mutate(motif = factor(motif, levels = motif_lookup$motif)) %>%
  group_by(pos) %>%
  mutate(mm_UtoC_corrected = (mm_UtoC - mm_UtoC[UNpsiAR_conc == 0])*100) %>%
  ungroup()

# --- Plots ---

p <- ggplot(plot_data_norm, aes(x = UNpsiAR_conc, y = mm_UtoC_corrected)) +
  geom_point(color = "dodgerblue", size = 3, alpha = 0.8) +
  geom_smooth(method = "lm", se = FALSE, color = "black", linetype = "dashed", linewidth = 0.8) +
  facet_wrap(~ motif, ncol = 4) +
  scale_y_continuous(limits = c(0, 100)) +
  scale_x_continuous(limits = c(0, 100)) +
  labs(
    title = "U-to-C Mismatch Rate by Motif",
    subtitle = "Corrected via Subtraction of UNpsiAR_0",
    x = "UNΨAR Concentration",
    y = "U-to-C Mismatch Rate (%)"
  ) +
  theme_minimal(base_family = "sans") +
  theme(
    plot.title = element_text(face = "bold", size = 20, hjust = 0.5),
    plot.subtitle = element_text(size = 12, hjust = 0.5),
    axis.title = element_text(size = 14),
    axis.text = element_text(size = 10),
    strip.text = element_text(face = "bold", size = 14),
    strip.background = element_rect(fill = "lightgrey", color = NA),
    panel.spacing = unit(1, "lines"),
    panel.border = element_rect(color = "black", fill = NA),
    legend.position = "none",
    aspect.ratio = 1
  )

ggsave("plots/calibration_curves_subtract.png", p, height = 6, width = 8, units = "in")
ggsave("plots/calibration_curves_subtract.pdf", p, height = 4, width = 6, units = "in")

# --- Formulas ---

# Calculate model coefficients for each motif
model_coeffs <- plot_data_norm %>%
  group_by(motif) %>%
  reframe(
    broom::tidy(lm(mm_UtoC_corrected ~ UNpsiAR_conc, data = pick(everything())))
  )

# Reshape data and create formula strings
formulas_df <- model_coeffs %>%
  select(motif, term, estimate) %>%
  pivot_wider(names_from = term, values_from = estimate) %>%
  dplyr::rename(
    intercept = "(Intercept)",
    slope = "UNpsiAR_conc"
  ) %>%
  mutate(
    equation_forward = paste0(
      "y = ", round(slope, 4), "x + ", round(intercept, 4)
    ),
    equation_deconv = paste0(
      "x = (y - (", round(intercept, 4), ")) / ", round(slope, 4)
    )
  )

# View and export the final dataframe
print(formulas_df)
write.csv(formulas_df, "derivation_formulas_subtract.csv", row.names = FALSE)
