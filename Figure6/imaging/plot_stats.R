library(dplyr)
library(stringr)
library(ggplot2)
library(scales)

# ---- Directory ----
active_dir <- "~/PUS7regulation2026/Figure6/imaging"

# ---- Colors ----
celltype293T <- "limegreen"
celltypeHepG2 <- "deepskyblue3"
cellline <- "#00a86b"

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
theme_boxplot <- function(base_size = 16) {
  theme_base_custom(base_size) +
    theme(legend.position = "none", panel.grid.major.x = element_blank(), panel.grid.minor = element_blank())
}
theme_bar <- function(base_size = 16) {
  theme_base_custom(base_size) +
    theme(legend.position = "none", panel.grid.major.x = element_blank(), panel.grid.minor = element_blank())
}

# ---- Datasets ----

all_data <- read.csv(file.path(active_dir, "rep2_data.csv"))

data <- all_data %>%
    filter(Group_ID %in% c("HEK293T_P3", "HepG2_P3", "HEK293T_P4", "HepG2_P4")) %>%
    mutate(labels = case_when(
        Group_ID == "HEK293T_P3" ~ "HEK293T PUS7 OE",
        Group_ID == "HepG2_P3" ~ "HepG2 PUS7 OE",
        Group_ID == "HEK293T_P4" ~ "HEK293T WT",
        Group_ID == "HepG2_P4" ~ "HepG2 WT"
    )) %>%
    mutate(labels = factor(labels, levels = c("HepG2 WT", "HEK293T WT", "HepG2 PUS7 OE", "HEK293T PUS7 OE"))) %>%
    mutate(color = case_when(
        str_detect(labels, "HEK293T") ~ celltype293T,
        str_detect(labels, "HepG2") ~ celltypeHepG2
    ))

# Extract replicate number from Image_Name (the number following 'rep')
data <- data %>%
  mutate(replicate = as.integer(str_extract(Image_Name, "(?<=rep)\\d+")))


# ---- Plot ----

# Downsample to standardize 
num <- 100

downsample_data <- data %>%
    group_by(labels) %>%
    slice_sample(n = num) %>%
    ungroup()

p <- ggplot(downsample_data, aes(x = labels, y = log2_ratio_nuc, fill = color)) +
    geom_violin(alpha = 0.7) +
    geom_boxplot(width = 0.2, size = 0.3) +
    scale_fill_identity() +
    theme_boxplot(base_size = 8) +
    labs(x = "Cell Line", y = "log2(Nuclear/Cytoplasmic Ratio)") +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    scale_y_continuous(breaks = pretty_breaks(n = 5))

plot_path <- file.path(active_dir, "violin_IF_quantification")
ggsave(p, filename = paste0(plot_path, ".png"), width = 2, height = 2.5)
ggsave(p, filename = paste0(plot_path, ".pdf"), width = 1.5, height = 1.5)


# ---- Statistical tests comparing HepG2 vs HEK293T (WT and OE) ----
# Use Wilcoxon rank-sum (non-parametric) as the primary test
# also report Welch's t-test for reference.
comparisons <- list(
    WT = c("HepG2 WT", "HEK293T WT"),
    OE = c("HepG2 PUS7 OE", "HEK293T PUS7 OE")
)

test_results <- lapply(names(comparisons), function(nm) {
    labs <- comparisons[[nm]]
    df <- downsample_data %>% filter(labels %in% labs)
    df$labels <- factor(df$labels, levels = labs)
    wil <- wilcox.test(log2_ratio_nuc ~ labels, data = df, exact = FALSE)
    ttt <- t.test(log2_ratio_nuc ~ labels, data = df)
    data.frame(
        comparison = nm,
        group1 = labs[1],
        group2 = labs[2],
        wilcox_p = wil$p.value,
        ttest_p = ttt$p.value,
        median1 = median(df$log2_ratio_nuc[df$labels == labs[1]], na.rm = TRUE),
        median2 = median(df$log2_ratio_nuc[df$labels == labs[2]], na.rm = TRUE),
        stringsAsFactors = FALSE
    )
}) %>% bind_rows()

test_results <- test_results %>%
    mutate(wilcox_p_adj = p.adjust(wilcox_p, method = "bonferroni"),
                 ttest_p_adj = p.adjust(ttest_p, method = "bonferroni"))

print(test_results)
write.csv(test_results, file = paste0(plot_path, "_stats.csv"), row.names = FALSE)