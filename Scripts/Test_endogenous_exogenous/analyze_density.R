#!/usr/bin/env Rscript

# ------------------------------------------------------------
# Gene density test against BUSCO distribution
# Outputs:
#  - Density_test_candidates.tsv: per-candidate coding_pct, position wrt 5-95%, empirical p-value
#  - Density_test_busco_summary.tsv: BUSCO q5, median, q95
#  - gene_density_distribution_annotated.svg: histogram with quantile shading and candidate lines
# Comments in English. No emojis.
# ------------------------------------------------------------

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 3) {
  stop("Usage: Rscript gene_density_test.R <density_file.tsv> <busco_list.txt> <candidates_list.txt>")
}

density_file    <- args[1]
busco_file      <- args[2]
candidates_file <- args[3]

suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(readr)
  library(RColorBrewer)
  library(grid)
  library(scales)
  library(svglite)
})

# ---------- I/O and data ----------
output_dir <- dirname(density_file)

df <- read.table(density_file, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
# Expect columns: scaffold, coding_pct in [0,100]
required_cols <- c("scaffold", "coding_pct")
missing_cols <- setdiff(required_cols, names(df))
if (length(missing_cols) > 0) {
  stop(sprintf("Missing required columns in density file: %s", paste(missing_cols, collapse = ", ")))
}

busco_scaff     <- readLines(busco_file)
candidate_scaff <- readLines(candidates_file)

df <- df %>% mutate(coding_pct = as.numeric(coding_pct))
df_busco <- df %>% filter(scaffold %in% busco_scaff, !is.na(coding_pct))
df_candidates <- df %>% filter(scaffold %in% candidate_scaff, !is.na(coding_pct))

if (nrow(df_busco) == 0) stop("No BUSCO scaffolds found in density file.")

# ---------- Quantiles and ECDF-based test ----------
qs <- quantile(df_busco$coding_pct, probs = c(0.05, 0.50, 0.95), na.rm = TRUE, names = TRUE, type = 7)
q5  <- unname(qs[1]); q50 <- unname(qs[2]); q95 <- unname(qs[3])

# Empirical two-sided p-values per candidate using ECDF against BUSCO
busco_vals <- sort(df_busco$coding_pct)

ecdf_empirical_p <- function(x, ref_sorted) {
  n <- length(ref_sorted)
  r_le <- sum(ref_sorted <= x)
  r_gt <- sum(ref_sorted >  x)
  p_lower <- r_le / n
  p_upper <- r_gt / n
  p <- 2 * min(p_lower, p_upper)
  p <- min(1, max(0, p))
  return(p)
}

if (nrow(df_candidates) > 0) {
  df_candidates <- df_candidates %>%
    mutate(
      position_5_95 = case_when(
        coding_pct < q5  ~ "Below 5%",
        coding_pct > q95 ~ "Above 95%",
        TRUE             ~ "Inside 5-95%"
      ),
      empirical_p = vapply(coding_pct, ecdf_empirical_p, numeric(1), ref_sorted = busco_vals)
    )
}

# ---------- Exports ----------
# TSV 1: candidates
cand_out_path <- file.path(output_dir, "Density_test_candidates.tsv")
if (nrow(df_candidates) > 0) {
  cand_out <- df_candidates %>%
    select(scaffold, coding_pct, position_5_95, empirical_p)
  write.table(cand_out, file = cand_out_path, sep = "\t", row.names = FALSE, quote = FALSE)
} else {
  cand_out <- data.frame(scaffold = character(), coding_pct = numeric(),
                         position_5_95 = character(), empirical_p = numeric())
  write.table(cand_out, file = cand_out_path, sep = "\t", row.names = FALSE, quote = FALSE)
}

# TSV 2: BUSCO summary (q5, median, q95)
summary_out_path <- file.path(output_dir, "Density_test_busco_summary.tsv")
summary_df <- data.frame(
  n_busco = nrow(df_busco),
  q5 = q5,
  median = q50,
  q95 = q95
)
write.table(summary_df, file = summary_out_path, sep = "\t", row.names = FALSE, quote = FALSE)

# ---------- Plot (SVG) ----------
n_candidates <- length(unique(df_candidates$scaffold))
if (n_candidates == 0) {
  scaffold_colors <- c()
} else if (n_candidates <= 8) {
  palette_colors <- brewer.pal(n_candidates, "Set2")
  scaffold_colors <- setNames(palette_colors, unique(df_candidates$scaffold))
} else {
  hues <- hue_pal()(n_candidates)
  scaffold_colors <- setNames(hues, unique(df_candidates$scaffold))
}

p <- ggplot(df_busco, aes(x = coding_pct)) +
  # Histogram of BUSCO density
  geom_histogram(
    binwidth = 1,
    fill = alpha("grey40", 0.5),
    color = "grey30",
    linewidth = 0.4
  ) +
  # Shaded 5-95% interval
  annotate("rect", xmin = q5, xmax = q95, ymin = -Inf, ymax = Inf,
           alpha = 0.15, fill = "grey60") +
  # Quantile lines
  geom_vline(xintercept = q5,  linetype = "dashed", color = "grey20") +
  geom_vline(xintercept = q50, linetype = "solid",  color = "grey20") +
  geom_vline(xintercept = q95, linetype = "dashed", color = "grey20") +
  # Candidate vertical lines
  { if (nrow(df_candidates) > 0)
      geom_vline(data = df_candidates,
                 aes(xintercept = coding_pct, color = scaffold),
                 linetype = "dashed", linewidth = 1)
    else NULL } +
  scale_color_manual(name = "Candidate scaffolds", values = scaffold_colors) +
  guides(color = guide_legend(override.aes = list(linetype = "solid", shape = NA))) +
  scale_x_continuous(breaks = seq(0, 100, 10), limits = c(0, 100)) +
  labs(
    title = "Gene density distribution on BUSCO scaffolds",
    subtitle = sprintf("BUSCO q5=%.2f, median=%.2f, q95=%.2f", q5, q50, q95),
    x = "Percentage of coding bases per scaffold",
    y = "Number of scaffolds"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    legend.position = "right",
    legend.title = element_text(face = "bold"),
    legend.key.height = unit(0.6, "cm"),
    legend.text = element_text(size = 10),
    plot.title = element_text(hjust = 0.5),
    panel.grid.major = element_line(color = "grey85"),
    panel.grid.minor = element_blank()
  )

svg_path <- file.path(output_dir, "gene_density_distribution_annotated.svg")
ggsave(filename = svg_path, plot = p, width = 10, height = 6, dpi = 300, device = svglite)

# ---------- Console report ----------
cat("\n--- BUSCO quantiles ---\n")
cat(sprintf("q5    : %.4f\n", q5))
cat(sprintf("median: %.4f\n", q50))
cat(sprintf("q95   : %.4f\n", q95))

if (nrow(df_candidates) > 0) {
  cat("\n--- Candidates (per scaffold) ---\n")
  print(cand_out)
}

cat(sprintf("\nOutputs written to:\n  %s\n  %s\n  %s\n", cand_out_path, summary_out_path, svg_path))

