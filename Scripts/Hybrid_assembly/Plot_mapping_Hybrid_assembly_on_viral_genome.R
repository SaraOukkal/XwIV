#!/beegfs/data/soft/R-4.5.3/bin/Rscript

# ============================================================
# XwIV vs hybrid assembly - alignment visualisation
#
# Plot 1:
#   All minimap2 alignments
#
# Plot 2:
#   Only hybrid scaffolds for which >= 90% of the complete
#   scaffold is covered by XwIV alignments
#
# No minimum alignment-length filter.
#
# Colour:
#   % of complete hybrid scaffold covered by XwIV alignments
#
# Opacity:
#   nucleotide identity of each individual alignment
#
# Output:
#   PDF plots + TSV summary tables
# ============================================================


# ============================================================
# 0. R LIBRARY AND REQUIRED PACKAGES
# ============================================================

user_lib <- "/beegfs/home/soukkal/R/x86_64-pc-linux-gnu-library/4.5"

dir.create(
  user_lib,
  recursive = TRUE,
  showWarnings = FALSE
)

.libPaths(c(
  user_lib,
  .libPaths()
))

required_packages <- c(
  "dplyr",
  "ggplot2",
  "scales"
)

for (pkg in required_packages) {

  if (!requireNamespace(pkg, quietly = TRUE)) {

    message("Installing missing package: ", pkg)

    install.packages(
      pkg,
      repos = "https://cloud.r-project.org",
      lib = user_lib
    )
  }

  library(
    pkg,
    character.only = TRUE
  )
}


# ============================================================
# 1. INPUT
# ============================================================

paf_file <- paste0(
  "/beegfs/data/soukkal/Thesis/Hymenoptera_Project/PDV/Results/",
  "Viral_genes/IVSPER_Analysis/Long_reads/07_Hybrid_assembly/",
  "Viral_search/Mapping_Hybrid_on_XwIV/XwIV_vs_Hybrid_assembly.paf"
)

outdir <- dirname(paf_file)

cat("\nInput PAF:\n")
cat(paf_file, "\n\n")

cat("Output directory:\n")
cat(outdir, "\n\n")


# ============================================================
# 2. READ PAF
# ============================================================

paf <- read.delim(
  paf_file,
  header = FALSE,
  sep = "\t",
  comment.char = "",
  fill = TRUE,
  stringsAsFactors = FALSE
)

if (ncol(paf) < 12) {
  stop("The PAF file contains fewer than the 12 mandatory columns.")
}

colnames(paf)[1:12] <- c(
  "scaffold",
  "scaffold_length",
  "qstart",
  "qend",
  "strand",
  "virus",
  "virus_length",
  "virus_start",
  "virus_end",
  "matches",
  "alignment_length",
  "mapq"
)


# ============================================================
# 3. CALCULATE ALIGNMENT STATISTICS
# ============================================================

paf <- paf %>%
  mutate(

    identity =
      100 * matches / alignment_length,

    query_span =
      qend - qstart,

    virus_span =
      virus_end - virus_start,

    individual_scaffold_coverage =
      100 * query_span / scaffold_length
  )


# ============================================================
# 4. BASIC SUMMARY
# ============================================================

cat("====================================================\n")
cat("COMPLETE PAF SUMMARY\n")
cat("====================================================\n\n")

cat(
  "Number of alignments: ",
  nrow(paf),
  "\n",
  sep = ""
)

cat(
  "Number of hybrid scaffolds: ",
  n_distinct(paf$scaffold),
  "\n",
  sep = ""
)

cat(
  "Viral genome length: ",
  paste(unique(paf$virus_length), collapse = ", "),
  " bp\n\n",
  sep = ""
)

cat("Alignment length summary:\n")
print(summary(paf$alignment_length))

cat("\nIdentity summary (%):\n")
print(summary(paf$identity))

cat("\n")


# ============================================================
# 5. FUNCTION TO MERGE OVERLAPPING INTERVALS
# ============================================================

merged_interval_length <- function(start, end) {

  intervals <- data.frame(
    start = start,
    end = end
  )

  intervals <- intervals[
    order(intervals$start, intervals$end),
  ]

  current_start <- intervals$start[1]
  current_end <- intervals$end[1]

  total_length <- 0

  if (nrow(intervals) > 1) {

    for (i in 2:nrow(intervals)) {

      if (intervals$start[i] <= current_end) {

        current_end <- max(
          current_end,
          intervals$end[i]
        )

      } else {

        total_length <-
          total_length +
          (current_end - current_start)

        current_start <- intervals$start[i]
        current_end <- intervals$end[i]
      }
    }
  }

  total_length <-
    total_length +
    (current_end - current_start)

  return(total_length)
}


# ============================================================
# 6. FUNCTION TO CALCULATE SCAFFOLD COVERAGE
# ============================================================

calculate_scaffold_coverage <- function(data) {

  data %>%
    group_by(scaffold) %>%
    summarise(

      scaffold_length =
        first(scaffold_length),

      number_alignments =
        n(),

      first_viral_position =
        min(virus_start),

      last_viral_position =
        max(virus_end),

      longest_alignment =
        max(alignment_length),

      viral_bp_on_scaffold =
        merged_interval_length(
          qstart,
          qend
        ),

      scaffold_viral_coverage =
        100 *
        viral_bp_on_scaffold /
        scaffold_length,

      mean_identity =
        weighted.mean(
          identity,
          alignment_length
        ),

      max_identity =
        max(identity),

      .groups = "drop"
    ) %>%
    arrange(
      first_viral_position,
      desc(longest_alignment)
    )
}


# ============================================================
# 7. CALCULATE COVERAGE USING ALL ALIGNMENTS
# ============================================================

coverage_all <-
  calculate_scaffold_coverage(paf)

cat("====================================================\n")
cat("SCAFFOLD COVERAGE SUMMARY\n")
cat("====================================================\n\n")

print(coverage_all)

cat("\n")

write.table(
  coverage_all,
  file = file.path(
    outdir,
    "XwIV_Hybrid_scaffold_coverage_ALL.tsv"
  ),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)


# ============================================================
# 8. IDENTIFY SCAFFOLDS >= 90% COVERED BY XwIV
# ============================================================

coverage_threshold <- 90

scaffolds_90 <- coverage_all %>%
  filter(
    scaffold_viral_coverage >= coverage_threshold
  ) %>%
  pull(scaffold)

paf_90 <- paf %>%
  filter(
    scaffold %in% scaffolds_90
  )

cat("====================================================\n")
cat("SCAFFOLDS >= 90% COVERED BY XwIV\n")
cat("====================================================\n\n")

cat(
  "Number of scaffolds >= ",
  coverage_threshold,
  "% covered: ",
  length(scaffolds_90),
  "\n\n",
  sep = ""
)

coverage_all %>%
  filter(
    scaffold_viral_coverage >= coverage_threshold
  ) %>%
  arrange(
    desc(scaffold_viral_coverage)
  ) %>%
  print()

cat("\n")

write.table(
  coverage_all %>%
    filter(
      scaffold_viral_coverage >= coverage_threshold
    ) %>%
    arrange(
      desc(scaffold_viral_coverage)
    ),
  file = file.path(
    outdir,
    "XwIV_Hybrid_scaffolds_coverage90.tsv"
  ),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)


# ============================================================
# 9. FUNCTION TO CREATE PLOT
# ============================================================

make_alignment_plot <- function(
  data,
  coverage_table,
  plot_name,
  plot_title,
  plot_subtitle
) {

  if (nrow(data) == 0) {

    warning(
      "No alignments available for ",
      plot_name
    )

    return(NULL)
  }

  scaffold_summary <- coverage_table %>%
    filter(
      scaffold %in% unique(data$scaffold)
    )

  data <- data %>%
    left_join(
      scaffold_summary %>%
        select(
          scaffold,
          viral_bp_on_scaffold,
          scaffold_viral_coverage
        ),
      by = "scaffold"
    )

  alignment_table <- data %>%
    select(
      scaffold,
      scaffold_length,
      qstart,
      qend,
      strand,
      virus,
      virus_length,
      virus_start,
      virus_end,
      matches,
      alignment_length,
      individual_scaffold_coverage,
      viral_bp_on_scaffold,
      scaffold_viral_coverage,
      identity,
      mapq
    ) %>%
    arrange(
      desc(alignment_length)
    )

  write.table(
    alignment_table,
    file = file.path(
      outdir,
      paste0(
        "XwIV_Hybrid_alignment_table_",
        plot_name,
        ".tsv"
      )
    ),
    sep = "\t",
    quote = FALSE,
    row.names = FALSE
  )

  scaffold_order <- scaffold_summary %>%
    arrange(
      first_viral_position,
      desc(longest_alignment)
    ) %>%
    pull(scaffold)

  data <- data %>%
    mutate(
      scaffold = factor(
        scaffold,
        levels = rev(scaffold_order)
      )
    )

  cat("\n====================================================\n")
  cat(plot_title, "\n")
  cat("====================================================\n\n")

  cat(
    "Number of alignments: ",
    nrow(data),
    "\n",
    sep = ""
  )

  cat(
    "Number of hybrid scaffolds: ",
    n_distinct(data$scaffold),
    "\n\n",
    sep = ""
  )

  p <- ggplot(data) +

    geom_segment(
      aes(
        x = 0,
        xend = virus_length,
        y = scaffold,
        yend = scaffold
      ),
      linewidth = 1.5,
      colour = "grey85"
    ) +

    geom_segment(
      aes(
        x = virus_start,
        xend = virus_end,
        y = scaffold,
        yend = scaffold,
        colour = scaffold_viral_coverage,
        alpha = identity
      ),
      linewidth = 4,
      lineend = "butt"
    ) +

    scale_colour_gradient(
      name =
        "Hybrid scaffold covered\nby XwIV (%)",
      low = "steelblue",
      high = "red",
      limits = c(0, 100)
    ) +

    scale_alpha_continuous(
      name =
        "Alignment\nidentity (%)",
      range = c(0.15, 1),
      limits = c(0, 100)
    ) +

    scale_x_continuous(
      labels =
        scales::label_number(
          scale = 1 / 1000,
          suffix = " kb"
        ),
      expand =
        expansion(
          mult = c(0.01, 0.01)
        )
    ) +

    labs(
      x = "Position on XwIV genome",
      y = "Hybrid scaffold",
      title = plot_title,
      subtitle = plot_subtitle
    ) +

    theme_classic(
      base_size = 12
    ) +

    theme(

      axis.text.y =
        element_text(
          size = 7
        ),

      axis.text.x =
        element_text(
          size = 10
        ),

      axis.title =
        element_text(
          size = 12
        ),

      plot.title =
        element_text(
          size = 14,
          face = "bold"
        ),

      plot.subtitle =
        element_text(
          size = 9
        ),

      legend.position =
        "right"
    )

  print(p)

  plot_height <- max(
    6,
    0.25 *
      n_distinct(data$scaffold)
  )

  pdf_file <- file.path(
    outdir,
    paste0(
      "XwIV_Hybrid_scaffold_alignment_",
      plot_name,
      ".pdf"
    )
  )

  ggsave(
    filename = pdf_file,
    plot = p,
    device = "pdf",
    width = 12,
    height = plot_height,
    limitsize = FALSE
  )

  cat(
    "\nPDF saved:\n",
    pdf_file,
    "\n",
    sep = ""
  )

  return(
    list(
      plot = p,
      summary = scaffold_summary,
      alignments = alignment_table
    )
  )
}


# ============================================================
# 10. PLOT 1 - ALL ALIGNMENTS
# ============================================================

result_all <- make_alignment_plot(

  data = paf,

  coverage_table = coverage_all,

  plot_name = "ALL",

  plot_title =
    "Hybrid assembly scaffolds aligned against the XwIV genome - all alignments",

  plot_subtitle =
    "No alignment length or identity filtering | Colour = % of hybrid scaffold covered by XwIV | Opacity = nucleotide identity"
)


# ============================================================
# 11. PLOT 2 - SCAFFOLDS >= 90% COVERED BY XwIV
# ============================================================

result_90 <- make_alignment_plot(

  data = paf_90,

  coverage_table = coverage_all,

  plot_name = "SCAFFOLD_COVERAGE90",

  plot_title =
    "Hybrid scaffolds with >= 90% of their sequence covered by XwIV",

  plot_subtitle =
    "Only scaffolds >= 90% covered by XwIV | Colour = % of hybrid scaffold covered by XwIV | Opacity = nucleotide identity"
)


# ============================================================
# 12. FINAL SUMMARY
# ============================================================

cat("\n====================================================\n")
cat("ANALYSIS COMPLETED\n")
cat("====================================================\n\n")

cat("ALL ALIGNMENTS\n")

cat(
  "  Alignments: ",
  nrow(paf),
  "\n",
  sep = ""
)

cat(
  "  Hybrid scaffolds: ",
  n_distinct(paf$scaffold),
  "\n\n",
  sep = ""
)

cat("SCAFFOLD COVERAGE >= 90%\n")

cat(
  "  Alignments: ",
  nrow(paf_90),
  "\n",
  sep = ""
)

cat(
  "  Hybrid scaffolds: ",
  n_distinct(paf_90$scaffold),
  "\n\n",
  sep = ""
)

cat("Generated plots:\n\n")

cat(
  file.path(
    outdir,
    "XwIV_Hybrid_scaffold_alignment_ALL.pdf"
  ),
  "\n"
)

cat(
  file.path(
    outdir,
    "XwIV_Hybrid_scaffold_alignment_SCAFFOLD_COVERAGE90.pdf"
  ),
  "\n"
)

cat("\nGenerated tables:\n\n")

cat(
  file.path(
    outdir,
    "XwIV_Hybrid_scaffold_coverage_ALL.tsv"
  ),
  "\n"
)

cat(
  file.path(
    outdir,
    "XwIV_Hybrid_scaffolds_coverage90.tsv"
  ),
  "\n"
)

cat(
  file.path(
    outdir,
    "XwIV_Hybrid_alignment_table_ALL.tsv"
  ),
  "\n"
)

cat(
  file.path(
    outdir,
    "XwIV_Hybrid_alignment_table_SCAFFOLD_COVERAGE90.tsv"
  ),
  "\n"
)

cat("\nDone.\n")
