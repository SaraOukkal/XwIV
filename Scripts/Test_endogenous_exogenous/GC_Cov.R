#!/usr/bin/env Rscript

# ------------------------------------------------------------
# GC vs depth per scaffold with BUSCO and candidate highlighting
# - BUSCO in grey
# - Exact colors for specified candidate scaffolds
# - Vertical line: BUSCO median GC
# - Horizontal line: BUSCO median depth
# - Y axis limited to [0, 100]
# - Points: fill alpha 0.5 with opaque outline of the same color
# Comments in English. No emojis.
# ------------------------------------------------------------

suppressPackageStartupMessages({
  pkgs <- c("Biostrings","dplyr","readr","ggplot2","data.table","svglite")
  missing <- pkgs[!vapply(pkgs, requireNamespace, logical(1), quietly = TRUE)]
  if (length(missing) > 0) {
    if ("Biostrings" %in% missing) {
      if (!requireNamespace("BiocManager", quietly = TRUE)) {
        install.packages("BiocManager", repos = "http://cran.us.r-project.org")
      }
      BiocManager::install("Biostrings")
      missing <- setdiff(missing, "Biostrings")
    }
    if (length(missing) > 0) {
      install.packages(missing, repos = "http://cran.us.r-project.org")
    }
  }
  library(Biostrings)
  library(dplyr)
  library(readr)
  library(ggplot2)
  library(data.table)
  library(svglite)
})

# ----------------------- Functions --------------------------

compute_gc <- function(fasta_file) {
  # Return per-scaffold Length and GC percent
  fa <- readDNAStringSet(fasta_file)
  gc_counts <- rowSums(letterFrequency(fa, letters = c("G","C")))
  data.frame(
    Scaffold = names(fa),
    Length   = as.numeric(width(fa)),
    GC       = round(100 * gc_counts / as.numeric(width(fa)), 4),
    stringsAsFactors = FALSE
  )
}

compute_depth_stats <- function(depth_file) {
  # depth_file is a BED-like file with 3 columns: scaffold, position, depth
  df <- fread(depth_file, col.names = c("Scaffold","Pos","Depth"), showProgress = FALSE)
  df <- df[!is.na(Depth)]
  df %>%
    group_by(Scaffold) %>%
    summarise(
      Median_depth = median(Depth),
      Mean_depth   = mean(Depth),
      .groups = "drop"
    )
}

load_busco_scaffolds <- function(full_table_tsv) {
  # BUSCO full_table.tsv: V2 = Status, V3 = Contig
  # We mark scaffolds that contain at least one Complete or Fragmented BUSCO
  read.delim(full_table_tsv, comment.char = "#", header = FALSE, stringsAsFactors = FALSE) %>%
    dplyr::filter(V2 %in% c("Complete","Fragmented")) %>%
    dplyr::transmute(Scaffold = V3) %>%
    dplyr::distinct()
}

points_with_outline <- function(data, x, y, color_map) {
  # Two-layer points: semi-transparent fill and opaque outline in the same color
  ggplot(data, aes(x = {{x}}, y = {{y}})) +
    geom_point(
      aes(fill = ColorKey, size = Length),
      shape = 21, color = NA, alpha = 0.5, stroke = 0
    ) +
    geom_point(
      aes(color = ColorKey, size = Length),
      shape = 21, fill = NA, alpha = 1, stroke = 0.5
    ) +
    scale_fill_manual(values = color_map, drop = FALSE, guide = "none") +
    scale_color_manual(values = color_map, drop = FALSE, name = "Scaffolds") +
    scale_size_continuous(range = c(1.5, 6), name = "Scaffold length (bp)") +
    coord_cartesian(ylim = c(0, 100))
}

base_theme <- theme_minimal(base_size = 16) +
  theme(
    axis.title.x = element_text(size = 18),
    axis.title.y = element_text(size = 18),
    axis.text.x  = element_text(size = 14),
    axis.text.y  = element_text(size = 14),
    legend.title = element_text(size = 16),
    legend.text  = element_text(size = 14),
    plot.title   = element_text(size = 18, face = "bold")
  )

# -------------------- Sample parameters ---------------------

samples <- list(
  sp1 = list(
    name = "sp1",
    depth = "/beegfs/data/soukkal/Thesis/Hymenoptera_Project/PDV/Results/Viral_genes/IVSPER_Analysis/xiphosomella_wirra_specimen1/Mapping/xiphosomella_wirra:DHJPAR0041296.per-base.bed",
    busco = "/beegfs/project/horizon/data/stats/busco/specimens/DHJPAR0041296/run_insecta_odb10/full_table.tsv",
    candidates = c("scaffold524|size39189", "scaffold187|size52917", "scaffold33401|size845"),
    fasta = "/beegfs/project/horizon/data/assembly/specimens/DHJPAR0041296/redundans/scaffolds.reduced.fa",
    outdir = "/beegfs/data/soukkal/Thesis/Hymenoptera_Project/PDV/Results/Viral_genes/IVSPER_Analysis/xiphosomella_wirra_specimen1/GC_Cov/"
  ),
  sp2 = list(
    name = "sp2",
    depth = "/beegfs/data/soukkal/Thesis/Hymenoptera_Project/PDV/Results/Viral_genes/IVSPER_Analysis/xiphosomella_wirra_specimen2/Mapping/bed/DHJPAR0036296.per-base.bed",
    busco = "/beegfs/project/horizon/data/stats/busco/specimens/DHJPAR0036296/run_insecta_odb10/full_table.tsv",
    candidates = c("scaffold3|size100523"),
    fasta = "/beegfs/project/horizon/data/assembly/specimens/DHJPAR0036296/redundans/scaffolds.reduced.fa",
    outdir = "/beegfs/data/soukkal/Thesis/Hymenoptera_Project/PDV/Results/Viral_genes/IVSPER_Analysis/xiphosomella_wirra_specimen2/GC_Cov/"
  )
)

# ------------------------ Main loop -------------------------

for (sp in samples) {
  cat("\nProcessing:", sp$name, "\n")
  tryCatch({
    dir.create(sp$outdir, showWarnings = FALSE, recursive = TRUE)

    cat("  Reading FASTA...\n")
    gc_tab <- compute_gc(sp$fasta)

    cat("  Reading depth file...\n")
    depth_tab <- compute_depth_stats(sp$depth)

    cat("  Reading BUSCO data...\n")
    busco_tab <- load_busco_scaffolds(sp$busco)

    # Merge and annotate scaffold status
    full_tab <- gc_tab %>%
      left_join(depth_tab, by = "Scaffold") %>%
      mutate(
        Status = case_when(
          Scaffold %in% busco_tab$Scaffold     ~ "BUSCO",
          Scaffold %in% sp$candidates          ~ "Candidate",
          TRUE                                 ~ "Other"
        ),
        # Exact mapping for candidate color keys
        ColorKey = case_when(
          Status == "BUSCO"                                ~ "BUSCO",
          Scaffold == "scaffold524|size39189"              ~ "scaffold524",
          Scaffold == "scaffold187|size52917"              ~ "scaffold187",
          Scaffold == "scaffold33401|size845"              ~ "scaffold33401",
          Scaffold == "scaffold3|size100523"               ~ "scaffold3",
          TRUE                                             ~ "Other"
        )
      )

    # Color map as requested
    color_map <- c(
      "BUSCO"         = "grey50",
      "Other"         = "grey80",
      "scaffold3"     = "gold",
      "scaffold524"   = "dodgerblue3",
      "scaffold33401" = "forestgreen",
      "scaffold187"   = "red3"
    )

    # BUSCO reference medians
    busco_only <- full_tab %>% filter(Status == "BUSCO")
    median_busco_depth <- suppressWarnings(median(busco_only$Median_depth, na.rm = TRUE))
    median_busco_gc    <- suppressWarnings(median(busco_only$GC,           na.rm = TRUE))

    # Fallbacks if no BUSCO points available
    if (!is.finite(median_busco_depth)) {
      median_busco_depth <- median(full_tab$Median_depth, na.rm = TRUE)
    }
    if (!is.finite(median_busco_gc)) {
      median_busco_gc <- median(full_tab$GC, na.rm = TRUE)
    }

    # ---------- Plot 1: all scaffolds ----------
    cat("  Plotting: all scaffolds...\n")
    p1 <- points_with_outline(full_tab, GC, Median_depth, color_map) +
      geom_hline(yintercept = median_busco_depth, linetype = "dashed", color = "grey30") +
      geom_vline(xintercept = median_busco_gc,    linetype = "dotted", color = "grey30") +
      labs(
        title = paste0(sp$name, ": GC vs median depth (all scaffolds)"),
        x = "% GC",
        y = "Median depth (x)"
      ) +
      base_theme

    ggsave(
      filename = file.path(sp$outdir, paste0(sp$name, "_GC_vs_depth_all.svg")),
      plot = p1, width = 10, height = 7, device = svglite
    )

    # ---------- Plot 2: BUSCO + candidates ----------
    subset_tab <- full_tab %>%
      dplyr::filter(ColorKey %in% c("BUSCO","scaffold3","scaffold524","scaffold33401","scaffold187"))

    p2 <- points_with_outline(subset_tab, GC, Median_depth, color_map) +
      geom_hline(yintercept = median_busco_depth, linetype = "dashed", color = "grey30") +
      geom_vline(xintercept = median_busco_gc,    linetype = "dotted", color = "grey30") +
      labs(
        title = paste0(sp$name, ": GC vs median depth (BUSCO + candidates)"),
        x = "% GC",
        y = "Median depth (x)"
      ) +
      base_theme

    ggsave(
      filename = file.path(sp$outdir, paste0(sp$name, "_GC_vs_depth_candidates_busco_only.svg")),
      plot = p2, width = 10, height = 7, device = svglite
    )

    cat("  Plots saved to:", sp$outdir, "\n")

  }, error = function(e) {
    cat("  ERROR during processing of", sp$name, ":", conditionMessage(e), "\n")
  })
}
