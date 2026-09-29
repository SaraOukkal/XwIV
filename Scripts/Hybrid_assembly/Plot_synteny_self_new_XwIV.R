library(ggplot2)
library(dplyr)
library(stringr)

# ============================================================
# PATHS
# ============================================================

BASE <- "/beegfs/data/soukkal/Thesis/Hymenoptera_Project/PDV/Results/Viral_genes/IVSPER_Analysis"

self_file <- file.path(
  BASE,
  "Long_reads/07_Hybrid_assembly/Self_alignment/XwIV_self_alignment_filtered.tsv"
)

orf_file <- file.path(
  BASE,
  "Long_reads/07_Hybrid_assembly/ORF_prediction/XwIV_hybrid_final_ORFs.tsv"
)

output_pdf <- file.path(
  BASE,
  "Long_reads/07_Hybrid_assembly/Self_alignment/XwIV_self_alignment_plot.pdf"
)


# ============================================================
# PARAMETERS
# ============================================================

genome_length <- 105371

# Vertical positions
y_top <- 2.0
y_bottom <- 0.7

# ORF arrow height
gene_height <- 0.28


# ============================================================
# LOAD SELF-ALIGNMENT TABLE
# ============================================================
#
# XwIV_self_alignment_filtered.tsv has no header.
#
# MMseqs columns:
#
# 1  query
# 2  target
# 3  qstart
# 4  qend
# 5  tstart
# 6  tend
# 7  alnlen
# 8  pident
# 9  evalue
# 10 bits
#
# ============================================================

self <- read.delim(
  self_file,
  header = FALSE,
  sep = "\t",
  stringsAsFactors = FALSE,
  check.names = FALSE
)

# Check that the expected 10 columns are present
if (ncol(self) != 10) {

  stop(
    paste0(
      "Unexpected number of columns in self-alignment file: ",
      ncol(self),
      ". Expected 10 columns."
    )
  )
}

colnames(self) <- c(
  "query",
  "target",
  "qstart",
  "qend",
  "tstart",
  "tend",
  "alnlen",
  "pident",
  "evalue",
  "bits"
)


cat("\nSelf-alignment table:\n")
cat("Rows:", nrow(self), "\n")
cat("Columns:\n")
print(colnames(self))


# ============================================================
# CONVERT COLUMNS TO NUMERIC
# ============================================================

self <- self %>%
  mutate(

    qstart = as.numeric(qstart),
    qend   = as.numeric(qend),

    tstart = as.numeric(tstart),
    tend   = as.numeric(tend),

    alnlen = as.numeric(alnlen),
    pident = as.numeric(pident),
    evalue = as.numeric(evalue),
    bits   = as.numeric(bits)
  )


# ============================================================
# NORMALIZE ALIGNMENT COORDINATES
# ============================================================
#
# MMseqs coordinates can be reported in either orientation.
#
# For plotting we normalize them so:
#
# start < end
#
# ============================================================

self <- self %>%
  mutate(

    q_start = pmin(qstart, qend),
    q_end   = pmax(qstart, qend),

    t_start = pmin(tstart, tend),
    t_end   = pmax(tstart, tend)
  )


# ============================================================
# REMOVE TRIVIAL SELF-ALIGNMENTS
# ============================================================
#
# Remove cases where exactly the same genomic region
# aligns to itself.
# ============================================================

self <- self %>%
  filter(
    !(
      q_start == t_start &
      q_end == t_end
    )
  )


# ============================================================
# REMOVE RECIPROCAL DUPLICATES
# ============================================================
#
# Self-search can contain:
#
# Region A -> Region B
#
# and
#
# Region B -> Region A
#
# These describe the same homology.
#
# We create a canonical representation of each pair so that
# only one copy is retained.
# ============================================================

self <- self %>%
  rowwise() %>%
  mutate(

    first_is_query =
      q_start < t_start |
      (
        q_start == t_start &
        q_end <= t_end
      ),

    region1_start = ifelse(
      first_is_query,
      q_start,
      t_start
    ),

    region1_end = ifelse(
      first_is_query,
      q_end,
      t_end
    ),

    region2_start = ifelse(
      first_is_query,
      t_start,
      q_start
    ),

    region2_end = ifelse(
      first_is_query,
      t_end,
      q_end
    )
  ) %>%
  ungroup() %>%

  arrange(
    region1_start,
    region1_end,
    region2_start,
    region2_end,
    evalue,
    desc(bits)
  ) %>%

  distinct(
    region1_start,
    region1_end,
    region2_start,
    region2_end,
    .keep_all = TRUE
  )


cat(
  "\nNon-redundant homologous region pairs:",
  nrow(self),
  "\n"
)


# ============================================================
# LOAD ORFs
# ============================================================

orfs <- read.delim(
  orf_file,
  header = TRUE,
  sep = "\t",
  stringsAsFactors = FALSE,
  check.names = FALSE
)


# ============================================================
# PREPARE ORFs
# ============================================================

orfs <- orfs %>%
  mutate(

    start_plot = pmin(start, end),
    end_plot   = pmax(start, end),

    # XwIV_ORF_25 -> 25
    label = str_extract(
      New_ORF,
      "[0-9]+$"
    )
  )


cat(
  "ORFs:",
  nrow(orfs),
  "\n"
)


# ============================================================
# FUNCTION TO CREATE ORF ARROWS
# ============================================================

make_gene_arrows <- function(df, y) {

  result <- list()

  for (i in seq_len(nrow(df))) {

    x1 <- df$start_plot[i]
    x2 <- df$end_plot[i]

    strand <- df$strand[i]

    gene_length <- x2 - x1 + 1

    # Arrow-head size
    head_length <- min(
      350,
      gene_length * 0.30
    )

    head_length <- min(
      head_length,
      gene_length
    )

    y_bottom_gene <- y - gene_height / 2
    y_top_gene    <- y + gene_height / 2
    y_mid_gene    <- y


    # ========================================================
    # FORWARD STRAND
    # ========================================================

    if (strand == "+") {

      if (gene_length <= 400) {

        coords <- data.frame(

          x = c(
            x1,
            x1,
            x2,
            x1
          ),

          y = c(
            y_bottom_gene,
            y_top_gene,
            y_mid_gene,
            y_bottom_gene
          )
        )

      } else {

        coords <- data.frame(

          x = c(
            x1,
            x2 - head_length,
            x2 - head_length,
            x2,
            x2 - head_length,
            x2 - head_length,
            x1
          ),

          y = c(
            y_bottom_gene,
            y_bottom_gene,
            y_bottom_gene,
            y_mid_gene,
            y_top_gene,
            y_top_gene,
            y_top_gene
          )
        )
      }


    # ========================================================
    # REVERSE STRAND
    # ========================================================

    } else {

      if (gene_length <= 400) {

        coords <- data.frame(

          x = c(
            x2,
            x2,
            x1,
            x2
          ),

          y = c(
            y_bottom_gene,
            y_top_gene,
            y_mid_gene,
            y_bottom_gene
          )
        )

      } else {

        coords <- data.frame(

          x = c(
            x2,
            x1 + head_length,
            x1 + head_length,
            x1,
            x1 + head_length,
            x1 + head_length,
            x2
          ),

          y = c(
            y_bottom_gene,
            y_bottom_gene,
            y_bottom_gene,
            y_mid_gene,
            y_top_gene,
            y_top_gene,
            y_top_gene
          )
        )
      }
    }


    coords$gene <- i

    result[[i]] <- coords
  }

  bind_rows(result)
}


# ============================================================
# CREATE ORF ARROWS FOR BOTH GENOME COPIES
# ============================================================

orf_top_poly <- make_gene_arrows(
  orfs,
  y_top
)

orf_bottom_poly <- make_gene_arrows(
  orfs,
  y_bottom
)


# ============================================================
# ORF LABEL POSITIONS
# ============================================================
#
# Top genome:
# labels above
#
# Bottom genome:
# labels below
#
# Three levels are used to reduce overlap.
# ============================================================

top_labels <- orfs %>%
  mutate(
    x = (start_plot + end_plot) / 2
  ) %>%
  arrange(x) %>%
  mutate(

    label_level = rep(
      0:2,
      length.out = n()
    ),

    label_y = y_top +
      gene_height / 2 +
      0.14 +
      label_level * 0.15
  )


bottom_labels <- orfs %>%
  mutate(
    x = (start_plot + end_plot) / 2
  ) %>%
  arrange(x) %>%
  mutate(

    label_level = rep(
      0:2,
      length.out = n()
    ),

    label_y = y_bottom -
      gene_height / 2 -
      0.14 -
      label_level * 0.15
  )


# ============================================================
# CREATE HOMOLOGY RIBBONS
# ============================================================
#
# The first homologous region is placed on the top genome.
#
# The second homologous region is placed on the bottom genome.
#
# Ribbon width therefore corresponds to the genomic extent
# of the MMseqs nucleotide alignment.
# ============================================================

ribbon_list <- list()

for (i in seq_len(nrow(self))) {

  ribbon_list[[i]] <- data.frame(

    x = c(
      self$region1_start[i],
      self$region1_end[i],
      self$region2_end[i],
      self$region2_start[i]
    ),

    y = c(
      y_top - gene_height / 2,
      y_top - gene_height / 2,
      y_bottom + gene_height / 2,
      y_bottom + gene_height / 2
    ),

    alignment = i,

    pident = self$pident[i],

    alnlen = self$alnlen[i]
  )
}


ribbons <- bind_rows(
  ribbon_list
)


# ============================================================
# BUILD PLOT
# ============================================================

p <- ggplot() +


  # ==========================================================
  # HOMOLOGOUS REGIONS
  # ==========================================================

  geom_polygon(
    data = ribbons,
    aes(
      x = x,
      y = y,
      group = alignment,
      fill = pident
    ),
    alpha = 0.50,
    colour = NA
  ) +


  # ==========================================================
  # NUCLEOTIDE IDENTITY SCALE
  # ==========================================================

  scale_fill_gradient(
    name = "Nucleotide\nidentity (%)",
    low = "grey85",
    high = "darkred"
  ) +


  # ==========================================================
  # GENOME BACKBONES
  # ==========================================================

  geom_segment(
    aes(
      x = 1,
      xend = genome_length,
      y = y_top,
      yend = y_top
    ),
    linewidth = 0.7
  ) +

  geom_segment(
    aes(
      x = 1,
      xend = genome_length,
      y = y_bottom,
      yend = y_bottom
    ),
    linewidth = 0.7
  ) +


  # ==========================================================
  # TOP ORFs
  # ==========================================================

  geom_polygon(
    data = orf_top_poly,
    aes(
      x = x,
      y = y,
      group = gene
    ),
    fill = "white",
    colour = "black",
    linewidth = 0.5
  ) +


  # ==========================================================
  # BOTTOM ORFs
  # ==========================================================

  geom_polygon(
    data = orf_bottom_poly,
    aes(
      x = x,
      y = y,
      group = gene
    ),
    fill = "white",
    colour = "black",
    linewidth = 0.5
  ) +


  # ==========================================================
  # TOP LABEL CONNECTORS
  # ==========================================================

  geom_segment(
    data = top_labels,
    aes(
      x = x,
      xend = x,
      y = y_top + gene_height / 2,
      yend = label_y - 0.04
    ),
    linewidth = 0.18
  ) +


  # ==========================================================
  # TOP ORF NUMBERS
  # ==========================================================

  geom_text(
    data = top_labels,
    aes(
      x = x,
      y = label_y,
      label = label
    ),
    size = 2.4
  ) +


  # ==========================================================
  # BOTTOM LABEL CONNECTORS
  # ==========================================================

  geom_segment(
    data = bottom_labels,
    aes(
      x = x,
      xend = x,
      y = y_bottom - gene_height / 2,
      yend = label_y + 0.04
    ),
    linewidth = 0.18
  ) +


  # ==========================================================
  # BOTTOM ORF NUMBERS
  # ==========================================================

  geom_text(
    data = bottom_labels,
    aes(
      x = x,
      y = label_y,
      label = label
    ),
    size = 2.4
  ) +


  # ==========================================================
  # GENOME LABELS
  # ==========================================================

  annotate(
    "text",
    x = -4000,
    y = y_top,
    label = "Hybrid XwIV",
    hjust = 1,
    fontface = "bold",
    size = 4
  ) +

  annotate(
    "text",
    x = -4000,
    y = y_bottom,
    label = "Hybrid XwIV",
    hjust = 1,
    fontface = "bold",
    size = 4
  ) +


  # ==========================================================
  # X AXIS
  # ==========================================================

  scale_x_continuous(
    name = "Genomic position (bp)",

    breaks = seq(
      0,
      110000,
      by = 10000
    ),

    labels = scales::label_comma(),

    expand = c(0, 0)
  ) +


  # ==========================================================
  # PLOT LIMITS
  # ==========================================================

  coord_cartesian(
    xlim = c(
      -14000,
      110000
    ),

    ylim = c(
      -0.15,
      2.85
    ),

    clip = "off"
  ) +


  # ==========================================================
  # THEME
  # ==========================================================

  theme_classic(
    base_size = 12
  ) +

  theme(

    axis.title.y = element_blank(),

    axis.text.y = element_blank(),

    axis.ticks.y = element_blank(),

    axis.line.y = element_blank(),

    legend.position = "right",

    plot.margin = margin(
      15,
      15,
      15,
      20
    )
  )


# ============================================================
# DISPLAY
# ============================================================

print(p)


# ============================================================
# SAVE
# ============================================================

ggsave(
  filename = output_pdf,
  plot = p,
  width = 16,
  height = 6
)


# ============================================================
# SUMMARY
# ============================================================

cat(
  "\nFigure written to:\n",
  output_pdf,
  "\n\n"
)

cat(
  "Genome length:",
  genome_length,
  "bp\n"
)

cat(
  "ORFs plotted:",
  nrow(orfs),
  "\n"
)

cat(
  "Non-redundant homologous region pairs:",
  nrow(self),
  "\n"
)

if (nrow(self) > 0) {

  cat(
    "Alignment length range:",
    min(self$alnlen, na.rm = TRUE),
    "-",
    max(self$alnlen, na.rm = TRUE),
    "bp\n"
  )

  cat(
    "Identity range:",
    round(min(self$pident, na.rm = TRUE), 2),
    "-",
    round(max(self$pident, na.rm = TRUE), 2),
    "%\n"
  )
}
