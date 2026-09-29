library(ggplot2)
library(dplyr)
library(stringr)

# ============================================================
# PATHS
# ============================================================

BASE <- "/beegfs/data/soukkal/Thesis/Hymenoptera_Project/PDV/Results/Viral_genes/IVSPER_Analysis"

old_file <- file.path(
  BASE,
  "xiphosomella_wirra_specimen2/ORF_prediction/scaffold3_ORFs_70AA_filtered.tsv"
)

new_file <- file.path(
  BASE,
  "Long_reads/07_Hybrid_assembly/ORF_prediction/XwIV_hybrid_final_ORFs.tsv"
)

output_pdf <- file.path(
  BASE,
  "Long_reads/07_Hybrid_assembly/ORF_prediction/XwIV_ORF_synteny.pdf"
)


# ============================================================
# PARAMETERS
# ============================================================

# In the previous scaffold3 assembly, the XwIV region
# started at position 9248.
old_viral_start <- 9248

# Genome lengths
old_genome_length <- 94426
new_genome_length <- 105371

# Vertical positions
y_old <- 2.0
y_new <- 0.8

# Height of ORF arrows
gene_height <- 0.32


# ============================================================
# LOAD PREVIOUS XwIV ORFs
# ============================================================

old <- read.delim(
  old_file,
  header = TRUE,
  sep = "\t",
  stringsAsFactors = FALSE,
  check.names = FALSE
)


# ============================================================
# PREPARE PREVIOUS XwIV ORFs
# ============================================================

old <- old %>%
  mutate(

    scaffold_start = pmin(start, end),
    scaffold_end   = pmax(start, end),

    # Position 9248 of scaffold3 becomes position 1 of XwIV
    start_plot = scaffold_start - old_viral_start + 1,
    end_plot   = scaffold_end   - old_viral_start + 1,

    # Common identifier
    old_ORF_ID = paste0("ORF", orf_number),

    # Number displayed on the plot
    label = as.character(orf_number)
  ) %>%

  filter(
    end_plot >= 1,
    start_plot <= old_genome_length
  ) %>%

  mutate(
    start_plot = pmax(start_plot, 1),
    end_plot   = pmin(end_plot, old_genome_length)
  )


# ============================================================
# LOAD NEW HYBRID XwIV ORFs
# ============================================================

new <- read.delim(
  new_file,
  header = TRUE,
  sep = "\t",
  stringsAsFactors = FALSE,
  check.names = FALSE
)


# ============================================================
# REVERSE THE NEW XwIV GENOME
# ============================================================
#
# The hybrid XwIV scaffold is in reverse orientation relative
# to the previous XwIV reconstruction.
#
# It is therefore reversed here for visualization.
# ============================================================

new <- new %>%
  mutate(

    original_start = pmin(start, end),
    original_end   = pmax(start, end),

    # Reverse coordinates
    start_plot = new_genome_length - original_end + 1,
    end_plot   = new_genome_length - original_start + 1,

    # Reverse ORF orientation
    plot_strand = ifelse(
      strand == "+",
      "-",
      "+"
    ),

    # XwIV_ORF_25 -> 25
    label = str_extract(
      New_ORF,
      "[0-9]+$"
    ),

    # Extract previous ORF identifier
    old_ORF_ID = ifelse(
      is.na(Previous_ORF) |
        Previous_ORF == "NA",
      NA_character_,
      str_extract(
        Previous_ORF,
        "ORF[0-9]+"
      )
    )
  )


# ============================================================
# IDENTIFY SHARED OLD ORFs
# ============================================================

shared_old_ORFs <- unique(
  new$old_ORF_ID[
    !is.na(new$old_ORF_ID)
  ]
)


# ============================================================
# CLASSIFY ORFs AS SHARED OR UNIQUE
# ============================================================
#
# Previous XwIV:
# unique = no corresponding ORF in hybrid XwIV
#
# Hybrid XwIV:
# unique = Previous_ORF is NA
# ============================================================

old <- old %>%
  mutate(
    status = ifelse(
      old_ORF_ID %in% shared_old_ORFs,
      "Shared",
      "Unique"
    )
  )


new <- new %>%
  mutate(
    status = ifelse(
      is.na(old_ORF_ID),
      "Unique",
      "Shared"
    )
  )


# ============================================================
# IDENTIFY HOMOLOGOUS ORF PAIRS
# ============================================================

links <- new %>%

  filter(!is.na(old_ORF_ID)) %>%

  inner_join(
    old %>%
      select(
        old_ORF_ID,
        old_start = start_plot,
        old_end = end_plot
      ),
    by = "old_ORF_ID"
  )


# ============================================================
# FUNCTION TO CREATE GENE ARROWS
# ============================================================

make_gene_arrows <- function(df, y, strand_column) {

  result <- list()

  for (i in seq_len(nrow(df))) {

    x1 <- df$start_plot[i]
    x2 <- df$end_plot[i]

    strand <- df[[strand_column]][i]

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

    y_bottom <- y - gene_height / 2
    y_top    <- y + gene_height / 2
    y_mid    <- y


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
            y_bottom,
            y_top,
            y_mid,
            y_bottom
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
            y_bottom,
            y_bottom,
            y_bottom,
            y_mid,
            y_top,
            y_top,
            y_top
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
            y_bottom,
            y_top,
            y_mid,
            y_bottom
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
            y_bottom,
            y_bottom,
            y_bottom,
            y_mid,
            y_top,
            y_top,
            y_top
          )
        )
      }
    }

    coords$gene <- i
    coords$status <- df$status[i]

    result[[i]] <- coords
  }

  bind_rows(result)
}


# ============================================================
# CREATE ORF ARROWS
# ============================================================

old_poly <- make_gene_arrows(
  old,
  y_old,
  "strand"
)

new_poly <- make_gene_arrows(
  new,
  y_new,
  "plot_strand"
)


# ============================================================
# LABEL POSITIONS
# ============================================================
#
# Labels are distributed over three vertical levels.
#
# Previous XwIV:
# above the genome.
#
# Hybrid XwIV:
# below the genome.
# ============================================================

old_labels <- old %>%
  mutate(
    x = (start_plot + end_plot) / 2
  ) %>%
  arrange(x) %>%
  mutate(

    label_level = rep(
      0:2,
      length.out = n()
    ),

    label_y = y_old +
      gene_height / 2 +
      0.15 +
      label_level * 0.16
  )


new_labels <- new %>%
  mutate(
    x = (start_plot + end_plot) / 2
  ) %>%
  arrange(x) %>%
  mutate(

    label_level = rep(
      0:2,
      length.out = n()
    ),

    label_y = y_new -
      gene_height / 2 -
      0.15 -
      label_level * 0.16
  )


# ============================================================
# CREATE HOMOLOGY RIBBONS
# ============================================================

ribbons <- links %>%
  rowwise() %>%
  do({

    data.frame(

      x = c(
        .$old_start,
        .$old_end,
        .$end_plot,
        .$start_plot
      ),

      y = c(
        y_old - gene_height / 2,
        y_old - gene_height / 2,
        y_new + gene_height / 2,
        y_new + gene_height / 2
      ),

      pair = .$old_ORF_ID
    )

  }) %>%
  ungroup()


# ============================================================
# BUILD PLOT
# ============================================================

p <- ggplot() +


  # ==========================================================
  # HOMOLOGY RIBBONS
  # ==========================================================

  geom_polygon(
    data = ribbons,
    aes(
      x = x,
      y = y,
      group = pair
    ),
    fill = "grey70",
    alpha = 0.30,
    linewidth = 0
  ) +


  # ==========================================================
  # GENOME BACKBONES
  # ==========================================================

  geom_segment(
    aes(
      x = 1,
      xend = old_genome_length,
      y = y_old,
      yend = y_old
    ),
    linewidth = 0.7
  ) +

  geom_segment(
    aes(
      x = 1,
      xend = new_genome_length,
      y = y_new,
      yend = y_new
    ),
    linewidth = 0.7
  ) +


  # ==========================================================
  # PREVIOUS XwIV ORFs
  # ==========================================================

  geom_polygon(
    data = old_poly,
    aes(
      x = x,
      y = y,
      group = gene,
      fill = status
    ),
    colour = "black",
    linewidth = 0.55
  ) +


  # ==========================================================
  # HYBRID XwIV ORFs
  # ==========================================================

  geom_polygon(
    data = new_poly,
    aes(
      x = x,
      y = y,
      group = gene,
      fill = status
    ),
    colour = "black",
    linewidth = 0.55
  ) +


  # ==========================================================
  # COLORS
  #
  # Shared ORFs = white
  # Unique ORFs = red
  # ==========================================================

  scale_fill_manual(
    name = NULL,
    values = c(
      "Shared" = "white",
      "Unique" = "red"
    ),
    breaks = c(
      "Shared",
      "Unique"
    ),
    labels = c(
      "Shared ORF",
      "Unique ORF"
    )
  ) +


  # ==========================================================
  # CONNECT OLD ORFs TO THEIR LABELS
  # ==========================================================

  geom_segment(
    data = old_labels,
    aes(
      x = x,
      xend = x,
      y = y_old + gene_height / 2,
      yend = label_y - 0.04
    ),
    linewidth = 0.20
  ) +


  # ==========================================================
  # OLD ORF NUMBERS
  # ==========================================================

  geom_text(
    data = old_labels,
    aes(
      x = x,
      y = label_y,
      label = label
    ),
    size = 2.6
  ) +


  # ==========================================================
  # CONNECT NEW ORFs TO THEIR LABELS
  # ==========================================================

  geom_segment(
    data = new_labels,
    aes(
      x = x,
      xend = x,
      y = y_new - gene_height / 2,
      yend = label_y + 0.04
    ),
    linewidth = 0.20
  ) +


  # ==========================================================
  # NEW ORF NUMBERS
  # ==========================================================

  geom_text(
    data = new_labels,
    aes(
      x = x,
      y = label_y,
      label = label
    ),
    size = 2.6
  ) +


  # ==========================================================
  # GENOME NAMES
  # ==========================================================

  annotate(
    "text",
    x = -5000,
    y = y_old,
    label = "Previous XwIV",
    hjust = 1,
    fontface = "bold",
    size = 4
  ) +

  annotate(
    "text",
    x = -5000,
    y = y_new,
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
      -15000,
      110000
    ),

    ylim = c(
      -0.05,
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

    # Put legend above figure
    legend.position = "top",

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
  "Previous XwIV ORFs:",
  nrow(old),
  "\n"
)

cat(
  "  Shared:",
  sum(old$status == "Shared"),
  "\n"
)

cat(
  "  Unique:",
  sum(old$status == "Unique"),
  "\n\n"
)

cat(
  "Hybrid XwIV ORFs:",
  nrow(new),
  "\n"
)

cat(
  "  Shared:",
  sum(new$status == "Shared"),
  "\n"
)

cat(
  "  Unique:",
  sum(new$status == "Unique"),
  "\n\n"
)

cat(
  "Homologous ORF pairs:",
  nrow(links),
  "\n"
)
