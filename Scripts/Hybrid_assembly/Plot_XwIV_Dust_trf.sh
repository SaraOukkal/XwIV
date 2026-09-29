#!/bin/bash

#SBATCH --job-name=XwIV_complexity_plot
#SBATCH --output=%x_%j.out
#SBATCH --error=%x_%j.err
#SBATCH --time=00:30:00
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --exclude=pbil-deb[14-27]

set -euo pipefail


# ============================================================
# PATHS
# ============================================================

OUTPUT_DIR="/beegfs/data/soukkal/Thesis/Hymenoptera_Project/PDV/Results/Viral_genes/IVSPER_Analysis/Long_reads/07_Hybrid_assembly/Terminal_complexity"

DUST="${OUTPUT_DIR}/XwIV_DUST_regions.tsv"

TRF_RAW="${OUTPUT_DIR}/XwIV_hybrid_genome.fa.2.7.7.80.10.50.500.dat"

TRF_TSV="${OUTPUT_DIR}/XwIV_TRF_regions.tsv"

GENOME="${OUTPUT_DIR}/XwIV_hybrid_genome.fa"
FAI="${GENOME}.fai"

R_SCRIPT="${OUTPUT_DIR}/plot_XwIV_complexity.R"

PLOT_REGIONS="${OUTPUT_DIR}/XwIV_low_complexity_tandem_repeats.pdf"

PLOT_WINDOWS="${OUTPUT_DIR}/XwIV_low_complexity_tandem_repeats_1kb_windows.pdf"

WINDOW_TABLE="${OUTPUT_DIR}/XwIV_complexity_1kb_windows.tsv"


# ============================================================
# CHECK FILES
# ============================================================

for FILE in "${DUST}" "${TRF_RAW}" "${FAI}"; do

    if [[ ! -s "${FILE}" ]]; then
        echo "ERROR: file missing or empty:"
        echo "${FILE}"
        exit 1
    fi

done


GENOME_LENGTH=$(cut -f2 "${FAI}" | head -n 1)

echo
echo "XwIV genome length: ${GENOME_LENGTH} bp"


# ============================================================
# CONVERT TRF OUTPUT
# ============================================================

echo
echo "============================================"
echo "Converting TRF output"
echo "============================================"


echo -e "start\tend\tlength\tperiod_size\tcopy_number\tconsensus_size\tpercent_matches\tpercent_indels\tscore\tA_percent\tC_percent\tG_percent\tT_percent\tentropy\tconsensus_sequence" > "${TRF_TSV}"


awk '
    /^[0-9]/ {

        start=$1
        end=$2
        repeat_length=end-start+1

        print start "\t" \
              end "\t" \
              repeat_length "\t" \
              $3 "\t" \
              $4 "\t" \
              $5 "\t" \
              $6 "\t" \
              $7 "\t" \
              $8 "\t" \
              $9 "\t" \
              $10 "\t" \
              $11 "\t" \
              $12 "\t" \
              $13 "\t" \
              $14
    }
' "${TRF_RAW}" >> "${TRF_TSV}"


echo
echo "Number of TRF repeats:"
awk 'NR>1' "${TRF_TSV}" | wc -l


# ============================================================
# CREATE R SCRIPT
# ============================================================

cat > "${R_SCRIPT}" <<EOF

library(ggplot2)
library(dplyr)

genome_length <- ${GENOME_LENGTH}

dust_file <- "${DUST}"
trf_file <- "${TRF_TSV}"

plot_regions_file <- "${PLOT_REGIONS}"
plot_windows_file <- "${PLOT_WINDOWS}"

window_table_file <- "${WINDOW_TABLE}"

window_size <- 1000


# ============================================================
# LOAD DATA
# ============================================================

dust <- read.delim(
    dust_file,
    header = TRUE,
    sep = "\t",
    stringsAsFactors = FALSE
)

trf <- read.delim(
    trf_file,
    header = TRUE,
    sep = "\t",
    stringsAsFactors = FALSE
)


# ============================================================
# BASIC CHECKS
# ============================================================

dust <- dust %>%
    filter(
        !is.na(start),
        !is.na(end),
        start <= end
    )

trf <- trf %>%
    filter(
        !is.na(start),
        !is.na(end),
        start <= end
    )


# ============================================================
# FIGURE 1
# EXACT POSITIONS OF DUST AND TRF REGIONS
# ============================================================

dust_plot <- dust %>%
    transmute(
        start = start,
        end = end,
        track = "Low complexity (DUST)"
    )


trf_plot <- trf %>%
    transmute(
        start = start,
        end = end,
        track = "Tandem repeats (TRF)"
    )


features <- bind_rows(
    dust_plot,
    trf_plot
)


features\$track <- factor(
    features\$track,
    levels = c(
        "Low complexity (DUST)",
        "Tandem repeats (TRF)"
    )
)


features <- features %>%
    mutate(
        y = case_when(
            track == "Low complexity (DUST)" ~ 2,
            track == "Tandem repeats (TRF)" ~ 1
        )
    )


backbone <- data.frame(
    y = c(2, 1)
)


p_regions <- ggplot() +

    geom_segment(
        data = backbone,
        aes(
            x = 1,
            xend = genome_length,
            y = y,
            yend = y
        ),
        linewidth = 0.5
    ) +

    geom_rect(
        data = features,
        aes(
            xmin = start,
            xmax = end,
            ymin = y - 0.22,
            ymax = y + 0.22,
            fill = track
        ),
        linewidth = 0
    ) +

    scale_x_continuous(
        limits = c(
            0,
            genome_length
        ),
        breaks = seq(
            0,
            100000,
            by = 10000
        ),
        labels = function(x) x / 1000,
        expand = c(
            0,
            0
        )
    ) +

    scale_y_continuous(
        limits = c(
            0.5,
            2.5
        ),
        breaks = c(
            1,
            2
        ),
        labels = c(
            "Tandem repeats (TRF)",
            "Low complexity (DUST)"
        )
    ) +

    labs(
        title = "Sequence complexity and tandem repeats across the XwIV genome",
        x = "Position on XwIV genome (kb)",
        y = NULL
    ) +

    guides(
        fill = "none"
    ) +

    theme_classic(
        base_size = 12
    ) +

    theme(
        axis.line.y = element_blank(),
        axis.ticks.y = element_blank(),

        axis.text.y = element_text(
            size = 11
        ),

        plot.title = element_text(
            face = "bold",
            hjust = 0.5
        ),

        plot.margin = margin(
            15,
            20,
            15,
            20
        )
    )


ggsave(
    plot_regions_file,
    p_regions,
    width = 12,
    height = 3.5
)


# ============================================================
# FUNCTION TO MERGE OVERLAPPING INTERVALS
# ============================================================

merge_intervals <- function(start, end) {

    if (length(start) == 0) {
        return(
            data.frame(
                start = numeric(0),
                end = numeric(0)
            )
        )
    }

    x <- data.frame(
        start = start,
        end = end
    )

    x <- x[
        order(
            x\$start,
            x\$end
        ),
    ]

    merged_start <- c()
    merged_end <- c()

    current_start <- x\$start[1]
    current_end <- x\$end[1]

    if (nrow(x) > 1) {

        for (i in 2:nrow(x)) {

            if (x\$start[i] <= current_end + 1) {

                current_end <- max(
                    current_end,
                    x\$end[i]
                )

            } else {

                merged_start <- c(
                    merged_start,
                    current_start
                )

                merged_end <- c(
                    merged_end,
                    current_end
                )

                current_start <- x\$start[i]
                current_end <- x\$end[i]
            }
        }
    }

    merged_start <- c(
        merged_start,
        current_start
    )

    merged_end <- c(
        merged_end,
        current_end
    )

    data.frame(
        start = merged_start,
        end = merged_end
    )
}


# ============================================================
# MERGE OVERLAPPING DUST / TRF INTERVALS
# ============================================================

dust_merged <- merge_intervals(
    dust\$start,
    dust\$end
)

trf_merged <- merge_intervals(
    trf\$start,
    trf\$end
)


# ============================================================
# CREATE 1-KB WINDOWS
# ============================================================

windows <- data.frame(
    start = seq(
        1,
        genome_length,
        by = window_size
    )
)

windows\$end <- pmin(
    windows\$start + window_size - 1,
    genome_length
)

windows\$window_length <- (
    windows\$end -
    windows\$start +
    1
)

windows\$midpoint <- (
    windows\$start +
    windows\$end
) / 2


# ============================================================
# FUNCTION TO CALCULATE UNIQUE BP COVERAGE
# ============================================================

calculate_coverage <- function(
    window_start,
    window_end,
    intervals
) {

    if (nrow(intervals) == 0) {
        return(0)
    }

    overlap_start <- pmax(
        window_start,
        intervals\$start
    )

    overlap_end <- pmin(
        window_end,
        intervals\$end
    )

    overlap_length <- pmax(
        0,
        overlap_end -
        overlap_start +
        1
    )

    sum(overlap_length)
}


# ============================================================
# CALCULATE COVERAGE PER WINDOW
# ============================================================

windows\$DUST_bp <- mapply(
    calculate_coverage,
    windows\$start,
    windows\$end,
    MoreArgs = list(
        intervals = dust_merged
    )
)


windows\$TRF_bp <- mapply(
    calculate_coverage,
    windows\$start,
    windows\$end,
    MoreArgs = list(
        intervals = trf_merged
    )
)


windows\$DUST_percent <- (
    windows\$DUST_bp /
    windows\$window_length
) * 100


windows\$TRF_percent <- (
    windows\$TRF_bp /
    windows\$window_length
) * 100


# ============================================================
# SAVE WINDOW TABLE
# ============================================================

write.table(
    windows,
    file = window_table_file,
    sep = "\t",
    quote = FALSE,
    row.names = FALSE
)


# ============================================================
# PREPARE WINDOW DATA FOR PLOT
# ============================================================

dust_windows <- windows %>%
    transmute(
        start = start,
        end = end,
        midpoint = midpoint,
        percent = DUST_percent,
        track = "Low complexity (DUST)"
    )


trf_windows <- windows %>%
    transmute(
        start = start,
        end = end,
        midpoint = midpoint,
        percent = TRF_percent,
        track = "Tandem repeats (TRF)"
    )


window_plot_data <- bind_rows(
    dust_windows,
    trf_windows
)


window_plot_data\$track <- factor(
    window_plot_data\$track,
    levels = c(
        "Low complexity (DUST)",
        "Tandem repeats (TRF)"
    )
)


# ============================================================
# FIGURE 2
# PERCENTAGE OF EACH 1-KB WINDOW COVERED
# ============================================================

p_windows <- ggplot(
    window_plot_data,
    aes(
        x = midpoint,
        y = percent,
        fill = track
    )
) +

    geom_col(
        width = window_size,
        linewidth = 0
    ) +

    facet_wrap(
        ~track,
        ncol = 1,
        scales = "fixed"
    ) +

    scale_x_continuous(
        limits = c(
            0,
            genome_length
        ),
        breaks = seq(
            0,
            100000,
            by = 10000
        ),
        labels = function(x) x / 1000,
        expand = c(
            0,
            0
        )
    ) +

    scale_y_continuous(
        limits = c(
            0,
            100
        ),
        breaks = seq(
            0,
            100,
            by = 20
        ),
        expand = c(
            0,
            0
        )
    ) +

    labs(
        title = "Distribution of low-complexity and tandem-repeat sequences across XwIV",
        x = "Position on XwIV genome (kb)",
        y = "Window covered (%)"
    ) +

    guides(
        fill = "none"
    ) +

    theme_classic(
        base_size = 12
    ) +

    theme(
        strip.background = element_blank(),

        strip.text = element_text(
            face = "bold",
            size = 11
        ),

        plot.title = element_text(
            face = "bold",
            hjust = 0.5
        ),

        panel.spacing = grid::unit(
            0.7,
            "lines"
        ),

        plot.margin = margin(
            15,
            20,
            15,
            20
        )
    )


ggsave(
    plot_windows_file,
    p_windows,
    width = 12,
    height = 6
)


# ============================================================
# SUMMARY
# ============================================================

cat("\n")
cat("============================================\n")
cat("SUMMARY\n")
cat("============================================\n")

cat("\nGenome length:", genome_length, "bp\n")

cat("\nDUST regions:", nrow(dust), "\n")

cat(
    "Unique DUST bp:",
    sum(
        dust_merged\$end -
        dust_merged\$start +
        1
    ),
    "\n"
)

cat(
    "Genome covered by DUST:",
    round(
        100 *
        sum(
            dust_merged\$end -
            dust_merged\$start +
            1
        ) /
        genome_length,
        3
    ),
    "%\n"
)


cat("\nTRF regions:", nrow(trf), "\n")

cat(
    "Unique TRF bp:",
    sum(
        trf_merged\$end -
        trf_merged\$start +
        1
    ),
    "\n"
)

cat(
    "Genome covered by TRF:",
    round(
        100 *
        sum(
            trf_merged\$end -
            trf_merged\$start +
            1
        ) /
        genome_length,
        3
    ),
    "%\n"
)


if (nrow(trf) > 0) {

    cat(
        "Longest individual TRF region:",
        max(trf\$length),
        "bp\n"
    )

    cat(
        "Largest TRF period:",
        max(trf\$period_size),
        "bp\n"
    )

    cat(
        "Maximum TRF copy number:",
        max(trf\$copy_number),
        "\n"
    )
}


cat("\nMaximum DUST coverage in a 1-kb window:")
cat(
    round(
        max(windows\$DUST_percent),
        2
    ),
    "%\n"
)

cat(
    "Window:",
    windows\$start[
        which.max(
            windows\$DUST_percent
        )
    ],
    "-",
    windows\$end[
        which.max(
            windows\$DUST_percent
        )
    ],
    "\n"
)


cat("\nMaximum TRF coverage in a 1-kb window:")
cat(
    round(
        max(windows\$TRF_percent),
        2
    ),
    "%\n"
)

cat(
    "Window:",
    windows\$start[
        which.max(
            windows\$TRF_percent
        )
    ],
    "-",
    windows\$end[
        which.max(
            windows\$TRF_percent
        )
    ],
    "\n"
)


cat("\nFigure 1:\n")
cat(plot_regions_file)
cat("\n")

cat("\nFigure 2:\n")
cat(plot_windows_file)
cat("\n")

cat("\nWindow table:\n")
cat(window_table_file)
cat("\n")


# ============================================================
# END R
# ============================================================

EOF


# ============================================================
# RUN R
# ============================================================

echo
echo "============================================"
echo "Generating figures"
echo "============================================"

Rscript "${R_SCRIPT}"


# ============================================================
# DONE
# ============================================================

echo
echo "============================================"
echo "DONE"
echo "============================================"

echo
echo "Converted TRF table:"
echo "${TRF_TSV}"

echo
echo "Window table:"
echo "${WINDOW_TABLE}"

echo
echo "Figure 1:"
echo "${PLOT_REGIONS}"

echo
echo "Figure 2:"
echo "${PLOT_WINDOWS}"

echo
