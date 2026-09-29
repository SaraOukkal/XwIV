#!/bin/bash

#SBATCH --job-name=XwIV_circularity
#SBATCH --output=/beegfs/data/soukkal/Thesis/Hymenoptera_Project/PDV/Results/Viral_genes/IVSPER_Analysis/Long_reads/07_Hybrid_assembly/Circularity/XwIV_circularity_%j.out
#SBATCH --error=/beegfs/data/soukkal/Thesis/Hymenoptera_Project/PDV/Results/Viral_genes/IVSPER_Analysis/Long_reads/07_Hybrid_assembly/Circularity/XwIV_circularity_%j.err
#SBATCH --cpus-per-task=8
#SBATCH --mem=16G
#SBATCH --time=02:00:00
#SBATCH --constraint="skylake|haswell|broadwell"
#SBATCH --exclude=pbil-deb27

set -euo pipefail


# ============================================================
# PARAMETERS
# ============================================================

CONTIG="jcf7180000126179"
ONT_CONTIG="contig_28997"

# Number of bp extracted from each extremity of XwIV
FLANK=10000


# ============================================================
# PATHS
# ============================================================

ASSEMBLY="/beegfs/data/soukkal/Thesis/Hymenoptera_Project/PDV/Results/Viral_genes/IVSPER_Analysis/Long_reads/07_Hybrid_assembly/MaSuRCA/CA.mr.99.17.15.0.02/primary.genome.scf.fasta"

ONT_ASSEMBLY="/beegfs/data/soukkal/Thesis/Hymenoptera_Project/PDV/Results/Viral_genes/IVSPER_Analysis/Long_reads/03_Long_read_assembly/Racon/Xiphosomella_wirra_ONT.racon2.fa"

ONT="/beegfs/data/soukkal/Thesis/Hymenoptera_Project/PDV/Results/Viral_genes/IVSPER_Analysis/Long_reads/01_Reads/Xiphosomella_wirra_ONT.fastq.gz"

OUTDIR="/beegfs/data/soukkal/Thesis/Hymenoptera_Project/PDV/Results/Viral_genes/IVSPER_Analysis/Long_reads/07_Hybrid_assembly/Circularity"

MINIMAP2="/beegfs/home/soukkal/miniconda3/bin/minimap2"
SAMTOOLS="/beegfs/data/soft/samtools-1.19/bin/samtools"

mkdir -p "${OUTDIR}"


# ============================================================
# OUTPUT FILES
# ============================================================

XWIV="${OUTDIR}/XwIV_hybrid.fa"

END_FA="${OUTDIR}/XwIV_END_${FLANK}bp.fa"
START_FA="${OUTDIR}/XwIV_START_${FLANK}bp.fa"

JUNCTION_FA="${OUTDIR}/XwIV_END_START_junction.fa"


# ONT READ OUTPUTS

READS_PAF="${OUTDIR}/ONT_reads_vs_XwIV_END_START.paf"
READS_TABLE="${OUTDIR}/ONT_reads_vs_XwIV_END_START.tsv"
READS_BAM="${OUTDIR}/ONT_reads_vs_XwIV_END_START.sorted.bam"

READS_R="${OUTDIR}/plot_XwIV_END_START_reads.R"
READS_PLOT="${OUTDIR}/XwIV_END_START_all_ONT_reads.pdf"


# ONT CONTIG OUTPUTS

ONT_CONTIG_FA="${OUTDIR}/${ONT_CONTIG}.fa"

CONTIG_PAF="${OUTDIR}/${ONT_CONTIG}_vs_XwIV_END_START.paf"
CONTIG_TABLE="${OUTDIR}/${ONT_CONTIG}_vs_XwIV_END_START.tsv"

CONTIG_R="${OUTDIR}/plot_${ONT_CONTIG}_vs_XwIV_END_START.R"
CONTIG_PLOT="${OUTDIR}/${ONT_CONTIG}_vs_XwIV_END_START.pdf"


# ============================================================
# CHECK INPUT FILES
# ============================================================

if [[ ! -s "${ASSEMBLY}" ]]; then
    echo "ERROR: Hybrid assembly not found:"
    echo "${ASSEMBLY}"
    exit 1
fi

if [[ ! -s "${ONT_ASSEMBLY}" ]]; then
    echo "ERROR: ONT-only assembly not found:"
    echo "${ONT_ASSEMBLY}"
    exit 1
fi

if [[ ! -s "${ONT}" ]]; then
    echo "ERROR: ONT reads not found:"
    echo "${ONT}"
    exit 1
fi


# ============================================================
# EXTRACT XwIV HYBRID SCAFFOLD
# ============================================================

echo
echo "============================================"
echo "Extracting hybrid XwIV scaffold"
echo "============================================"

"${SAMTOOLS}" faidx "${ASSEMBLY}"

"${SAMTOOLS}" faidx \
    "${ASSEMBLY}" \
    "${CONTIG}" \
    > "${XWIV}"

"${SAMTOOLS}" faidx "${XWIV}"

CONTIG_LENGTH=$(cut -f2 "${XWIV}.fai")

echo "Hybrid XwIV scaffold: ${CONTIG}"
echo "Length: ${CONTIG_LENGTH} bp"


# ============================================================
# EXTRACT TERMINAL REGIONS
# ============================================================

END_START=$((CONTIG_LENGTH - FLANK + 1))

echo
echo "XwIV END:"
echo "${END_START}-${CONTIG_LENGTH}"
echo
echo "XwIV START:"
echo "1-${FLANK}"


"${SAMTOOLS}" faidx \
    "${XWIV}" \
    "${CONTIG}:${END_START}-${CONTIG_LENGTH}" \
    > "${END_FA}"

"${SAMTOOLS}" faidx \
    "${XWIV}" \
    "${CONTIG}:1-${FLANK}" \
    > "${START_FA}"


# ============================================================
# CREATE ARTIFICIAL END -> START REFERENCE
#
# 1                    10000 10001                 20000
# |-----------------------| |-------------------------|
#        XwIV END                   XwIV START
#
# ============================================================

echo
echo "Creating artificial END -> START reference..."

{
    echo ">XwIV_END_START"

    grep -v "^>" "${END_FA}" | tr -d '\n'
    grep -v "^>" "${START_FA}" | tr -d '\n'

    echo
} > "${JUNCTION_FA}"

"${SAMTOOLS}" faidx "${JUNCTION_FA}"

echo "Artificial reference:"
echo "${JUNCTION_FA}"
echo
echo "Junction:"
echo "${FLANK} | $((FLANK + 1))"


# ============================================================
# PART 1
# MAP ALL RAW ONT READS
#
# NO FILTERING.
# ============================================================

echo
echo "============================================"
echo "Mapping all raw ONT reads"
echo "============================================"

"${MINIMAP2}" \
    -x map-ont \
    -t 8 \
    "${JUNCTION_FA}" \
    "${ONT}" \
    > "${READS_PAF}"


# ============================================================
# BAM FOR IGV
# ============================================================

echo
echo "Generating ONT BAM for IGV..."

"${MINIMAP2}" \
    -ax map-ont \
    -t 8 \
    "${JUNCTION_FA}" \
    "${ONT}" \
    | "${SAMTOOLS}" sort \
        -@ 8 \
        -o "${READS_BAM}"

"${SAMTOOLS}" index "${READS_BAM}"


# ============================================================
# CONVERT READ PAF TO TABLE
#
# Every minimap2 alignment is retained.
# ============================================================

echo
echo "Creating ONT read alignment table..."

awk '
BEGIN {
    OFS="\t"

    print "read_ID",
          "read_length",
          "read_start",
          "read_end",
          "strand",
          "target_start",
          "target_end",
          "matches",
          "alignment_length",
          "identity_percent",
          "MAPQ"
}
{
    identity=100*$10/$11

    print $1,
          $2,
          $3,
          $4,
          $5,
          $8+1,
          $9,
          $10,
          $11,
          identity,
          $12
}' "${READS_PAF}" > "${READS_TABLE}"


# ============================================================
# PLOT ALL RAW ONT READ ALIGNMENTS
#
# Multiple alignments belonging to the same read are placed
# on the same horizontal line.
#
# NO FILTERING.
# ============================================================

cat > "${READS_R}" <<EOF

library(ggplot2)

input_file <- "${READS_TABLE}"
output_file <- "${READS_PLOT}"

junction <- ${FLANK}
reference_length <- ${FLANK} * 2

dat <- read.delim(
    input_file,
    header = TRUE,
    sep = "\t",
    stringsAsFactors = FALSE,
    check.names = FALSE
)

if (nrow(dat) == 0) {

    pdf(output_file, width = 12, height = 4)

    plot.new()

    text(
        0.5,
        0.5,
        "No ONT reads mapped to the XwIV terminal reference",
        cex = 1.2
    )

    dev.off()

    quit(save = "no")
}


# ------------------------------------------------------------
# ONE Y POSITION PER READ
# ------------------------------------------------------------

read_order <- unique(
    dat\$read_ID[
        order(
            dat\$target_start,
            dat\$target_end
        )
    ]
)

dat\$read_ID <- factor(
    dat\$read_ID,
    levels = read_order
)

dat\$read_y <- as.numeric(dat\$read_ID)


# ------------------------------------------------------------
# PLOT
# ------------------------------------------------------------

p <- ggplot(dat) +

    annotate(
        "rect",
        xmin = 1,
        xmax = junction,
        ymin = -Inf,
        ymax = Inf,
        alpha = 0.08
    ) +

    annotate(
        "rect",
        xmin = junction,
        xmax = reference_length,
        ymin = -Inf,
        ymax = Inf,
        alpha = 0.03
    ) +

    geom_segment(
        aes(
            x = target_start,
            xend = target_end,
            y = read_y,
            yend = read_y
        ),
        linewidth = 0.7,
        lineend = "round"
    ) +

    geom_vline(
        xintercept = junction,
        linetype = "dashed",
        linewidth = 0.8
    ) +

    annotate(
        "text",
        x = junction / 2,
        y = length(read_order) + 2,
        label = "XwIV END",
        fontface = "bold",
        size = 4
    ) +

    annotate(
        "text",
        x = junction + junction / 2,
        y = length(read_order) + 2,
        label = "XwIV START",
        fontface = "bold",
        size = 4
    ) +

    scale_x_continuous(
        limits = c(1, reference_length),
        breaks = seq(
            0,
            reference_length,
            by = 2500
        ),
        expand = c(0, 0)
    ) +

    scale_y_continuous(
        breaks = NULL,
        expand = expansion(
            mult = c(0.01, 0.05)
        )
    ) +

    labs(
        title = "ONT reads mapped to the XwIV END-to-START reference",
        subtitle = paste0(
            length(read_order),
            " mapped reads; ",
            nrow(dat),
            " alignments; no filtering"
        ),
        x = "XwIV END (10 kb) | XwIV START (10 kb)",
        y = "ONT reads"
    ) +

    theme_classic(base_size = 12) +

    theme(
        axis.line.y = element_blank(),
        axis.ticks.y = element_blank()
    )


plot_height <- max(
    5,
    min(
        14,
        3 + length(read_order) * 0.03
    )
)

ggsave(
    output_file,
    p,
    width = 12,
    height = plot_height
)

cat(
    "\nMapped reads:",
    length(read_order),
    "\n"
)

cat(
    "Alignments:",
    nrow(dat),
    "\n"
)

cat(
    "Plot:",
    output_file,
    "\n"
)

EOF

echo
echo "Generating raw ONT read plot..."

Rscript "${READS_R}"


# ============================================================
# PART 2
# EXTRACT contig_28997 FROM THE ONT-ONLY ASSEMBLY
# ============================================================

echo
echo "============================================"
echo "Extracting ${ONT_CONTIG}"
echo "============================================"

"${SAMTOOLS}" faidx "${ONT_ASSEMBLY}"

"${SAMTOOLS}" faidx \
    "${ONT_ASSEMBLY}" \
    "${ONT_CONTIG}" \
    > "${ONT_CONTIG_FA}"

"${SAMTOOLS}" faidx "${ONT_CONTIG_FA}"

ONT_CONTIG_LENGTH=$(cut -f2 "${ONT_CONTIG_FA}.fai")

echo "${ONT_CONTIG} length: ${ONT_CONTIG_LENGTH} bp"


# ============================================================
# MAP contig_28997 AGAINST THE SAME ARTIFICIAL JUNCTION
#
# asm5 is appropriate here because we are aligning assembled
# nucleotide sequences expected to be highly similar.
#
# NO alignment filtering is performed afterwards.
# ============================================================

echo
echo "============================================"
echo "Mapping ${ONT_CONTIG} against XwIV END -> START"
echo "============================================"

"${MINIMAP2}" \
    -x asm5 \
    -c \
    "${JUNCTION_FA}" \
    "${ONT_CONTIG_FA}" \
    > "${CONTIG_PAF}"


# ============================================================
# CONVERT CONTIG PAF TO TABLE
# ============================================================

awk '
BEGIN {
    OFS="\t"

    print "query",
          "query_length",
          "query_start",
          "query_end",
          "strand",
          "target",
          "target_length",
          "target_start",
          "target_end",
          "matches",
          "alignment_length",
          "identity_percent",
          "MAPQ"
}
{
    identity=100*$10/$11

    print $1,
          $2,
          $3+1,
          $4,
          $5,
          $6,
          $7,
          $8+1,
          $9,
          $10,
          $11,
          identity,
          $12
}' "${CONTIG_PAF}" > "${CONTIG_TABLE}"


# ============================================================
# PLOT contig_28997 AGAINST ARTIFICIAL JUNCTION
#
# All minimap2 alignments are displayed.
# ============================================================

cat > "${CONTIG_R}" <<EOF

library(ggplot2)

input_file <- "${CONTIG_TABLE}"
output_file <- "${CONTIG_PLOT}"

junction <- ${FLANK}
reference_length <- ${FLANK} * 2

dat <- read.delim(
    input_file,
    header = TRUE,
    sep = "\t",
    stringsAsFactors = FALSE,
    check.names = FALSE
)

if (nrow(dat) == 0) {

    pdf(output_file, width = 12, height = 4)

    plot.new()

    text(
        0.5,
        0.5,
        "${ONT_CONTIG} did not map to the artificial XwIV junction",
        cex = 1.2
    )

    dev.off()

    quit(save = "no")
}


# ------------------------------------------------------------
# PLOT
# ------------------------------------------------------------

p <- ggplot(dat) +

    annotate(
        "rect",
        xmin = 1,
        xmax = junction,
        ymin = -Inf,
        ymax = Inf,
        alpha = 0.08
    ) +

    annotate(
        "rect",
        xmin = junction,
        xmax = reference_length,
        ymin = -Inf,
        ymax = Inf,
        alpha = 0.03
    ) +

    # Reference backbone
    geom_segment(
        aes(
            x = 1,
            xend = reference_length,
            y = 0,
            yend = 0
        ),
        linewidth = 0.5
    ) +

    # All alignments of contig_28997
    geom_segment(
        aes(
            x = target_start,
            xend = target_end,
            y = seq_len(nrow(dat)),
            yend = seq_len(nrow(dat))
        ),
        linewidth = 2,
        lineend = "round"
    ) +

    # Junction
    geom_vline(
        xintercept = junction,
        linetype = "dashed",
        linewidth = 0.8
    ) +

    annotate(
        "text",
        x = junction / 2,
        y = nrow(dat) + 1,
        label = "XwIV END",
        fontface = "bold",
        size = 4
    ) +

    annotate(
        "text",
        x = junction + junction / 2,
        y = nrow(dat) + 1,
        label = "XwIV START",
        fontface = "bold",
        size = 4
    ) +

    scale_x_continuous(
        limits = c(1, reference_length),
        breaks = seq(
            0,
            reference_length,
            by = 2500
        ),
        expand = c(0, 0)
    ) +

    scale_y_continuous(
        breaks = seq_len(nrow(dat)),
        labels = paste0(
            "alignment ",
            seq_len(nrow(dat))
        ),
        expand = expansion(
            mult = c(0.10, 0.15)
        )
    ) +

    labs(
        title = "${ONT_CONTIG} mapped to the XwIV END-to-START reference",
        subtitle = paste0(
            nrow(dat),
            " minimap2 alignment(s); no filtering"
        ),
        x = "XwIV END (10 kb) | XwIV START (10 kb)",
        y = "${ONT_CONTIG}"
    ) +

    theme_classic(base_size = 12)


ggsave(
    output_file,
    p,
    width = 12,
    height = max(4, 3 + nrow(dat) * 0.5)
)

cat(
    "\n${ONT_CONTIG} alignments:",
    nrow(dat),
    "\n"
)

cat(
    "Plot:",
    output_file,
    "\n"
)

EOF


echo
echo "Generating ${ONT_CONTIG} plot..."

Rscript "${CONTIG_R}"


# ============================================================
# SUMMARY
# ============================================================

N_READ_ALIGNMENTS=$(awk 'END{print NR-1}' "${READS_TABLE}")

N_READS=$(awk '
NR > 1 {
    reads[$1]=1
}
END {
    print length(reads)
}' "${READS_TABLE}")

N_CONTIG_ALIGNMENTS=$(awk 'END{print NR-1}' "${CONTIG_TABLE}")


echo
echo "============================================"
echo "XwIV TERMINAL MAPPING COMPLETE"
echo "============================================"
echo
echo "Hybrid XwIV:"
echo "  ${CONTIG}"
echo "  Length = ${CONTIG_LENGTH} bp"
echo
echo "Artificial reference:"
echo "  END   = ${END_START}-${CONTIG_LENGTH}"
echo "  START = 1-${FLANK}"
echo "  Junction = ${FLANK} | $((FLANK + 1))"
echo
echo "--------------------------------------------"
echo "RAW ONT READS"
echo "--------------------------------------------"
echo
echo "Mapped reads:"
echo "${N_READS}"
echo
echo "Alignments:"
echo "${N_READ_ALIGNMENTS}"
echo
echo "Table:"
echo "${READS_TABLE}"
echo
echo "BAM:"
echo "${READS_BAM}"
echo
echo "Plot:"
echo "${READS_PLOT}"
echo
echo "--------------------------------------------"
echo "${ONT_CONTIG}"
echo "--------------------------------------------"
echo
echo "Length:"
echo "${ONT_CONTIG_LENGTH} bp"
echo
echo "Alignments:"
echo "${N_CONTIG_ALIGNMENTS}"
echo
echo "Table:"
echo "${CONTIG_TABLE}"
echo
echo "Plot:"
echo "${CONTIG_PLOT}"
echo
echo "============================================"
