#!/bin/bash

#SBATCH --job-name=XwIV_compmap
#SBATCH --cpus-per-task=16
#SBATCH --mem=64G
#SBATCH --time=24:00:00
#SBATCH --output=XwIV_competitive_mapping.%j.out
#SBATCH --error=XwIV_competitive_mapping.%j.err

set -euo pipefail


# ============================================================
# CONFIGURATION
# ============================================================

THREADS=16

ASSEMBLY="/beegfs/data/soukkal/Thesis/Hymenoptera_Project/PDV/Results/Viral_genes/IVSPER_Analysis/Long_reads/07_Hybrid_assembly/MaSuRCA/CA.mr.99.17.15.0.02/primary.genome.scf.fasta"

VIRAL_CONTIG="jcf7180000126179"

ONT="/beegfs/data/soukkal/Thesis/Hymenoptera_Project/PDV/Results/Viral_genes/IVSPER_Analysis/Long_reads/01_Reads/Xiphosomella_wirra_ONT.fastq.gz"

ILLUMINA_R1="/beegfs/data/soukkal/Thesis/Hymenoptera_Project/PDV/Results/Viral_genes/IVSPER_Analysis/Long_reads/07_Hybrid_assembly/Reads/Xwirra_combined_R1.fastq.gz"
ILLUMINA_R2="/beegfs/data/soukkal/Thesis/Hymenoptera_Project/PDV/Results/Viral_genes/IVSPER_Analysis/Long_reads/07_Hybrid_assembly/Reads/Xwirra_combined_R2.fastq.gz"

MINIMAP2="/beegfs/home/soukkal/miniconda3/bin/minimap2"
SAMTOOLS="/beegfs/data/soft/samtools-1.19/bin/samtools"

BASE="/beegfs/data/soukkal/Thesis/Hymenoptera_Project/PDV/Results/Viral_genes/IVSPER_Analysis/Long_reads/07_Hybrid_assembly/Rotated_junction_mapping"

OUTDIR="${BASE}/Competitive_mapping"

REFDIR="${OUTDIR}/Reference"
ONTDIR="${OUTDIR}/ONT"
ILLUMINADIR="${OUTDIR}/Illumina"
PLOTDIR="${OUTDIR}/Plots"

mkdir -p "${REFDIR}" "${ONTDIR}" "${ILLUMINADIR}" "${PLOTDIR}"

ROTATED_ASSEMBLY="${REFDIR}/Xwirra_hybrid_XwIV_rotated.fa"
ROTATION_INFO="${REFDIR}/XwIV_rotation_info.tsv"

ONT_BAM="${ONTDIR}/Xwirra_ONT_vs_complete_hybrid_rotated.sorted.bam"
ILLUMINA_BAM="${ILLUMINADIR}/Xwirra_Illumina_vs_complete_hybrid_rotated.sorted.bam"


# ============================================================
# CHECK INPUT
# ============================================================

for FILE in "${ASSEMBLY}" "${ONT}" "${ILLUMINA_R1}" "${ILLUMINA_R2}"
do
    if [[ ! -s "${FILE}" ]]; then
        echo "ERROR: missing or empty file: ${FILE}"
        exit 1
    fi
done


# ============================================================
# INDEX ORIGINAL ASSEMBLY
# ============================================================

"${SAMTOOLS}" faidx "${ASSEMBLY}"

"${SAMTOOLS}" faidx "${ASSEMBLY}" "${VIRAL_CONTIG}" > "${REFDIR}/XwIV_original.fa"

VIRAL_LENGTH=$(awk '!/^>/ {gsub(/[[:space:]]/,""); n+=length($0)} END {print n}' "${REFDIR}/XwIV_original.fa")

echo "XwIV length: ${VIRAL_LENGTH} bp"


# ============================================================
# ROTATE XwIV
# ============================================================

LEFT_HALF_END=$(( VIRAL_LENGTH / 2 ))
RIGHT_HALF_START=$(( LEFT_HALF_END + 1 ))

FIRST_PIECE_LENGTH=$(( VIRAL_LENGTH - RIGHT_HALF_START + 1 ))

JUNCTION_POSITION="${FIRST_PIECE_LENGTH}"

echo "Rotation point: ${LEFT_HALF_END}|${RIGHT_HALF_START}"
echo "END -> START junction after rotated position: ${JUNCTION_POSITION}"


"${SAMTOOLS}" faidx "${ASSEMBLY}" "${VIRAL_CONTIG}:${RIGHT_HALF_START}-${VIRAL_LENGTH}" | grep -v "^>" | tr -d '\n' > "${REFDIR}/XwIV_part2.tmp"

"${SAMTOOLS}" faidx "${ASSEMBLY}" "${VIRAL_CONTIG}:1-${LEFT_HALF_END}" | grep -v "^>" | tr -d '\n' > "${REFDIR}/XwIV_part1.tmp"

{
    echo ">${VIRAL_CONTIG}"
    cat "${REFDIR}/XwIV_part2.tmp"
    cat "${REFDIR}/XwIV_part1.tmp"
    echo
} > "${REFDIR}/XwIV_rotated.fa"

rm -f "${REFDIR}/XwIV_part1.tmp" "${REFDIR}/XwIV_part2.tmp"


# ============================================================
# CREATE COMPLETE COMPETITIVE REFERENCE
# ============================================================

python - "${ASSEMBLY}" "${REFDIR}/XwIV_rotated.fa" "${VIRAL_CONTIG}" "${ROTATED_ASSEMBLY}" <<'PYTHON'

import sys

assembly_file, rotated_file, target, output_file = sys.argv[1:]


def read_fasta(filename):

    name = None
    seq = []

    with open(filename) as handle:

        for line in handle:

            line = line.rstrip()

            if not line:
                continue

            if line.startswith(">"):

                if name is not None:
                    yield name, "".join(seq)

                name = line[1:].split()[0]
                seq = []

            else:

                seq.append(line)

        if name is not None:
            yield name, "".join(seq)


rotated = list(
    read_fasta(rotated_file)
)

if len(rotated) != 1:
    raise RuntimeError(
        "Rotated FASTA must contain exactly one sequence."
    )

rotated_name, rotated_seq = rotated[0]

if rotated_name != target:
    raise RuntimeError(
        "Rotated sequence name does not match target contig."
    )

found = 0

with open(output_file, "w") as out:

    for name, sequence in read_fasta(assembly_file):

        if name == target:

            sequence = rotated_seq
            found += 1

        out.write(
            f">{name}\n"
        )

        for i in range(
            0,
            len(sequence),
            80
        ):

            out.write(
                sequence[i:i+80]
                +
                "\n"
            )

if found != 1:

    raise RuntimeError(
        f"Expected one {target}; found {found}."
    )

print(
    f"Replaced {target} with rotated sequence."
)

PYTHON


# ============================================================
# SAVE ROTATION INFORMATION
# ============================================================

{
    printf "parameter\tvalue\n"
    printf "contig\t%s\n" "${VIRAL_CONTIG}"
    printf "genome_length\t%s\n" "${VIRAL_LENGTH}"
    printf "left_half_end\t%s\n" "${LEFT_HALF_END}"
    printf "right_half_start\t%s\n" "${RIGHT_HALF_START}"
    printf "junction_after_position\t%s\n" "${JUNCTION_POSITION}"
} > "${ROTATION_INFO}"


# ============================================================
# INDEX COMPETITIVE REFERENCE
# ============================================================

"${SAMTOOLS}" faidx "${ROTATED_ASSEMBLY}"


# ============================================================
# ONT COMPETITIVE MAPPING
# ============================================================

echo "Mapping ONT reads..."

"${MINIMAP2}" -ax map-ont -t "${THREADS}" "${ROTATED_ASSEMBLY}" "${ONT}" | "${SAMTOOLS}" sort -@ "${THREADS}" -o "${ONT_BAM}" -

"${SAMTOOLS}" index -@ "${THREADS}" "${ONT_BAM}"


# ============================================================
# ILLUMINA COMPETITIVE MAPPING
# ============================================================

echo "Mapping Illumina reads..."

"${MINIMAP2}" -ax sr -t "${THREADS}" "${ROTATED_ASSEMBLY}" "${ILLUMINA_R1}" "${ILLUMINA_R2}" | "${SAMTOOLS}" sort -@ "${THREADS}" -o "${ILLUMINA_BAM}" -

"${SAMTOOLS}" index -@ "${THREADS}" "${ILLUMINA_BAM}"


# ============================================================
# FLAGSTAT
# ============================================================

"${SAMTOOLS}" flagstat "${ONT_BAM}" > "${ONTDIR}/Xwirra_ONT.flagstat.txt"

"${SAMTOOLS}" flagstat "${ILLUMINA_BAM}" > "${ILLUMINADIR}/Xwirra_Illumina.flagstat.txt"


# ============================================================
# PRIMARY ALIGNMENTS ON XwIV
#
# -F 2308 removes:
#   4    unmapped
#   256  secondary
#   2048 supplementary
# ============================================================

"${SAMTOOLS}" view -bh -F 2308 "${ONT_BAM}" "${VIRAL_CONTIG}" > "${ONTDIR}/XwIV_primary.bam"

"${SAMTOOLS}" index "${ONTDIR}/XwIV_primary.bam"


"${SAMTOOLS}" view -bh -F 2308 "${ILLUMINA_BAM}" "${VIRAL_CONTIG}" > "${ILLUMINADIR}/XwIV_primary.bam"

"${SAMTOOLS}" index "${ILLUMINADIR}/XwIV_primary.bam"


# ============================================================
# DEPTH
#
# Restrict explicitly to XwIV.
# ============================================================

"${SAMTOOLS}" depth -aa -r "${VIRAL_CONTIG}" "${ONTDIR}/XwIV_primary.bam" > "${ONTDIR}/XwIV_primary.depth.tsv"

"${SAMTOOLS}" depth -aa -r "${VIRAL_CONTIG}" "${ILLUMINADIR}/XwIV_primary.bam" > "${ILLUMINADIR}/XwIV_primary.depth.tsv"


# ============================================================
# PLOTS AND PAIRED-END ANALYSIS
# ============================================================

python - "${ONTDIR}/XwIV_primary.bam" "${ILLUMINADIR}/XwIV_primary.bam" "${ONTDIR}/XwIV_primary.depth.tsv" "${ILLUMINADIR}/XwIV_primary.depth.tsv" "${SAMTOOLS}" "${VIRAL_CONTIG}" "${JUNCTION_POSITION}" "${VIRAL_LENGTH}" "${PLOTDIR}" <<'PYTHON'

import sys
import os
import re
import subprocess
import statistics

import numpy as np
import matplotlib.pyplot as plt


# ============================================================
# INPUT
# ============================================================

ONT_BAM = sys.argv[1]
ILLUMINA_BAM = sys.argv[2]

ONT_DEPTH_FILE = sys.argv[3]
ILLUMINA_DEPTH_FILE = sys.argv[4]

SAMTOOLS = sys.argv[5]
CONTIG = sys.argv[6]

JUNCTION = int(
    sys.argv[7]
)

GENOME_LENGTH = int(
    sys.argv[8]
)

OUTDIR = sys.argv[9]

ZOOM_FLANK = 500

REGION_START = max(
    1,
    JUNCTION - ZOOM_FLANK
)

REGION_END = min(
    GENOME_LENGTH,
    JUNCTION + ZOOM_FLANK
)


# ============================================================
# READ DEPTH
# ============================================================

def read_depth(filename):

    positions = []
    depths = []

    with open(filename) as handle:

        for line in handle:

            fields = line.rstrip().split(
                "\t"
            )

            if len(fields) < 3:
                continue

            # Only analyse XwIV
            if fields[0] != CONTIG:
                continue

            positions.append(
                int(fields[1])
            )

            depths.append(
                int(fields[2])
            )

    return (
        np.asarray(positions),
        np.asarray(depths)
    )


ont_pos, ont_depth = read_depth(
    ONT_DEPTH_FILE
)

ill_pos, ill_depth = read_depth(
    ILLUMINA_DEPTH_FILE
)


# ============================================================
# DEPTH STATISTICS
# ============================================================

ont_median = np.median(
    ont_depth
)

ill_median = np.median(
    ill_depth
)

ont_max = np.max(
    ont_depth
)

ill_max = np.max(
    ill_depth
)


print()

print(
    "DEPTH"
)

print(
    f"ONT median: {ont_median:.1f}x"
)

print(
    f"ONT maximum: {ont_max}x"
)

print(
    f"Illumina median: {ill_median:.1f}x"
)

print(
    f"Illumina maximum: {ill_max}x"
)


# ============================================================
# FIGURE 1
# WHOLE-GENOME DEPTH
# ============================================================

fig, axes = plt.subplots(
    2,
    1,
    figsize=(
        16,
        8
    ),
    sharex=True
)


# ------------------------------------------------------------
# ONT
# ------------------------------------------------------------

axes[0].plot(
    ont_pos,
    ont_depth,
    color="#2166AC",
    linewidth=0.65
)

axes[0].axvline(
    JUNCTION + 0.5,
    color="black",
    linestyle="--",
    linewidth=1.5
)

axes[0].axhline(
    ont_median,
    color="grey",
    linestyle=":",
    linewidth=1.2,
    label=f"Median = {ont_median:.1f}×"
)

axes[0].set_ylabel(
    "ONT depth"
)

axes[0].set_title(
    "ONT long reads",
    fontweight="bold"
)

axes[0].legend(
    frameon=False
)


# ------------------------------------------------------------
# ILLUMINA
# ------------------------------------------------------------

axes[1].plot(
    ill_pos,
    ill_depth,
    color="#D73027",
    linewidth=0.65
)

axes[1].axvline(
    JUNCTION + 0.5,
    color="black",
    linestyle="--",
    linewidth=1.5
)

axes[1].axhline(
    ill_median,
    color="grey",
    linestyle=":",
    linewidth=1.2,
    label=f"Median = {ill_median:.1f}×"
)

axes[1].set_ylabel(
    "Illumina depth"
)

axes[1].set_xlabel(
    "Position on rotated XwIV genome (bp)"
)

axes[1].set_title(
    "Illumina short reads",
    fontweight="bold"
)

axes[1].legend(
    frameon=False
)


for ax in axes:

    ax.set_xlim(
        1,
        GENOME_LENGTH
    )

    ax.grid(
        axis="x",
        alpha=0.15
    )


axes[0].text(
    JUNCTION + 0.5,
    1.02,
    "END → START",
    transform=axes[0].get_xaxis_transform(),
    ha="center",
    va="bottom",
    fontsize=9
)


fig.suptitle(
    "Read depth across the rotated XwIV genome",
    fontsize=15,
    fontweight="bold"
)


plt.tight_layout(
    rect=[
        0,
        0,
        1,
        0.96
    ]
)


for ext in [
    "pdf",
    "png",
    "svg"
]:

    kwargs = {
        "bbox_inches": "tight"
    }

    if ext == "png":

        kwargs[
            "dpi"
        ] = 300

    plt.savefig(
        os.path.join(
            OUTDIR,
            f"XwIV_competitive_mapping_whole_genome_depth.{ext}"
        ),
        **kwargs
    )


plt.close()


# ============================================================
# CIGAR PARSING
# ============================================================

CIGAR_RE = re.compile(
    r"(\d+)([MIDNSHP=X])"
)


def reference_end(
    start,
    cigar
):

    ref_length = 0

    for length, operation in CIGAR_RE.findall(
        cigar
    ):

        if operation in {
            "M",
            "D",
            "N",
            "=",
            "X"
        }:

            ref_length += int(
                length
            )

    return (
        start
        +
        ref_length
        -
        1
    )


# ============================================================
# PARSE SAM ALIGNMENT
# ============================================================

def parse_sam_line(line):

    fields = line.split(
        "\t"
    )

    if len(fields) < 11:
        return None

    flag = int(
        fields[1]
    )

    start = int(
        fields[3]
    )

    cigar = fields[5]

    if cigar == "*":
        return None

    return {

        "qname":
            fields[0],

        "flag":
            flag,

        "rname":
            fields[2],

        "start":
            start,

        "end":
            reference_end(
                start,
                cigar
            ),

        "mapq":
            int(fields[4]),

        "cigar":
            cigar,

        "rnext":
            fields[6],

        "pnext":
            int(fields[7]),

        "tlen":
            int(fields[8]),

        "read1":
            bool(
                flag & 64
            ),

        "read2":
            bool(
                flag & 128
            ),

        "reverse":
            bool(
                flag & 16
            ),

        "proper_pair":
            bool(
                flag & 2
            )
    }


# ============================================================
# GET PRIMARY ALIGNMENTS
#
# -F 2308 removes:
# unmapped + secondary + supplementary
# ============================================================

def get_primary_alignments(
    bam,
    region=None
):

    command = [
        SAMTOOLS,
        "view",
        "-F",
        "2308",
        bam
    ]

    if region is not None:

        command.append(
            region
        )

    result = subprocess.run(
        command,
        capture_output=True,
        text=True,
        check=True
    )

    output = []

    for line in result.stdout.splitlines():

        aln = parse_sam_line(
            line
        )

        if aln is not None:

            output.append(
                aln
            )

    return output


REGION = (
    f"{CONTIG}:"
    f"{REGION_START}-"
    f"{REGION_END}"
)


# ============================================================
# ONT ALIGNMENTS
# ============================================================

ont_alignments = get_primary_alignments(
    ONT_BAM,
    REGION
)


for aln in ont_alignments:

    aln[
        "crosses_junction"
    ] = (
        aln["start"] <= JUNCTION
        and
        aln["end"] >= JUNCTION + 1
    )


ont_crossing = sum(
    aln[
        "crosses_junction"
    ]
    for aln in ont_alignments
)


# ============================================================
# STACK ONT ALIGNMENTS
# ============================================================

def stack_alignments(
    alignments
):

    ordered = sorted(
        alignments,
        key=lambda x: (
            x["start"],
            x["end"]
        )
    )

    levels_end = []
    output = []

    for aln in ordered:

        level = None

        for i, previous_end in enumerate(
            levels_end
        ):

            if aln["start"] > previous_end:

                level = i

                levels_end[
                    i
                ] = aln["end"]

                break

        if level is None:

            level = len(
                levels_end
            )

            levels_end.append(
                aln["end"]
            )

        x = aln.copy()

        x[
            "level"
        ] = level

        output.append(
            x
        )

    return (
        output,
        len(levels_end)
    )


ont_plot, ont_levels = stack_alignments(
    ont_alignments
)


# ============================================================
# ILLUMINA:
# IDENTIFY READ NAMES PRESENT IN THE ZOOM REGION
# ============================================================

illumina_region_alignments = get_primary_alignments(
    ILLUMINA_BAM,
    REGION
)


qnames_in_region = {
    aln["qname"]
    for aln in illumina_region_alignments
}


# ============================================================
# RETRIEVE BOTH PRIMARY MATES FROM THE COMPLETE XwIV CONTIG
# ============================================================

all_xwiv_alignments = get_primary_alignments(
    ILLUMINA_BAM,
    CONTIG
)


all_by_qname = {}


for aln in all_xwiv_alignments:

    if aln[
        "qname"
    ] not in qnames_in_region:

        continue

    all_by_qname.setdefault(
        aln["qname"],
        []
    ).append(
        aln
    )


# ============================================================
# BUILD COMPLETE ILLUMINA PAIRS
# ============================================================

pairs = []


for qname, reads in all_by_qname.items():

    r1_candidates = [
        r
        for r in reads
        if r["read1"]
    ]

    r2_candidates = [
        r
        for r in reads
        if r["read2"]
    ]


    # Require exactly one primary R1
    # and exactly one primary R2 on XwIV.

    if (
        len(r1_candidates) != 1
        or
        len(r2_candidates) != 1
    ):

        continue


    r1 = r1_candidates[0]
    r2 = r2_candidates[0]


    fragment_start = min(
        r1["start"],
        r2["start"]
    )

    fragment_end = max(
        r1["end"],
        r2["end"]
    )


    # Only retain fragments overlapping
    # the displayed ±500-bp region.

    if (
        fragment_end < REGION_START
        or
        fragment_start > REGION_END
    ):

        continue


    # --------------------------------------------------------
    # Determine genomic left/right mate
    # --------------------------------------------------------

    if r1["start"] <= r2["start"]:

        left = r1
        right = r2

    else:

        left = r2
        right = r1


    # --------------------------------------------------------
    # Does an actual sequenced read cross the junction?
    # --------------------------------------------------------

    r1_crosses = (
        r1["start"] <= JUNCTION
        and
        r1["end"] >= JUNCTION + 1
    )

    r2_crosses = (
        r2["start"] <= JUNCTION
        and
        r2["end"] >= JUNCTION + 1
    )

    read_crosses = (
        r1_crosses
        or
        r2_crosses
    )


    # --------------------------------------------------------
    # Does the inferred paired-end fragment span junction?
    # --------------------------------------------------------

    spans_junction = (
        fragment_start <= JUNCTION
        and
        fragment_end >= JUNCTION + 1
    )


    # --------------------------------------------------------
    # Expected inward-facing paired-end orientation:
    #
    #       --->          <---
    # --------------------------------------------------------

    inward_facing = (
        not left[
            "reverse"
        ]
        and
        right[
            "reverse"
        ]
    )


    # --------------------------------------------------------
    # Proper-pair flag must be present for both mates
    # --------------------------------------------------------

    proper_pair = (
        r1[
            "proper_pair"
        ]
        and
        r2[
            "proper_pair"
        ]
    )


    # --------------------------------------------------------
    # Concordant pair
    # --------------------------------------------------------

    concordant = (
        proper_pair
        and
        inward_facing
    )


    # --------------------------------------------------------
    # Template length
    # --------------------------------------------------------

    insert_size = max(
        abs(
            r1["tlen"]
        ),
        abs(
            r2["tlen"]
        )
    )


    pairs.append(
        {

            "qname":
                qname,

            "r1":
                r1,

            "r2":
                r2,

            "left":
                left,

            "right":
                right,

            "fragment_start":
                fragment_start,

            "fragment_end":
                fragment_end,

            "insert_size":
                insert_size,

            "proper_pair":
                proper_pair,

            "inward_facing":
                inward_facing,

            "concordant":
                concordant,

            "spans_junction":
                spans_junction,

            "r1_crosses":
                r1_crosses,

            "r2_crosses":
                r2_crosses,

            "read_crosses":
                read_crosses,

            "min_mapq":
                min(
                    r1["mapq"],
                    r2["mapq"]
                )
        }
    )


# ============================================================
# ILLUMINA PAIR STATISTICS
# ============================================================

spanning_pairs = [
    pair
    for pair in pairs
    if pair[
        "spans_junction"
    ]
]


concordant_spanning_pairs = [
    pair
    for pair in spanning_pairs
    if pair[
        "concordant"
    ]
]


illumina_read_crossing_pairs = [
    pair
    for pair in pairs
    if pair[
        "read_crosses"
    ]
]


print()

print(
    "============================================================"
)

print(
    "ILLUMINA PAIR ANALYSIS"
)

print(
    "============================================================"
)

print()

print(
    f"Junction: "
    f"{JUNCTION}|{JUNCTION + 1}"
)

print(
    f"Plot region: "
    f"{REGION_START}-{REGION_END}"
)

print()

print(
    f"Complete primary Illumina pairs in plotting region: "
    f"{len(pairs)}"
)

print()

print(
    f"Pairs whose inferred fragment spans junction: "
    f"{len(spanning_pairs)}"
)

print(
    f"  Proper pairs: "
    f"{sum(p['proper_pair'] for p in spanning_pairs)}"
)

print(
    f"  Inward-facing pairs: "
    f"{sum(p['inward_facing'] for p in spanning_pairs)}"
)

print(
    f"  Proper + inward-facing pairs: "
    f"{len(concordant_spanning_pairs)}"
)

print(
    f"Pairs with an actual read crossing junction: "
    f"{len(illumina_read_crossing_pairs)}"
)


# ============================================================
# STATISTICS FOR CONCORDANT SPANNING PAIRS
# ============================================================

if concordant_spanning_pairs:

    insert_sizes = [
        pair[
            "insert_size"
        ]
        for pair in concordant_spanning_pairs
    ]

    min_mapqs = [
        pair[
            "min_mapq"
        ]
        for pair in concordant_spanning_pairs
    ]

    print()

    print(
        "CONCORDANT SPANNING PAIRS"
    )

    print(
        f"Insert size: "
        f"min={min(insert_sizes)} bp, "
        f"median={statistics.median(insert_sizes):.1f} bp, "
        f"max={max(insert_sizes)} bp"
    )

    print(
        f"Minimum mate MAPQ: "
        f"min={min(min_mapqs)}, "
        f"median={statistics.median(min_mapqs):.1f}, "
        f"max={max(min_mapqs)}"
    )


# ============================================================
# SAVE ILLUMINA PAIR TABLES
# ============================================================

columns = [

    "read_name",

    "R1_start",
    "R1_end",
    "R1_strand",
    "R1_MAPQ",
    "R1_CIGAR",

    "R2_start",
    "R2_end",
    "R2_strand",
    "R2_MAPQ",
    "R2_CIGAR",

    "fragment_start",
    "fragment_end",

    "insert_size",

    "proper_pair",
    "inward_facing",
    "concordant",

    "spans_junction",

    "R1_crosses_junction",
    "R2_crosses_junction",

    "read_crosses_junction"
]


def pair_to_row(pair):

    r1 = pair["r1"]
    r2 = pair["r2"]

    return [

        pair[
            "qname"
        ],

        r1[
            "start"
        ],

        r1[
            "end"
        ],

        "-" if r1[
            "reverse"
        ] else "+",

        r1[
            "mapq"
        ],

        r1[
            "cigar"
        ],

        r2[
            "start"
        ],

        r2[
            "end"
        ],

        "-" if r2[
            "reverse"
        ] else "+",

        r2[
            "mapq"
        ],

        r2[
            "cigar"
        ],

        pair[
            "fragment_start"
        ],

        pair[
            "fragment_end"
        ],

        pair[
            "insert_size"
        ],

        pair[
            "proper_pair"
        ],

        pair[
            "inward_facing"
        ],

        pair[
            "concordant"
        ],

        pair[
            "spans_junction"
        ],

        pair[
            "r1_crosses"
        ],

        pair[
            "r2_crosses"
        ],

        pair[
            "read_crosses"
        ]
    ]


for filename, selected_pairs in [

    (
        "Illumina_pairs_junction_500bp.tsv",
        pairs
    ),

    (
        "Illumina_pairs_spanning_junction_detailed.tsv",
        spanning_pairs
    )
]:

    with open(
        os.path.join(
            OUTDIR,
            filename
        ),
        "w"
    ) as out:

        out.write(
            "\t".join(
                columns
            )
            +
            "\n"
        )

        for pair in selected_pairs:

            out.write(
                "\t".join(
                    str(x)
                    for x in pair_to_row(
                        pair
                    )
                )
                +
                "\n"
            )


# ============================================================
# STACK ILLUMINA PAIRS
#
# Each paired-end fragment is treated as one graphical object.
# ============================================================

def stack_pairs(
    pairs
):

    ordered = sorted(
        pairs,
        key=lambda x: (
            x[
                "fragment_start"
            ],
            x[
                "fragment_end"
            ]
        )
    )

    levels_end = []
    output = []


    for pair in ordered:

        level = None

        for i, previous_end in enumerate(
            levels_end
        ):

            if (
                pair[
                    "fragment_start"
                ]
                >
                previous_end + 5
            ):

                level = i

                levels_end[
                    i
                ] = pair[
                    "fragment_end"
                ]

                break


        if level is None:

            level = len(
                levels_end
            )

            levels_end.append(
                pair[
                    "fragment_end"
                ]
            )


        x = pair.copy()

        x[
            "level"
        ] = level

        output.append(
            x
        )


    return (
        output,
        len(levels_end)
    )


ill_plot, ill_levels = stack_pairs(
    pairs
)


# ============================================================
# FIGURE 2
# JUNCTION ±500 BP
# ============================================================

fig, axes = plt.subplots(
    2,
    1,
    figsize=(
        14,
        9
    ),
    sharex=True
)


# ============================================================
# ONT PANEL
#
# Solid blue = actual aligned read sequence.
# ============================================================

for aln in ont_plot:

    x1 = max(
        aln[
            "start"
        ],
        REGION_START
    )

    x2 = min(
        aln[
            "end"
        ],
        REGION_END
    )

    if x2 < x1:
        continue


    if aln[
        "crosses_junction"
    ]:

        read_color = "black"
        linewidth = 2.4

    else:

        read_color = "#2166AC"
        linewidth = 1.2


    axes[0].plot(

        [
            x1,
            x2
        ],

        [
            aln[
                "level"
            ],
            aln[
                "level"
            ]
        ],

        color=read_color,

        linewidth=linewidth,

        solid_capstyle="butt",

        zorder=3
    )


axes[0].axvline(
    JUNCTION + 0.5,
    color="black",
    linestyle="--",
    linewidth=1.8,
    zorder=10
)


axes[0].set_xlim(
    REGION_START,
    REGION_END
)


axes[0].set_ylim(
    -1,
    max(
        1,
        ont_levels
    )
)


axes[0].set_yticks([])


axes[0].set_ylabel(
    "ONT"
)


axes[0].set_title(
    (
        f"ONT — "
        f"{ont_crossing} primary read(s) crossing junction"
    ),
    fontsize=11,
    fontweight="bold"
)


# ============================================================
# ILLUMINA PANEL
#
# Red solid segments:
#     actual sequenced portions (R1 and R2).
#
# Grey dashed connector:
#     inferred unsequenced interval between mates.
#
# One pair = one graphical row.
# ============================================================

for pair in ill_plot:

    y = pair[
        "level"
    ]

    left = pair[
        "left"
    ]

    right = pair[
        "right"
    ]


    # --------------------------------------------------------
    # DASHED CONNECTOR BETWEEN MATES
    #
    # This interval is NOT directly sequenced.
    # --------------------------------------------------------

    connector_start = left[
        "end"
    ]

    connector_end = right[
        "start"
    ]


    if connector_end > connector_start:

        x1 = max(
            connector_start,
            REGION_START
        )

        x2 = min(
            connector_end,
            REGION_END
        )


        if x2 >= x1:

            axes[1].plot(

                [
                    x1,
                    x2
                ],

                [
                    y,
                    y
                ],

                color="0.70",

                linewidth=0.7,

                linestyle=(
                    0,
                    (
                        3,
                        3
                    )
                ),

                zorder=1
            )


    # --------------------------------------------------------
    # ACTUAL R1 / R2 ALIGNMENTS
    # --------------------------------------------------------

    for mate in [

        pair[
            "r1"
        ],

        pair[
            "r2"
        ]
    ]:

        x1 = max(
            mate[
                "start"
            ],
            REGION_START
        )

        x2 = min(
            mate[
                "end"
            ],
            REGION_END
        )


        if x2 < x1:
            continue


        mate_crosses = (
            mate[
                "start"
            ]
            <=
            JUNCTION
            and
            mate[
                "end"
            ]
            >=
            JUNCTION + 1
        )


        if mate_crosses:

            mate_color = "black"
            linewidth = 3.2

        else:

            mate_color = "#D73027"
            linewidth = 3.0


        axes[1].plot(

            [
                x1,
                x2
            ],

            [
                y,
                y
            ],

            color=mate_color,

            linewidth=linewidth,

            solid_capstyle="butt",

            zorder=3
        )


axes[1].axvline(
    JUNCTION + 0.5,
    color="black",
    linestyle="--",
    linewidth=1.8,
    zorder=10
)


axes[1].set_xlim(
    REGION_START,
    REGION_END
)


axes[1].set_ylim(
    -1,
    max(
        1,
        ill_levels
    )
)


axes[1].set_yticks([])


axes[1].set_ylabel(
    "Illumina"
)


axes[1].set_xlabel(
    "Position on rotated XwIV genome (bp)"
)


axes[1].set_title(

    (
        "Illumina paired-end — "
        f"{len(illumina_read_crossing_pairs)} "
        "read(s) crossing junction; "
        f"{len(concordant_spanning_pairs)} "
        "concordant pair(s) spanning junction"
    ),

    fontsize=11,

    fontweight="bold"
)


# ============================================================
# JUNCTION LABEL
# ============================================================

axes[0].text(

    JUNCTION + 0.5,

    1.02,

    "END → START",

    transform=axes[0].get_xaxis_transform(),

    ha="center",

    va="bottom",

    fontsize=10,

    fontweight="bold"
)


# ============================================================
# FINAL FORMATTING
# ============================================================

fig.suptitle(
    "Read alignments across the XwIV terminal junction",
    fontsize=15,
    fontweight="bold"
)


plt.tight_layout(
    rect=[
        0,
        0,
        1,
        0.95
    ]
)


# ============================================================
# SAVE FIGURE
# ============================================================

for ext in [
    "pdf",
    "png",
    "svg"
]:

    kwargs = {
        "bbox_inches":
            "tight"
    }

    if ext == "png":

        kwargs[
            "dpi"
        ] = 300


    plt.savefig(

        os.path.join(
            OUTDIR,
            f"XwIV_competitive_mapping_junction_500bp.{ext}"
        ),

        **kwargs
    )


plt.close()


# ============================================================
# FINISHED
# ============================================================

print()

print(
    "============================================================"
)

print(
    "PLOT OUTPUT"
)

print(
    "============================================================"
)

print()

print(
    "Junction figure:"
)

print(
    os.path.join(
        OUTDIR,
        "XwIV_competitive_mapping_junction_500bp.pdf"
    )
)

print()

print(
    "Illumina pairs:"
)

print(
    os.path.join(
        OUTDIR,
        "Illumina_pairs_junction_500bp.tsv"
    )
)

print()

print(
    "Illumina spanning pairs:"
)

print(
    os.path.join(
        OUTDIR,
        "Illumina_pairs_spanning_junction_detailed.tsv"
    )
)

PYTHON


# ============================================================
# FINISHED
# ============================================================

echo
echo "============================================================"
echo "COMPETITIVE MAPPING COMPLETE"
echo "============================================================"
echo

echo "Competitive reference:"
echo "${ROTATED_ASSEMBLY}"

echo
echo "ONT BAM:"
echo "${ONT_BAM}"

echo
echo "Illumina BAM:"
echo "${ILLUMINA_BAM}"

echo
echo "Whole-genome depth figure:"
echo "${PLOTDIR}/XwIV_competitive_mapping_whole_genome_depth.pdf"

echo
echo "Junction mapping figure:"
echo "${PLOTDIR}/XwIV_competitive_mapping_junction_500bp.pdf"

echo
echo "Illumina paired-end table:"
echo "${PLOTDIR}/Illumina_pairs_junction_500bp.tsv"

echo
echo "Illumina pairs spanning junction:"
echo "${PLOTDIR}/Illumina_pairs_spanning_junction_detailed.tsv"

echo
echo "Done."
