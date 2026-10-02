#!/usr/bin/env python3

import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt


# ============================================================
# CONFIGURATION
# ============================================================

BASE = (
    "/beegfs/data/soukkal/Thesis/Hymenoptera_Project/PDV/Results/"
    "Viral_genes/IVSPER_Analysis/Long_reads/07_Hybrid_assembly/"
    "Rotated_junction_mapping/Competitive_mapping"
)

VIRAL_CONTIG = "jcf7180000126179"

ONT_DEPTH_FILE = os.path.join(
    BASE,
    "ONT",
    "XwIV_primary.depth.tsv"
)

ILLUMINA_DEPTH_FILE = os.path.join(
    BASE,
    "Illumina",
    "XwIV_primary.depth.tsv"
)

ROTATION_INFO = os.path.join(
    BASE,
    "Reference",
    "XwIV_rotation_info.tsv"
)

OUTDIR = os.path.join(
    BASE,
    "Zero_depth_analysis"
)

os.makedirs(
    OUTDIR,
    exist_ok=True
)


# ============================================================
# OUTPUT FILES
# ============================================================

ONT_RUNS_FILE = os.path.join(
    OUTDIR,
    "ONT_zero_depth_runs.tsv"
)

ILLUMINA_RUNS_FILE = os.path.join(
    OUTDIR,
    "Illumina_zero_depth_runs.tsv"
)

JOINT_RUNS_FILE = os.path.join(
    OUTDIR,
    "ONT_Illumina_joint_zero_depth_runs.tsv"
)

SUMMARY_FILE = os.path.join(
    OUTDIR,
    "XwIV_zero_depth_summary.txt"
)

POSITION_TABLE = os.path.join(
    OUTDIR,
    "XwIV_base_by_base_depth_status.tsv"
)

PLOT_PREFIX = os.path.join(
    OUTDIR,
    "XwIV_zero_depth_genome_wide"
)

ZERO_ONLY_PREFIX = os.path.join(
    OUTDIR,
    "XwIV_zero_depth_tracks"
)


# ============================================================
# READ ROTATION INFORMATION
# ============================================================

rotation = {}

with open(ROTATION_INFO) as handle:

    header = handle.readline()

    for line in handle:

        fields = line.rstrip().split("\t")

        if len(fields) < 2:
            continue

        rotation[fields[0]] = fields[1]


GENOME_LENGTH = int(
    rotation["genome_length"]
)

JUNCTION = int(
    rotation["junction_after_position"]
)


print()
print("============================================================")
print("XwIV ZERO-DEPTH ANALYSIS")
print("============================================================")
print()

print(
    f"Viral contig: {VIRAL_CONTIG}"
)

print(
    f"Genome length: {GENOME_LENGTH:,} bp"
)

print(
    f"END -> START junction: "
    f"{JUNCTION:,}|{JUNCTION + 1:,}"
)


# ============================================================
# READ DEPTH FILE
#
# IMPORTANT:
# The depth files may contain positions from other scaffolds
# because the BAM header contains the complete hybrid assembly.
#
# Only lines corresponding to the XwIV contig are retained.
# ============================================================

def read_depth(filename):

    depth = np.zeros(
        GENOME_LENGTH,
        dtype=np.int64
    )

    seen = np.zeros(
        GENOME_LENGTH,
        dtype=bool
    )

    n_lines_total = 0
    n_lines_viral = 0

    with open(filename) as handle:

        for line in handle:

            if not line.strip():
                continue

            fields = line.rstrip().split("\t")

            if len(fields) < 3:
                continue

            n_lines_total += 1

            # ------------------------------------------------
            # Keep ONLY the XwIV contig
            # ------------------------------------------------

            if fields[0] != VIRAL_CONTIG:
                continue

            n_lines_viral += 1

            position = int(
                fields[1]
            )

            value = int(
                fields[2]
            )

            if (
                position < 1
                or
                position > GENOME_LENGTH
            ):

                raise RuntimeError(
                    f"Invalid XwIV position {position} "
                    f"in {filename}. Expected positions "
                    f"between 1 and {GENOME_LENGTH}."
                )

            if seen[position - 1]:

                raise RuntimeError(
                    f"Position {position} occurs more than once "
                    f"for {VIRAL_CONTIG} in {filename}."
                )

            depth[
                position - 1
            ] = value

            seen[
                position - 1
            ] = True


    print(
        f"  Total depth lines read: "
        f"{n_lines_total:,}"
    )

    print(
        f"  XwIV lines retained: "
        f"{n_lines_viral:,}"
    )


    # --------------------------------------------------------
    # All 105,371 XwIV positions should be represented
    # because depth was generated with -aa.
    # --------------------------------------------------------

    missing = int(
        np.sum(
            ~seen
        )
    )

    if missing != 0:

        missing_positions = (
            np.where(~seen)[0] + 1
        )

        preview = ", ".join(
            str(x)
            for x in missing_positions[:20]
        )

        raise RuntimeError(
            f"{filename}: {missing:,} XwIV positions "
            f"are absent from the depth file. "
            f"First missing positions: {preview}"
        )


    if n_lines_viral != GENOME_LENGTH:

        raise RuntimeError(
            f"{filename}: expected exactly "
            f"{GENOME_LENGTH:,} lines for {VIRAL_CONTIG}, "
            f"but found {n_lines_viral:,}."
        )


    return depth


# ============================================================
# LOAD DEPTH
# ============================================================

print()
print("Loading ONT depth...")

ont_depth = read_depth(
    ONT_DEPTH_FILE
)

print()
print("Loading Illumina depth...")

illumina_depth = read_depth(
    ILLUMINA_DEPTH_FILE
)


# ============================================================
# ZERO-DEPTH MASKS
# ============================================================

ont_zero = (
    ont_depth == 0
)

illumina_zero = (
    illumina_depth == 0
)

joint_zero = (
    ont_zero
    &
    illumina_zero
)


# ============================================================
# FIND CONSECUTIVE ZERO-DEPTH RUNS
# ============================================================

def find_runs(mask):

    """
    Return consecutive True runs as:
    start, end, length

    Coordinates are 1-based and inclusive.
    """

    runs = []

    in_run = False
    start = None

    for i, value in enumerate(
        mask,
        start=1
    ):

        if value and not in_run:

            start = i
            in_run = True

        elif not value and in_run:

            end = i - 1

            runs.append(
                (
                    start,
                    end,
                    end - start + 1
                )
            )

            in_run = False


    if in_run:

        end = len(mask)

        runs.append(
            (
                start,
                end,
                end - start + 1
            )
        )


    return runs


ont_runs = find_runs(
    ont_zero
)

illumina_runs = find_runs(
    illumina_zero
)

joint_runs = find_runs(
    joint_zero
)


# ============================================================
# DOES A ZERO-DEPTH RUN SPAN THE JUNCTION?
#
# The artificial junction lies BETWEEN positions:
#
#       JUNCTION | JUNCTION + 1
#
# A zero-depth run is therefore considered to span the
# junction only if BOTH flanking positions belong to the
# same continuous zero-depth run.
# ============================================================

def run_contains_junction(
    start,
    end
):

    return (
        start <= JUNCTION
        and
        end >= JUNCTION + 1
    )


# ============================================================
# CONVERT RUNS TO DATAFRAME
# ============================================================

def runs_to_dataframe(
    runs,
    dataset
):

    rows = []

    for i, (
        start,
        end,
        length
    ) in enumerate(
        runs,
        start=1
    ):

        rows.append(
            {
                "dataset": dataset,
                "run": i,
                "start": start,
                "end": end,
                "length_bp": length,
                "contains_END_START_junction":
                    run_contains_junction(
                        start,
                        end
                    )
            }
        )


    return pd.DataFrame(
        rows,
        columns=[
            "dataset",
            "run",
            "start",
            "end",
            "length_bp",
            "contains_END_START_junction"
        ]
    )


ont_df = runs_to_dataframe(
    ont_runs,
    "ONT"
)

illumina_df = runs_to_dataframe(
    illumina_runs,
    "Illumina"
)

joint_df = runs_to_dataframe(
    joint_runs,
    "ONT+Illumina"
)


# ============================================================
# RANK ZERO-DEPTH RUNS BY LENGTH
# ============================================================

def add_length_rank(df):

    df = df.copy()

    if len(df) == 0:

        df[
            "length_rank"
        ] = pd.Series(
            dtype="Int64"
        )

        return df


    df[
        "length_rank"
    ] = (
        df[
            "length_bp"
        ]
        .rank(
            method="min",
            ascending=False
        )
        .astype(int)
    )

    return df


ont_df = add_length_rank(
    ont_df
)

illumina_df = add_length_rank(
    illumina_df
)

joint_df = add_length_rank(
    joint_df
)


# ============================================================
# SAVE ZERO-DEPTH RUN TABLES
# ============================================================

ont_df.to_csv(
    ONT_RUNS_FILE,
    sep="\t",
    index=False
)

illumina_df.to_csv(
    ILLUMINA_RUNS_FILE,
    sep="\t",
    index=False
)

joint_df.to_csv(
    JOINT_RUNS_FILE,
    sep="\t",
    index=False
)


# ============================================================
# BASE-BY-BASE TABLE
# ============================================================

positions = np.arange(
    1,
    GENOME_LENGTH + 1
)


position_df = pd.DataFrame(
    {
        "contig": VIRAL_CONTIG,
        "position": positions,
        "ONT_depth": ont_depth,
        "Illumina_depth": illumina_depth,
        "ONT_zero": ont_zero,
        "Illumina_zero": illumina_zero,
        "joint_zero": joint_zero
    }
)


position_df.to_csv(
    POSITION_TABLE,
    sep="\t",
    index=False
)


# ============================================================
# SUMMARY FUNCTION
# ============================================================

def dataset_summary(
    name,
    depth,
    zero_mask,
    run_df
):

    zero_positions = int(
        np.sum(
            zero_mask
        )
    )

    zero_percent = (
        zero_positions
        /
        GENOME_LENGTH
        *
        100
    )

    number_runs = len(
        run_df
    )

    if number_runs > 0:

        longest = int(
            run_df[
                "length_bp"
            ].max()
        )

    else:

        longest = 0


    junction_rows = run_df[
        run_df[
            "contains_END_START_junction"
        ]
        ==
        True
    ]


    lines = []

    lines.append(
        name
    )

    lines.append(
        "-" * len(name)
    )

    lines.append(
        f"Zero-depth positions: "
        f"{zero_positions:,} / "
        f"{GENOME_LENGTH:,} "
        f"({zero_percent:.4f}%)"
    )

    lines.append(
        f"Number of zero-depth runs: "
        f"{number_runs:,}"
    )

    lines.append(
        f"Longest zero-depth run: "
        f"{longest:,} bp"
    )


    if depth is not None:

        lines.append(
            f"Depth at position {JUNCTION:,} "
            f"(immediately before junction): "
            f"{depth[JUNCTION - 1]}"
        )

        lines.append(
            f"Depth at position {JUNCTION + 1:,} "
            f"(immediately after junction): "
            f"{depth[JUNCTION]}"
        )


    if len(
        junction_rows
    ) == 0:

        lines.append(
            "No continuous zero-depth run spans "
            "both sides of the END -> START junction."
        )

    else:

        for _, row in junction_rows.iterrows():

            lines.append(
                "Zero-depth run spanning junction: "
                f"{int(row['start']):,}-"
                f"{int(row['end']):,}"
            )

            lines.append(
                "Junction zero-depth run length: "
                f"{int(row['length_bp']):,} bp"
            )

            lines.append(
                "Junction zero-depth run rank by length: "
                f"{int(row['length_rank'])} / "
                f"{number_runs}"
            )


    return lines


# ============================================================
# BUILD SUMMARY
# ============================================================

summary_lines = []

summary_lines.append(
    "XwIV ZERO-DEPTH ANALYSIS"
)

summary_lines.append(
    "=" * 60
)

summary_lines.append(
    ""
)

summary_lines.append(
    f"Viral contig: {VIRAL_CONTIG}"
)

summary_lines.append(
    f"Genome length: {GENOME_LENGTH:,} bp"
)

summary_lines.append(
    f"Artificial END -> START junction: "
    f"{JUNCTION:,}|{JUNCTION + 1:,}"
)

summary_lines.append(
    ""
)


summary_lines.extend(
    dataset_summary(
        "ONT",
        ont_depth,
        ont_zero,
        ont_df
    )
)

summary_lines.append(
    ""
)


summary_lines.extend(
    dataset_summary(
        "Illumina",
        illumina_depth,
        illumina_zero,
        illumina_df
    )
)

summary_lines.append(
    ""
)


summary_lines.extend(
    dataset_summary(
        "ONT + Illumina simultaneously",
        None,
        joint_zero,
        joint_df
    )
)


# ============================================================
# LIST ALL JOINT ZERO-DEPTH RUNS SORTED BY LENGTH
# ============================================================

summary_lines.append(
    ""
)

summary_lines.append(
    "JOINT ZERO-DEPTH RUNS SORTED BY LENGTH"
)

summary_lines.append(
    "-" * 60
)


if len(
    joint_df
) == 0:

    summary_lines.append(
        "No joint zero-depth runs detected."
    )

else:

    sorted_joint = joint_df.sort_values(
        [
            "length_bp",
            "start"
        ],
        ascending=[
            False,
            True
        ]
    )

    for _, row in sorted_joint.iterrows():

        marker = ""

        if row[
            "contains_END_START_junction"
        ]:

            marker = "  <-- END -> START junction"

        summary_lines.append(
            f"{int(row['start']):,}-"
            f"{int(row['end']):,}: "
            f"{int(row['length_bp']):,} bp"
            f"{marker}"
        )


# ============================================================
# WRITE SUMMARY
# ============================================================

with open(
    SUMMARY_FILE,
    "w"
) as out:

    out.write(
        "\n".join(
            summary_lines
        )
    )

    out.write(
        "\n"
    )


print()
print(
    "\n".join(
        summary_lines
    )
)


# ============================================================
# FIGURE 1
#
# Genome-wide depth + zero-depth tracks.
# ============================================================

fig, axes = plt.subplots(
    3,
    1,
    figsize=(
        16,
        8
    ),
    sharex=True,
    gridspec_kw={
        "height_ratios": [
            2,
            2,
            1.2
        ]
    }
)


# ------------------------------------------------------------
# ONT depth
# ------------------------------------------------------------

axes[0].plot(
    positions,
    ont_depth,
    linewidth=0.6
)

axes[0].axvline(
    JUNCTION + 0.5,
    linestyle="--",
    linewidth=1.5
)

axes[0].set_ylabel(
    "ONT depth"
)

axes[0].set_title(
    "ONT long reads"
)


# ------------------------------------------------------------
# Illumina depth
# ------------------------------------------------------------

axes[1].plot(
    positions,
    illumina_depth,
    linewidth=0.6
)

axes[1].axvline(
    JUNCTION + 0.5,
    linestyle="--",
    linewidth=1.5
)

axes[1].set_ylabel(
    "Illumina depth"
)

axes[1].set_title(
    "Illumina short reads"
)


# ------------------------------------------------------------
# Zero-depth tracks
# ------------------------------------------------------------

ax = axes[2]


def draw_zero_runs(
    axis,
    run_df,
    y
):

    # Genome baseline
    axis.hlines(
        y,
        1,
        GENOME_LENGTH,
        linewidth=3,
        alpha=0.20
    )

    # Zero-depth regions
    for _, row in run_df.iterrows():

        axis.hlines(
            y,
            int(
                row["start"]
            ),
            int(
                row["end"]
            ),
            linewidth=6
        )


draw_zero_runs(
    ax,
    ont_df,
    2
)

draw_zero_runs(
    ax,
    illumina_df,
    1
)

draw_zero_runs(
    ax,
    joint_df,
    0
)


ax.axvline(
    JUNCTION + 0.5,
    linestyle="--",
    linewidth=1.5
)


ax.set_yticks(
    [
        0,
        1,
        2
    ]
)

ax.set_yticklabels(
    [
        "Both",
        "Illumina",
        "ONT"
    ]
)

ax.set_ylim(
    -0.6,
    2.6
)

ax.set_ylabel(
    "Zero depth"
)

ax.set_xlabel(
    "Position on rotated XwIV genome (bp)"
)

ax.set_title(
    "Genome-wide zero-depth regions"
)


# ------------------------------------------------------------
# Common formatting
# ------------------------------------------------------------

for axis in axes:

    axis.set_xlim(
        1,
        GENOME_LENGTH
    )

    axis.grid(
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
    "Genome-wide assessment of zero read coverage in XwIV",
    fontsize=15
)


plt.tight_layout(
    rect=[
        0,
        0,
        1,
        0.96
    ]
)


for extension in [
    "pdf",
    "png",
    "svg"
]:

    kwargs = {
        "bbox_inches": "tight"
    }

    if extension == "png":

        kwargs[
            "dpi"
        ] = 300

    plt.savefig(
        f"{PLOT_PREFIX}.{extension}",
        **kwargs
    )


plt.close()


# ============================================================
# FIGURE 2
#
# Zero-depth tracks only.
# ============================================================

fig, ax = plt.subplots(
    figsize=(
        16,
        3.5
    )
)


draw_zero_runs(
    ax,
    ont_df,
    2
)

draw_zero_runs(
    ax,
    illumina_df,
    1
)

draw_zero_runs(
    ax,
    joint_df,
    0
)


ax.axvline(
    JUNCTION + 0.5,
    linestyle="--",
    linewidth=1.8
)


ax.text(
    JUNCTION + 0.5,
    1.02,
    "END → START",
    transform=ax.get_xaxis_transform(),
    ha="center",
    va="bottom",
    fontsize=9
)


ax.set_yticks(
    [
        0,
        1,
        2
    ]
)

ax.set_yticklabels(
    [
        "Both",
        "Illumina",
        "ONT"
    ]
)

ax.set_ylim(
    -0.6,
    2.6
)

ax.set_xlim(
    1,
    GENOME_LENGTH
)

ax.set_xlabel(
    "Position on rotated XwIV genome (bp)"
)

ax.set_title(
    "Zero-depth regions across the XwIV genome"
)

ax.grid(
    axis="x",
    alpha=0.15
)


plt.tight_layout()


for extension in [
    "pdf",
    "png",
    "svg"
]:

    kwargs = {
        "bbox_inches": "tight"
    }

    if extension == "png":

        kwargs[
            "dpi"
        ] = 300

    plt.savefig(
        f"{ZERO_ONLY_PREFIX}.{extension}",
        **kwargs
    )


plt.close()


# ============================================================
# FINISHED
# ============================================================

print()
print("============================================================")
print("OUTPUT")
print("============================================================")
print()

print(
    f"Summary:\n{SUMMARY_FILE}"
)

print()

print(
    f"ONT zero-depth runs:\n{ONT_RUNS_FILE}"
)

print()

print(
    f"Illumina zero-depth runs:\n{ILLUMINA_RUNS_FILE}"
)

print()

print(
    f"Joint zero-depth runs:\n{JOINT_RUNS_FILE}"
)

print()

print(
    f"Base-by-base table:\n{POSITION_TABLE}"
)

print()

print(
    f"Genome-wide figure:\n{PLOT_PREFIX}.pdf"
)

print()

print(
    f"Zero-depth track figure:\n{ZERO_ONLY_PREFIX}.pdf"
)

print()
print("Done.")
