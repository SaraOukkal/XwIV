#!/usr/bin/env python3

import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

from scipy.stats import mannwhitneyu


# ============================================================
# CONFIGURATION
# ============================================================

BASE = (
    "/beegfs/data/soukkal/Thesis/Hymenoptera_Project/PDV/Results/"
    "Viral_genes/IVSPER_Analysis/Long_reads/07_Hybrid_assembly"
)

ORF_FILE = os.path.join(
    BASE,
    "ORF_prediction",
    "XwIV_hybrid_final_ORFs.tsv"
)

COMPLEXITY_DIR = os.path.join(
    BASE,
    "Terminal_complexity"
)

DUST_FILE = os.path.join(
    COMPLEXITY_DIR,
    "XwIV_DUST_regions.tsv"
)

TRF_FILE = os.path.join(
    COMPLEXITY_DIR,
    "XwIV_TRF_regions.tsv"
)

OUTDIR = os.path.join(
    COMPLEXITY_DIR,
    "Terminal_region_statistics"
)

os.makedirs(
    OUTDIR,
    exist_ok=True
)


# ============================================================
# GENOME DEFINITION
# ============================================================

GENOME_LENGTH = 105371

WINDOW_SIZE = 1000

# Terminal regions defined from structural analyses
LEFT_TERMINAL_START = 1
LEFT_TERMINAL_END = 9691

RIGHT_TERMINAL_START = 102655
RIGHT_TERMINAL_END = 105371


# ============================================================
# OUTPUT FILES
# ============================================================

WINDOW_TABLE = os.path.join(
    OUTDIR,
    "XwIV_1kb_window_statistics.tsv"
)

TEST_TABLE = os.path.join(
    OUTDIR,
    "XwIV_terminal_vs_core_statistical_tests.tsv"
)

SUMMARY_TABLE = os.path.join(
    OUTDIR,
    "XwIV_terminal_vs_core_summary.tsv"
)

PLOT_PREFIX = os.path.join(
    OUTDIR,
    "XwIV_terminal_vs_core_comparison"
)


# ============================================================
# HELPER FUNCTIONS
# ============================================================

def merge_intervals(intervals):
    """
    Merge overlapping genomic intervals.

    Coordinates are treated as 1-based and inclusive.
    """

    if not intervals:
        return []

    intervals = sorted(
        intervals,
        key=lambda x: (x[0], x[1])
    )

    merged = [
        list(intervals[0])
    ]

    for start, end in intervals[1:]:

        previous = merged[-1]

        if start <= previous[1] + 1:

            previous[1] = max(
                previous[1],
                end
            )

        else:

            merged.append(
                [start, end]
            )

    return [
        tuple(x)
        for x in merged
    ]


def overlap_length(
    start1,
    end1,
    start2,
    end2
):
    """
    Inclusive overlap length between two intervals.
    """

    start = max(
        start1,
        start2
    )

    end = min(
        end1,
        end2
    )

    if end < start:
        return 0

    return end - start + 1


def covered_bases(
    window_start,
    window_end,
    intervals
):
    """
    Number of unique bases in a window covered by intervals.

    Intervals should already be merged.
    """

    total = 0

    for start, end in intervals:

        if end < window_start:
            continue

        if start > window_end:
            break

        total += overlap_length(
            window_start,
            window_end,
            start,
            end
        )

    return total


# ============================================================
# LOAD ORFs
# ============================================================

print()
print("============================================================")
print("LOADING ORFs")
print("============================================================")
print()

orfs = pd.read_csv(
    ORF_FILE,
    sep="\t"
)

required_orf_columns = {
    "start",
    "end"
}

if not required_orf_columns.issubset(
    orfs.columns
):

    raise RuntimeError(
        "ORF file must contain columns 'start' and 'end'. "
        f"Columns found: {list(orfs.columns)}"
    )


orf_intervals = []

for _, row in orfs.iterrows():

    start = int(
        row["start"]
    )

    end = int(
        row["end"]
    )

    # Coordinates are genomic coordinates:
    # start < end even for minus-strand ORFs.
    start, end = sorted(
        [start, end]
    )

    orf_intervals.append(
        (
            start,
            end
        )
    )


orf_intervals = merge_intervals(
    orf_intervals
)


print(
    f"ORFs loaded: {len(orfs)}"
)

print(
    f"Merged coding intervals: {len(orf_intervals)}"
)


# ============================================================
# LOAD DUST REGIONS
# ============================================================

print()
print("============================================================")
print("LOADING DUST REGIONS")
print("============================================================")
print()


dust = pd.read_csv(
    DUST_FILE,
    sep="\t"
)

print(
    f"DUST columns: {list(dust.columns)}"
)


# Identify start/end columns robustly

dust_columns_lower = {
    str(c).lower(): c
    for c in dust.columns
}

if (
    "start" not in dust_columns_lower
    or
    "end" not in dust_columns_lower
):

    raise RuntimeError(
        "Could not identify start/end columns in DUST file."
    )


dust_start_col = dust_columns_lower[
    "start"
]

dust_end_col = dust_columns_lower[
    "end"
]


dust_intervals = []

for _, row in dust.iterrows():

    start = int(
        row[dust_start_col]
    )

    end = int(
        row[dust_end_col]
    )

    start, end = sorted(
        [start, end]
    )

    dust_intervals.append(
        (
            start,
            end
        )
    )


dust_intervals = merge_intervals(
    dust_intervals
)


print(
    f"Merged DUST intervals: {len(dust_intervals)}"
)


# ============================================================
# LOAD TRF REGIONS
# ============================================================

print()
print("============================================================")
print("LOADING TRF REGIONS")
print("============================================================")
print()


trf = pd.read_csv(
    TRF_FILE,
    sep="\t"
)

print(
    f"TRF columns: {list(trf.columns)}"
)


trf_columns_lower = {
    str(c).lower(): c
    for c in trf.columns
}

if (
    "start" not in trf_columns_lower
    or
    "end" not in trf_columns_lower
):

    raise RuntimeError(
        "Could not identify start/end columns in TRF file."
    )


trf_start_col = trf_columns_lower[
    "start"
]

trf_end_col = trf_columns_lower[
    "end"
]


trf_intervals = []

for _, row in trf.iterrows():

    start = int(
        row[trf_start_col]
    )

    end = int(
        row[trf_end_col]
    )

    start, end = sorted(
        [start, end]
    )

    trf_intervals.append(
        (
            start,
            end
        )
    )


trf_intervals = merge_intervals(
    trf_intervals
)


print(
    f"Merged TRF intervals: {len(trf_intervals)}"
)


# ============================================================
# GENERATE COMPLETE 1-kb WINDOWS
# ============================================================

print()
print("============================================================")
print("GENERATING 1-kb WINDOWS")
print("============================================================")
print()


windows = []

window_number = 0

for start in range(
    1,
    GENOME_LENGTH + 1,
    WINDOW_SIZE
):

    end = min(
        start + WINDOW_SIZE - 1,
        GENOME_LENGTH
    )

    length = end - start + 1

    window_number += 1


    # --------------------------------------------------------
    # Classify window
    # --------------------------------------------------------

    if length < WINDOW_SIZE:

        region = "excluded_incomplete"

    elif (
        start >= LEFT_TERMINAL_START
        and
        end <= LEFT_TERMINAL_END
    ):

        region = "terminal"

    elif (
        start >= RIGHT_TERMINAL_START
        and
        end <= RIGHT_TERMINAL_END
    ):

        region = "terminal"

    elif (
        start > LEFT_TERMINAL_END
        and
        end < RIGHT_TERMINAL_START
    ):

        region = "core"

    else:

        # Window overlaps a terminal/core boundary.
        region = "excluded_boundary"


    # --------------------------------------------------------
    # Coverage
    # --------------------------------------------------------

    coding_bp = covered_bases(
        start,
        end,
        orf_intervals
    )

    dust_bp = covered_bases(
        start,
        end,
        dust_intervals
    )

    trf_bp = covered_bases(
        start,
        end,
        trf_intervals
    )


    windows.append(
        {
            "window": window_number,
            "start": start,
            "end": end,
            "length": length,
            "region": region,

            "coding_bp": coding_bp,
            "coding_percent": (
                coding_bp
                /
                length
                *
                100
            ),

            "low_complexity_bp": dust_bp,
            "low_complexity_percent": (
                dust_bp
                /
                length
                *
                100
            ),

            "tandem_repeat_bp": trf_bp,
            "tandem_repeat_percent": (
                trf_bp
                /
                length
                *
                100
            )
        }
    )


windows = pd.DataFrame(
    windows
)


windows.to_csv(
    WINDOW_TABLE,
    sep="\t",
    index=False
)


# ============================================================
# WINDOWS USED FOR TESTS
# ============================================================

test_windows = windows[
    windows["region"].isin(
        [
            "terminal",
            "core"
        ]
    )
].copy()


terminal = test_windows[
    test_windows["region"]
    ==
    "terminal"
].copy()


core = test_windows[
    test_windows["region"]
    ==
    "core"
].copy()


print(
    f"Terminal 1-kb windows: {len(terminal)}"
)

print(
    f"Core 1-kb windows: {len(core)}"
)

print()

print("Excluded windows:")

print(
    windows["region"].value_counts()
)


# ============================================================
# STATISTICAL TESTS
# ============================================================

print()
print("============================================================")
print("STATISTICAL TESTS")
print("============================================================")
print()


variables = [
    (
        "ORF density",
        "coding_percent"
    ),
    (
        "Low complexity",
        "low_complexity_percent"
    ),
    (
        "Tandem repeats",
        "tandem_repeat_percent"
    )
]


results = []


for name, column in variables:

    terminal_values = terminal[
        column
    ].to_numpy()

    core_values = core[
        column
    ].to_numpy()


    statistic, pvalue = mannwhitneyu(
        terminal_values,
        core_values,
        alternative="two-sided"
    )


    result = {
        "variable": name,

        "terminal_n": len(
            terminal_values
        ),

        "core_n": len(
            core_values
        ),

        "terminal_mean_percent": np.mean(
            terminal_values
        ),

        "terminal_median_percent": np.median(
            terminal_values
        ),

        "terminal_Q1_percent": np.percentile(
            terminal_values,
            25
        ),

        "terminal_Q3_percent": np.percentile(
            terminal_values,
            75
        ),

        "core_mean_percent": np.mean(
            core_values
        ),

        "core_median_percent": np.median(
            core_values
        ),

        "core_Q1_percent": np.percentile(
            core_values,
            25
        ),

        "core_Q3_percent": np.percentile(
            core_values,
            75
        ),

        "mann_whitney_U": statistic,

        "p_value": pvalue
    }

    results.append(
        result
    )


results = pd.DataFrame(
    results
)


# ============================================================
# BENJAMINI-HOCHBERG FDR
# ============================================================

pvalues = results[
    "p_value"
].to_numpy()


order = np.argsort(
    pvalues
)

ranked = pvalues[
    order
]


m = len(
    ranked
)


adjusted = np.empty(
    m,
    dtype=float
)


previous = 1.0


for i in range(
    m - 1,
    -1,
    -1
):

    rank = i + 1

    value = (
        ranked[i]
        *
        m
        /
        rank
    )

    value = min(
        value,
        previous,
        1.0
    )

    adjusted[i] = value

    previous = value


fdr = np.empty(
    m,
    dtype=float
)

fdr[
    order
] = adjusted


results[
    "FDR_BH"
] = fdr


results.to_csv(
    TEST_TABLE,
    sep="\t",
    index=False
)


# ============================================================
# SUMMARY TABLE
# ============================================================

summary_rows = []


for name, column in variables:

    for region_name, dataframe in [
        (
            "terminal",
            terminal
        ),
        (
            "core",
            core
        )
    ]:

        values = dataframe[
            column
        ].to_numpy()

        summary_rows.append(
            {
                "variable": name,
                "region": region_name,
                "n_windows": len(
                    values
                ),
                "mean_percent": np.mean(
                    values
                ),
                "median_percent": np.median(
                    values
                ),
                "Q1_percent": np.percentile(
                    values,
                    25
                ),
                "Q3_percent": np.percentile(
                    values,
                    75
                ),
                "min_percent": np.min(
                    values
                ),
                "max_percent": np.max(
                    values
                )
            }
        )


summary = pd.DataFrame(
    summary_rows
)


summary.to_csv(
    SUMMARY_TABLE,
    sep="\t",
    index=False
)


# ============================================================
# PRINT RESULTS
# ============================================================

print()

for _, row in results.iterrows():

    print(
        row["variable"]
    )

    print(
        "  Terminal median: "
        f"{row['terminal_median_percent']:.2f}%"
    )

    print(
        "  Core median:     "
        f"{row['core_median_percent']:.2f}%"
    )

    print(
        "  Mann-Whitney U:  "
        f"{row['mann_whitney_U']:.2f}"
    )

    print(
        "  P-value:         "
        f"{row['p_value']:.6g}"
    )

    print(
        "  FDR:             "
        f"{row['FDR_BH']:.6g}"
    )

    print()


# ============================================================
# PLOT
# ============================================================

fig, axes = plt.subplots(
    1,
    3,
    figsize=(
        13,
        5
    )
)


for ax, (
    name,
    column
) in zip(
    axes,
    variables
):

    terminal_values = terminal[
        column
    ].to_numpy()

    core_values = core[
        column
    ].to_numpy()


    # Boxplots
    ax.boxplot(
        [
            core_values,
            terminal_values
        ],
        labels=[
            "Core",
            "Terminal"
        ],
        showfliers=False
    )


    # Individual windows
    rng = np.random.default_rng(
        42
    )

    for x, values in [
        (
            1,
            core_values
        ),
        (
            2,
            terminal_values
        )
    ]:

        jitter = rng.normal(
            0,
            0.045,
            len(values)
        )

        ax.scatter(
            np.full(
                len(values),
                x
            )
            +
            jitter,
            values,
            s=18,
            alpha=0.65,
            zorder=3
        )


    result = results[
        results["variable"]
        ==
        name
    ].iloc[0]


    ax.set_title(
        name
    )

    ax.set_ylabel(
        "% of 1-kb window"
    )

    ax.text(
        0.5,
        1.02,
        (
            f"P = {result['p_value']:.3g}; "
            f"FDR = {result['FDR_BH']:.3g}"
        ),
        transform=ax.transAxes,
        ha="center",
        va="bottom",
        fontsize=9
    )


fig.suptitle(
    "Comparison of XwIV terminal regions with the genome core",
    fontsize=14
)


plt.tight_layout(
    rect=[
        0,
        0,
        1,
        0.94
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
# FINISHED
# ============================================================

print()
print("============================================================")
print("OUTPUT")
print("============================================================")
print()

print(
    f"Window data:\n{WINDOW_TABLE}"
)

print()

print(
    f"Statistical tests:\n{TEST_TABLE}"
)

print()

print(
    f"Summary:\n{SUMMARY_TABLE}"
)

print()

print(
    f"Figure:\n{PLOT_PREFIX}.pdf"
)

print()
print("Done.")
