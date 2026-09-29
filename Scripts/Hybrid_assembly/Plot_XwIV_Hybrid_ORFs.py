#!/usr/bin/env python3

import os
import re
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.patches import Patch
from dna_features_viewer import GraphicFeature, GraphicRecord


# ============================================================
# CONFIGURATION
# ============================================================

INPUT_DIR = (
    "/beegfs/data/soukkal/Thesis/Hymenoptera_Project/PDV/Results/"
    "Viral_genes/IVSPER_Analysis/Long_reads/07_Hybrid_assembly/"
    "ORF_prediction"
)

ORF_TABLE = os.path.join(
    INPUT_DIR,
    "XwIV_hybrid_final_ORFs.tsv"
)

ANNOTATION_TABLE = os.path.join(
    INPUT_DIR,
    "ORF_Annotation.tsv"
)

OUTPUT_PREFIX = os.path.join(
    INPUT_DIR,
    "XwIV_hybrid_ORF_annotation"
)

OUTPUT_SUMMARY = os.path.join(
    INPUT_DIR,
    "XwIV_hybrid_ORF_annotation_plot_table.tsv"
)

GENOME_NAME = "XwIV"
GENOME_LENGTH = 105371


# ============================================================
# COLORS
# ============================================================

COLOR_IVSPER = "#8E44AD"
COLOR_VIRUS = "#3498DB"
COLOR_EUKARYOTE = "#E67E22"
COLOR_PROKARYOTE = "#F1C40F"
COLOR_OTHER = "#BDBDBD"


# ============================================================
# FUNCTIONS
# ============================================================

def extract_orf_number(value):
    """
    Extract numeric ORF identifier from values such as:
    XwIV_ORF_12
    ORF12
    ORF_12
    12
    """

    if pd.isna(value):
        return None

    value = str(value).strip()

    match = re.search(r"(\d+)$", value)

    if match:
        return int(match.group(1))

    return None


def clean_text(value):
    """
    Normalize text and remove whitespace / invisible characters.
    """

    if pd.isna(value):
        return ""

    value = str(value)

    # Remove common invisible Unicode characters
    value = value.replace("\u200b", "")
    value = value.replace("\u200c", "")
    value = value.replace("\u200d", "")
    value = value.replace("\ufeff", "")
    value = value.replace("\xa0", " ")

    value = value.strip()

    if value.lower() in {
        "",
        "na",
        "nan",
        "none",
        "-"
    }:
        return ""

    return value


def normalize_strand(value):
    """
    Convert different strand representations to +1 or -1.
    """

    if pd.isna(value):
        raise ValueError("Missing strand value")

    value = str(value)

    # Remove whitespace and invisible characters
    value = value.replace("\u200b", "")
    value = value.replace("\u200c", "")
    value = value.replace("\u200d", "")
    value = value.replace("\ufeff", "")
    value = value.replace("\xa0", "")
    value = value.strip()

    # Normalize Unicode minus/dash characters
    value = value.replace("−", "-")
    value = value.replace("–", "-")
    value = value.replace("—", "-")

    if value in {
        "+",
        "1",
        "+1"
    }:
        return 1

    if value in {
        "-",
        "-1"
    }:
        return -1

    raise ValueError(
        f"Unknown strand value: {repr(value)}"
    )


def classify_annotation(type_value):
    """
    Assign plotting category according to annotation priority.

    Priority:
    1. IVSPER / IV segment gene Diadegma
    2. Virus
    3. Eukaryote
    4. Bacteria / Archaea
    5. Other / unannotated
    """

    if pd.isna(type_value):
        return "Unannotated"

    type_value = str(type_value).strip()

    if not type_value:
        return "Unannotated"

    annotations = [
        x.strip().lower()
        for x in type_value.split("/")
        if x.strip()
    ]

    # Highest priority
    if any(
        ("ivsper" in x)
        or ("iv segment gene diadegma" in x)
        for x in annotations
    ):
        return "IVSPER / IV segment gene"

    # Virus
    if any(
        "virus" in x
        for x in annotations
    ):
        return "Virus"

    # Eukaryote
    if any(
        "eukaryote" in x
        for x in annotations
    ):
        return "Eukaryote"

    # Bacteria / Archaea
    if any(
        ("bacteria" in x)
        or ("archea" in x)
        or ("archaea" in x)
        for x in annotations
    ):
        return "Bacteria / Archaea"

    return "Other / unannotated"


def category_to_color(category):

    color_dict = {
        "IVSPER / IV segment gene": COLOR_IVSPER,
        "Virus": COLOR_VIRUS,
        "Eukaryote": COLOR_EUKARYOTE,
        "Bacteria / Archaea": COLOR_PROKARYOTE,
        "Other / unannotated": COLOR_OTHER,
        "Unannotated": COLOR_OTHER
    }

    return color_dict.get(
        category,
        COLOR_OTHER
    )


# ============================================================
# CHECK INPUT FILES
# ============================================================

for input_file in [
    ORF_TABLE,
    ANNOTATION_TABLE
]:

    if not os.path.isfile(input_file):
        raise FileNotFoundError(
            f"Input file not found: {input_file}"
        )


# ============================================================
# READ ORF TABLE
# ============================================================

orfs = pd.read_csv(
    ORF_TABLE,
    sep="\t",
    dtype=str
)


required_orf_columns = {
    "New_ORF",
    "start",
    "end",
    "strand"
}


missing = required_orf_columns - set(orfs.columns)

if missing:
    raise ValueError(
        "Missing columns in XwIV_hybrid_final_ORFs.tsv: "
        + ", ".join(sorted(missing))
    )


orfs["start"] = pd.to_numeric(
    orfs["start"],
    errors="raise"
)

orfs["end"] = pd.to_numeric(
    orfs["end"],
    errors="raise"
)


# ============================================================
# EXTRACT ORF NUMBERS
# ============================================================

orfs["ORF_number"] = orfs["New_ORF"].apply(
    extract_orf_number
)


if orfs["ORF_number"].isna().any():

    bad = orfs.loc[
        orfs["ORF_number"].isna(),
        "New_ORF"
    ].tolist()

    raise ValueError(
        "Could not extract ORF numbers from: "
        + ", ".join(map(str, bad))
    )


orfs["ORF_number"] = orfs["ORF_number"].astype(int)


# ============================================================
# NORMALIZE STRAND
# ============================================================

strand_values = []

for _, row in orfs.iterrows():

    try:

        strand_values.append(
            normalize_strand(
                row["strand"]
            )
        )

    except ValueError:

        raise ValueError(
            f"Unknown strand value for "
            f"{row['New_ORF']}: "
            f"{repr(row['strand'])}"
        )


orfs["strand_numeric"] = strand_values


# ============================================================
# READ ANNOTATION TABLE
# ============================================================

annotations = pd.read_csv(
    ANNOTATION_TABLE,
    sep="\t",
    dtype=str
)


required_annotation_columns = {
    "ORF",
    "Gene name",
    "Type"
}


missing = required_annotation_columns - set(
    annotations.columns
)

if missing:
    raise ValueError(
        "Missing columns in ORF_Annotation.tsv: "
        + ", ".join(sorted(missing))
    )


annotations["ORF_number"] = annotations["ORF"].apply(
    extract_orf_number
)


if annotations["ORF_number"].isna().any():

    bad = annotations.loc[
        annotations["ORF_number"].isna(),
        "ORF"
    ].tolist()

    raise ValueError(
        "Could not extract ORF numbers from annotation table: "
        + ", ".join(map(str, bad))
    )


annotations["ORF_number"] = (
    annotations["ORF_number"]
    .astype(int)
)


# ============================================================
# CHECK DUPLICATE ANNOTATIONS
# ============================================================

duplicates = annotations[
    annotations["ORF_number"].duplicated(
        keep=False
    )
]


if not duplicates.empty:

    duplicate_numbers = sorted(
        duplicates["ORF_number"]
        .unique()
        .tolist()
    )

    raise ValueError(
        "Multiple annotation rows were found for ORF(s): "
        + ", ".join(map(str, duplicate_numbers))
    )


# ============================================================
# MERGE ORFs WITH ANNOTATIONS
# ============================================================

df = orfs.merge(
    annotations,
    on="ORF_number",
    how="left",
    suffixes=(
        "",
        "_annotation"
    )
)


# ============================================================
# ASSIGN ANNOTATION CATEGORY
# ============================================================

df["Annotation_category"] = df["Type"].apply(
    classify_annotation
)


df["Color"] = df[
    "Annotation_category"
].apply(
    category_to_color
)


# ============================================================
# ASSIGN LABELS
# ============================================================

def make_label(row):

    gene_name = clean_text(
        row.get(
            "Gene name",
            ""
        )
    )

    if gene_name:
        return gene_name

    return str(
        int(row["ORF_number"])
    )


df["Plot_label"] = df.apply(
    make_label,
    axis=1
)


# ============================================================
# SORT BY GENOMIC POSITION
# ============================================================

df["plot_start"] = df[
    ["start", "end"]
].min(axis=1)

df["plot_end"] = df[
    ["start", "end"]
].max(axis=1)


df = df.sort_values(
    by=[
        "plot_start",
        "plot_end"
    ]
).reset_index(drop=True)


# ============================================================
# CHECK COORDINATES
# ============================================================

if (
    (df["plot_start"] < 0).any()
    or
    (df["plot_end"] > GENOME_LENGTH).any()
):

    raise ValueError(
        "At least one ORF falls outside the expected "
        f"XwIV genome length ({GENOME_LENGTH} bp)."
    )


# ============================================================
# BUILD GRAPHIC FEATURES
# ============================================================

features = []


for _, row in df.iterrows():

    feature = GraphicFeature(
        start=int(row["plot_start"]),
        end=int(row["plot_end"]),
        strand=int(row["strand_numeric"]),
        color=row["Color"],
        label=row["Plot_label"]
    )

    features.append(feature)


# ============================================================
# CREATE GRAPHIC RECORD
# ============================================================

record = GraphicRecord(
    sequence_length=GENOME_LENGTH,
    features=features
)


# ============================================================
# PLOT
# ============================================================

ax, _ = record.plot(
    figure_width=18
)


ax.set_title(
    "XwIV genome – predicted ORFs",
    fontsize=14,
    fontweight="bold",
    pad=20
)


ax.set_xlabel(
    "Genomic position (bp)"
)


# ============================================================
# LEGEND
# ============================================================

legend_elements = [

    Patch(
        facecolor=COLOR_IVSPER,
        edgecolor="black",
        label="IVSPER / IV segment gene"
    ),

    Patch(
        facecolor=COLOR_VIRUS,
        edgecolor="black",
        label="Virus"
    ),

    Patch(
        facecolor=COLOR_EUKARYOTE,
        edgecolor="black",
        label="Eukaryote"
    ),

    Patch(
        facecolor=COLOR_PROKARYOTE,
        edgecolor="black",
        label="Bacteria / Archaea"
    ),

    Patch(
        facecolor=COLOR_OTHER,
        edgecolor="black",
        label="Other / unannotated"
    )
]


ax.legend(
    handles=legend_elements,
    loc="upper center",
    bbox_to_anchor=(
        0.5,
        -0.15
    ),
    ncol=5,
    frameon=False,
    fontsize=9
)


plt.tight_layout()


# ============================================================
# EXPORT FIGURES
# ============================================================

plt.savefig(
    OUTPUT_PREFIX + ".pdf",
    bbox_inches="tight"
)

plt.savefig(
    OUTPUT_PREFIX + ".svg",
    bbox_inches="tight"
)

plt.savefig(
    OUTPUT_PREFIX + ".png",
    dpi=300,
    bbox_inches="tight"
)

plt.close()


# ============================================================
# EXPORT PLOTTING TABLE
# ============================================================

output_columns = [
    "ORF_number",
    "New_ORF",
    "start",
    "end",
    "strand",
    "strand_numeric",
    "Plot_label",
    "Gene name",
    "Type",
    "Annotation_category",
    "Color"
]


output_columns = [
    x
    for x in output_columns
    if x in df.columns
]


df[
    output_columns
].to_csv(
    OUTPUT_SUMMARY,
    sep="\t",
    index=False
)


# ============================================================
# SUMMARY
# ============================================================

print()
print("============================================")
print("XwIV ORF ANNOTATION PLOT")
print("============================================")

print(
    f"Total ORFs plotted: {len(df)}"
)

print(
    f"ORFs present in annotation table: "
    f"{df['ORF'].notna().sum()}"
)

print()
print("Categories:")

category_counts = (
    df["Annotation_category"]
    .value_counts()
)

for category, count in category_counts.items():

    print(
        f"  {category}: {count}"
    )


print()
print("Figures:")
print(
    OUTPUT_PREFIX + ".pdf"
)
print(
    OUTPUT_PREFIX + ".svg"
)
print(
    OUTPUT_PREFIX + ".png"
)

print()
print("Plot annotation table:")
print(
    OUTPUT_SUMMARY
)

print()
print("Done.")
