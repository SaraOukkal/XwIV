#!/usr/bin/env python3

import pandas as pd
import os
import matplotlib.pyplot as plt
import seaborn as sns
from numpy import histogram
from dna_features_viewer import GraphicRecord, GraphicFeature
from matplotlib.patches import Patch

# === CONFIGURATION ===
input_dir = "/beegfs/data/soukkal/Thesis/Hymenoptera_Project/PDV/Results/Viral_genes/IVSPER_Analysis/xiphosomella_wirra_specimen2/ORF_prediction/"
tsv_all = os.path.join(input_dir, "scaffold3_ORFs_70AA.tsv")
tsv_filtered = os.path.join(input_dir, "scaffold3_ORFs_70AA_filtered.tsv")
annotation_path = os.path.join(input_dir, "ORFs_annotation.tsv")

# === COULEURS ===
type_colors = {
    "IVSPER": "#8A2BE2",         # violet foncé
    "IV_Circles": "#CBA8F5",     # violet clair
    "Eukaryote": "#FFD700",      # jaune
    "Virus": "#1E90FF",          # bleu
    "Prokaryote": "#DC143C"      # rouge
}
default_color = "#A9A9A9"   # gris

color_all = "#888888"
color_filtered = "#6495ED"
color_removed = "#FF6347"

# === LECTURE DES DONNÉES ===
orf_df = pd.read_csv(tsv_all, sep="\t")
filtered_df = pd.read_csv(tsv_filtered, sep="\t")
annotation_df = pd.read_csv(annotation_path, sep="\t")
annotation_df["ORFs"] = annotation_df["ORFs"].astype(str)

orf_df["strand"] = orf_df.apply(lambda row: "+" if row["end"] >= row["start"] else "-", axis=1)
orf_df["length"] = abs(orf_df["end"] - orf_df["start"]) + 1
filtered_df["length"] = abs(filtered_df["end"] - filtered_df["start"]) + 1
removed_df = orf_df[~orf_df["orf_number"].isin(filtered_df["orf_number"])].copy()

# === DICTIONNAIRE D'ANNOTATIONS ===
annot_dict = annotation_df.set_index("ORFs")[["Protein", "Type"]].to_dict("index")

# === PLOT 1 : Visualisation ORFs filtrés avec annotation ===
features = []
for _, row in filtered_df.iterrows():
    orf_id = str(row["orf_number"])
    strand = 1 if row["strand"] == "+" else -1
    entry = annot_dict.get(orf_id, None)

    if entry is not None:
        protein = entry.get("Protein", "")
        if isinstance(protein, float):
            protein = ""
        else:
            protein = str(protein).strip()
        label = f"{orf_id} {protein}" if protein else f"{orf_id}"
        orf_type = entry.get("Type", "")
    else:
        label = None
        orf_type = ""

    color = type_colors.get(orf_type, default_color)

    features.append(
        GraphicFeature(
            start=row["start"],
            end=row["end"],
            strand=strand,
            color=color,
            label=label
        )
    )

record = GraphicRecord(sequence_length=100523, features=features)
ax, _ = record.plot(figure_width=12)
xticks = list(range(0, 100001, 10000))
ax.set_xticks(xticks)
ax.set_xticklabels(xticks, rotation=45)
ax.set_title("Annotated ORFs on scaffold3")

legend_elements = [
    Patch(facecolor=color, edgecolor="black", label=typ)
    for typ, color in type_colors.items()
]
ax.legend(handles=legend_elements, loc="upper right", title="Type")

plt.tight_layout()
plt.savefig(os.path.join(input_dir, "scaffold3_ORFs_visualisation.svg"))
plt.close()

# === PLOT 2 : Distributions de longueur ===
fig, axes = plt.subplots(3, 1, figsize=(6, 10), sharex=True)

xmax = pd.concat([orf_df["length"], filtered_df["length"], removed_df["length"]]).max()
bin_edges = list(range(0, int(xmax + 200), 200))

all_counts, _ = histogram(orf_df["length"], bins=bin_edges)
filtered_counts, _ = histogram(filtered_df["length"], bins=bin_edges)
removed_counts, _ = histogram(removed_df["length"], bins=bin_edges)
ymax = max(all_counts.max(), filtered_counts.max(), removed_counts.max())

for ax, data, title, fill_color in zip(
    axes,
    [orf_df, filtered_df, removed_df],
    ["All ORFs", "Filtered ORFs", "Removed ORFs"],
    [color_all, color_filtered, color_removed]
):
    sns.histplot(
        data=data,
        x="length",
        bins=bin_edges,
        color=fill_color,
        edgecolor=fill_color,
        alpha=0.5,
        linewidth=0.7,
        ax=ax
    )
    ax.set_xlim(0, bin_edges[-1])
    ax.set_ylim(0, ymax + 1)
    ax.set_title(title)
    ax.set_ylabel("Count")
    ax.set_xlabel("ORF length (nt)")

plt.tight_layout()
plt.savefig(os.path.join(input_dir, "scaffold3_ORFs_length_distributions.svg"))
plt.close()

