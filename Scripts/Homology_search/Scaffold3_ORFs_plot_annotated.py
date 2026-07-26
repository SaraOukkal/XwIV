import pandas as pd
from dna_features_viewer import GraphicFeature, GraphicRecord
import matplotlib.pyplot as plt
import os

# === CONFIGURATION ===
ORF_TABLE = "/beegfs/data/soukkal/Thesis/Hymenoptera_Project/PDV/Results/Viral_genes/IVSPER_Analysis/xiphosomella_wirra_specimen2/ORF_prediction/scaffold3_ORFs_70AA.tsv"
OUTPUT_DIR = "/beegfs/data/soukkal/Thesis/Hymenoptera_Project/PDV/Results/Viral_genes/IVSPER_Analysis/xiphosomella_wirra_specimen2/ORF_prediction/"
SCAFFOLD_NAME = "scaffold3"
SCAFFOLD_LENGTH = 100523

# === LISTES DE CLASSIFICATION ===
ivspers = {f"ORF{x}" for x in [217, 249, 168, 48, 18, 259, 39, 202, 142, 54, 68, 70, 119, 222, 232, 271, 223, 248]}
ivspers_weak = {f"ORF{x}" for x in [183, 35, 235, 284, 258, 188]}
viraux = {f"ORF{x}" for x in [54, 190, 249, 77, 147, 100, 124, 120, 131]}
insectes = {f"ORF{x}" for x in [271, 205]}

# === COULEURS ===
color_map = {
    "ivspers": "#8A2BE2",        # violet
    "ivspers_weak": "#D8BFD8",   # violet clair
    "viraux": "#1E90FF",         # bleu
    "insectes": "#FFA500",       # orange
    "default": "#C0C0C0"         # gris
}

# === LECTURE DU TABLEAU ===
df = pd.read_csv(ORF_TABLE, sep="\t")

# === CONSTRUCTION DES FEATURES ===
features = []
for _, row in df.iterrows():
    label = f"ORF{row['orf_number']}"
    if label in ivspers:
        color = color_map["ivspers"]
        display_label = label
    elif label in ivspers_weak:
        color = color_map["ivspers_weak"]
        display_label = label
    elif label in viraux:
        color = color_map["viraux"]
        display_label = label
    elif label in insectes:
        color = color_map["insectes"]
        display_label = label
    else:
        color = color_map["default"]
        display_label = None  # Pas de trait vertical parasite

    features.append(GraphicFeature(
        start=min(row["start"], row["end"]),
        end=max(row["start"], row["end"]),
        strand=1 if row["start"] < row["end"] else -1,
        color=color,
        label=display_label
    ))

# === VISUALISATION ===
record = GraphicRecord(sequence_length=SCAFFOLD_LENGTH, features=features)
ax, _ = record.plot(figure_width=12)

# Graduation tous les 10kb
ax.set_title(SCAFFOLD_NAME)
ax.set_xticks(range(0, SCAFFOLD_LENGTH + 1, 10000))
ax.set_xticklabels([f"{x//1000} kb" for x in range(0, SCAFFOLD_LENGTH + 1, 10000)])
plt.tight_layout()

# === EXPORT DES FIGURES ===
out_base = os.path.join(OUTPUT_DIR, f"{SCAFFOLD_NAME}_ORFs_annotation_clean")
plt.savefig(f"{out_base}.png", dpi=300)
plt.savefig(f"{out_base}.svg")
plt.close()

