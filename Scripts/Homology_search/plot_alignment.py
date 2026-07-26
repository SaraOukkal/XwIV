from Bio import AlignIO
from dna_features_viewer import GraphicFeature, GraphicRecord
import matplotlib.pyplot as plt
import os

# === CONFIGURATION ===
ALIGNMENT_FILE = "/beegfs/data/soukkal/Thesis/Hymenoptera_Project/PDV/Results/Viral_genes/IVSPER_Analysis/xiphosomella_wirra_specimen2/clustalo_alignment/alignment_with_rc.aln"
OUTPUT_DIR = "/beegfs/data/soukkal/Thesis/Hymenoptera_Project/PDV/Results/Viral_genes/IVSPER_Analysis/xiphosomella_wirra_specimen2/clustalo_alignment/"
SCAFFOLD_REF = "scaffold3|size100523"
SCAFFOLD_LENGTH = 100523
COLORS = {
    "scaffold3|size100523": "#FFD700",    # jaune
    "scaffold187|size52917": "#FF0000",   # rouge
    "scaffold33401|size845_RC": "#00FF00",# vert
    "scaffold524|size39189": "#0000FF",   # bleu
}

# === CHARGEMENT DE L'ALIGNEMENT ===
alignment = AlignIO.read(ALIGNMENT_FILE, "fasta")

# Map nom -> séquence alignée
seqs = {record.id: str(record.seq) for record in alignment}
assert SCAFFOLD_REF in seqs, f"{SCAFFOLD_REF} non trouvé dans l'alignement."

# Position de référence (scaffold3)
ref_seq = seqs[SCAFFOLD_REF]
ref_pos = 0
ref_coords = []
for base in ref_seq:
    if base != "-":
        ref_coords.append(ref_pos)
        ref_pos += 1
    else:
        ref_coords.append(None)

# Extraction des blocs alignés
features = []

for name, aln_seq in seqs.items():
    if name == SCAFFOLD_REF:
        continue

    start = None
    for i, base in enumerate(aln_seq):
        if base != "-" and ref_coords[i] is not None:
            if start is None:
                start = ref_coords[i]
        else:
            if start is not None:
                end = ref_coords[i-1] + 1
                features.append(GraphicFeature(start=start, end=end, strand=0,
                                               color=COLORS.get(name, "#999999"), label=name))
                start = None
    if start is not None:
        end = ref_coords[-1] + 1
        features.append(GraphicFeature(start=start, end=end, strand=0,
                                       color=COLORS.get(name, "#999999"), label=name))

# === Ajout de la bande de référence ===
features.insert(0, GraphicFeature(start=0, end=SCAFFOLD_LENGTH, strand=0,
                                  color=COLORS[SCAFFOLD_REF], label="scaffold3 (référence)"))

# === PLOT ===
record = GraphicRecord(sequence_length=SCAFFOLD_LENGTH, features=features)
ax, _ = record.plot(figure_width=12)
ax.set_title("Alignement des scaffolds sur scaffold3")
plt.tight_layout()

# === EXPORT ===
plt.savefig(os.path.join(OUTPUT_DIR, "scaffold3_alignment_blocks.png"))
plt.savefig(os.path.join(OUTPUT_DIR, "scaffold3_alignment_blocks.svg"))
plt.close()

