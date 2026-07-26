#!/usr/bin/env python3

# Stacked alignment blocks on a fixed viral-genome reference, with per-contig colors and lane stacking.
# Groups (top to bottom): sp1 scaffolds, sp1 contigs, sp2 scaffolds, sp2 contigs.
# Reference drawn in grey on top. Labels are moved to a legend to avoid clutter.

from Bio import AlignIO
from dna_features_viewer import GraphicFeature, GraphicRecord
import matplotlib.pyplot as plt
import matplotlib as mpl
import os
from collections import defaultdict, OrderedDict

# =========================
# CONFIG
# =========================
ALIGNMENT_FILE = "/beegfs/data/soukkal/Thesis/Hymenoptera_Project/PDV/Results/Viral_genes/IVSPER_Analysis/Viral_Genome/Alignements/alignment_on_viral_genome.aln"
OUTPUT_DIR = "/beegfs/data/soukkal/Thesis/Hymenoptera_Project/PDV/Results/Viral_genes/IVSPER_Analysis/Viral_Genome/Alignements/"
os.makedirs(OUTPUT_DIR, exist_ok=True)

# Force reference header
REF_NAME = "scaffold3_circularized_manual"

# Known scaffold IDs for grouping and coloring (edit as needed)
SP1_SCAFFOLDS_SET = {
    "scaffold187|size52917",
    "scaffold33401|size845",
    "scaffold33401|size845_RC",
    "scaffold524|size39189",
}
SP2_SCAFFOLDS_SET = {
    "scaffold3|size100523"  # if present in the MSA apart from reference
}

# Colors
COLOR_REF = "#B0B0B0"  # grey for reference
SCAFFOLD_COLORS = {
    "scaffold187|size52917": "#FF0000",
    "scaffold33401|size845": "#00AA00",
    "scaffold33401|size845_RC": "#00AA00",
    "scaffold524|size39189": "#0000FF",
    "scaffold3|size100523": "#FFD700",
}
# Contig color palettes (cycled)
CONTIG_CMAP = mpl.cm.get_cmap("tab20")  # distinct qualitative palette
CONTIG_MAX_COLORS = 20  # tab20 has 20 distinct colors

# Visual parameters
BASE_FIGURE_WIDTH = 16
PER_ROW_HEIGHT = 1.4
LANE_PADDING = 10  # bp between features in the same lane
FEATURE_EDGE_COLOR = None  # set to "black" if you want borders around blocks

TITLE = "Alignment to viral genome (stacked lanes; unique colors per contig)"

# =========================
# LOAD ALIGNMENT
# =========================
alignment = AlignIO.read(ALIGNMENT_FILE, "fasta")
seqs = {rec.id: str(rec.seq) for rec in alignment}

if REF_NAME not in seqs:
    raise ValueError(f"Reference '{REF_NAME}' not found in alignment: {ALIGNMENT_FILE}")

ref_seq = seqs[REF_NAME]

# Map alignment columns to reference ungapped coordinates
ref_coords = []
ref_pos = 0
for b in ref_seq:
    if b != "-":
        ref_coords.append(ref_pos)
        ref_pos += 1
    else:
        ref_coords.append(None)

REF_LEN = ref_pos
if REF_LEN == 0:
    raise ValueError("Reference ungapped length is zero. Check the alignment.")

# =========================
# GROUPING AND COLORING
# =========================
def base_name(name: str) -> str:
    return name[:-3] if name.endswith("_RC") else name

def is_sp1_contig(name: str) -> bool:
    return name.startswith("sp1_contig|")

def is_sp2_contig(name: str) -> bool:
    return name.startswith("sp2_contig|")

def is_sp1_scaffold(name: str) -> bool:
    return base_name(name) in SP1_SCAFFOLDS_SET

def is_sp2_scaffold(name: str) -> bool:
    return base_name(name) in SP2_SCAFFOLDS_SET

def group_of(name: str) -> str:
    if is_sp1_scaffold(name):
        return "sp1_scaffolds"
    if is_sp1_contig(name):
        return "sp1_contigs"
    if is_sp2_scaffold(name):
        return "sp2_scaffolds"
    if is_sp2_contig(name):
        return "sp2_contigs"
    return "misc"

# Assign a stable color per contig name
contig_color_map = {}

def color_for(name: str) -> str:
    # Scaffolds: fixed colors if known
    bn = base_name(name)
    if bn in SCAFFOLD_COLORS or name in SCAFFOLD_COLORS:
        return SCAFFOLD_COLORS.get(name, SCAFFOLD_COLORS.get(bn, "#999999"))
    # Contigs: unique color per full contig id
    if is_sp1_contig(name) or is_sp2_contig(name):
        if name not in contig_color_map:
            idx = len(contig_color_map) % CONTIG_MAX_COLORS
            rgba = CONTIG_CMAP(idx / max(1, CONTIG_MAX_COLORS - 1))
            contig_color_map[name] = mpl.colors.to_hex(rgba)
        return contig_color_map[name]
    # Fallback
    return "#999999"

# =========================
# EXTRACT ALIGNED BLOCKS
# =========================
features_by_group = defaultdict(list)

def add_blocks_for(name: str, aln_seq: str):
    start = None
    for i, base in enumerate(aln_seq):
        if base != "-" and ref_coords[i] is not None:
            if start is None:
                start = ref_coords[i]
        else:
            if start is not None:
                end = ref_coords[i - 1] + 1
                features_by_group[group_of(name)].append(
                    GraphicFeature(
                        start=start, end=end, strand=0,
                        color=color_for(name), label=None,  # label suppressed to avoid clutter
                        linecolor=FEATURE_EDGE_COLOR
                    )
                )
                start = None
    if start is not None:
        end = ref_coords[-1] + 1
        features_by_group[group_of(name)].append(
            GraphicFeature(
                start=start, end=end, strand=0,
                color=color_for(name), label=None,
                linecolor=FEATURE_EDGE_COLOR
            )
        )

for name, aln_seq in seqs.items():
    if name == REF_NAME:
        continue
    add_blocks_for(name, aln_seq)

# =========================
# LANE ASSIGNMENT
# =========================
def assign_lanes(features, padding=0):
    feats = sorted(features, key=lambda f: (f.start, f.end))
    lanes = []
    lane_ends = []
    for f in feats:
        placed = False
        for li, last_end in enumerate(lane_ends):
            if f.start >= last_end + padding:
                lanes[li].append(f)
                lane_ends[li] = f.end
                placed = True
                break
        if not placed:
            lanes.append([f])
            lane_ends.append(f.end)
    return lanes

group_order = ["sp1_scaffolds", "sp1_contigs", "sp2_scaffolds", "sp2_contigs"]
if "misc" in features_by_group:
    group_order.append("misc")

lanes_by_group = {}
total_axes = 1  # one for the reference
for g in group_order:
    lanes = assign_lanes(features_by_group.get(g, []), padding=LANE_PADDING)
    lanes_by_group[g] = lanes
    total_axes += max(1, len(lanes))

# =========================
# LEGEND (for contigs)
# =========================
# Build legend entries for contigs only, in display order by group then by name
legend_items = OrderedDict()
for g in group_order:
    for feats in lanes_by_group[g]:
        for feat in feats:
            # We need the name to retrieve the color mapping; rebuild from color map
            # dna_features_viewer loses the original name; instead, reconstruct legend
            # from contig_color_map which we filled during color_for calls:
            pass

# The above "pass" is just a placeholder; we need names. We can collect names earlier.

# Collect names per group for legend (only contigs)
names_per_group = {g: [] for g in group_order}
for name in seqs:
    if name == REF_NAME:
        continue
    g = group_of(name)
    if g in ("sp1_contigs", "sp2_contigs"):
        names_per_group[g].append(name)

# Stable sort
for g in names_per_group:
    names_per_group[g] = sorted(set(names_per_group[g]))

# Build legend entries
legend_labels = []
legend_handles = []
for g in ["sp1_contigs", "sp2_contigs"]:
    if names_per_group[g]:
        # group header as a dummy transparent handle
        legend_labels.append(f"{g}")
        legend_handles.append(mpl.lines.Line2D([], [], linestyle="none"))
        for nm in names_per_group[g]:
            legend_labels.append(nm)
            legend_handles.append(
                mpl.patches.Patch(facecolor=color_for(nm), edgecolor="none")
            )

# =========================
# PLOTTING
# =========================
fig_height = max(8, PER_ROW_HEIGHT * total_axes)
fig, axes = plt.subplots(
    nrows=total_axes, ncols=1,
    figsize=(BASE_FIGURE_WIDTH, fig_height),
    squeeze=False
)
axes = [ax for axrow in axes for ax in axrow]

ax_idx = 0

# Reference on top
ref_record = GraphicRecord(sequence_length=REF_LEN, features=[
    GraphicFeature(start=0, end=REF_LEN, strand=0, color=COLOR_REF, label=None)
])
ref_record.plot(ax=axes[ax_idx], figure_width=BASE_FIGURE_WIDTH)
axes[ax_idx].set_xlim(0, REF_LEN)
axes[ax_idx].set_ylabel("reference")
axes[ax_idx].set_title(TITLE)
axes[ax_idx].set_yticks([])
ax_idx += 1

def group_label(ax, name):
    labels = {
        "sp1_scaffolds": "sp1 scaffolds",
        "sp1_contigs": "sp1 contigs",
        "sp2_scaffolds": "sp2 scaffolds",
        "sp2_contigs": "sp2 contigs",
        "misc": "other"
    }
    ax.set_ylabel(labels.get(name, name))

# Draw groups
for g in group_order:
    lanes = lanes_by_group[g]
    if len(lanes) == 0:
        empty_rec = GraphicRecord(sequence_length=REF_LEN, features=[])
        empty_rec.plot(ax=axes[ax_idx], figure_width=BASE_FIGURE_WIDTH)
        axes[ax_idx].set_xlim(0, REF_LEN)
        axes[ax_idx].set_yticks([])
        group_label(axes[ax_idx], g)
        ax_idx += 1
        continue
    for li, lane in enumerate(lanes):
        rec = GraphicRecord(sequence_length=REF_LEN, features=lane)
        rec.plot(ax=axes[ax_idx], figure_width=BASE_FIGURE_WIDTH)
        axes[ax_idx].set_xlim(0, REF_LEN)
        axes[ax_idx].set_yticks([])
        if li == 0:
            group_label(axes[ax_idx], g)
        ax_idx += 1

axes[-1].set_xlabel("Position on viral genome (bp)")

# Add legend on the right if there are contigs
if legend_labels:
    fig.legend(
        legend_handles, legend_labels,
        loc="center left", bbox_to_anchor=(1.005, 0.5),
        frameon=False, title="Sequences"
    )
    plt.subplots_adjust(right=0.80)  # leave room for legend

plt.tight_layout()

out_png = os.path.join(OUTPUT_DIR, "viral_genome_alignment_stacked_per_contig_color.png")
out_svg = os.path.join(OUTPUT_DIR, "viral_genome_alignment_stacked_per_contig_color.svg")
plt.savefig(out_png, dpi=300)
plt.savefig(out_svg)
plt.close()

# Console summary for sanity checks
print(f"Saved: {out_png}")
print(f"Saved: {out_svg}")
print(f"Reference: {REF_NAME} length {REF_LEN} bp")
for g in group_order:
    total_feats = sum(len(l) for l in lanes_by_group[g])
    print(f"{g}: lanes={len(lanes_by_group[g])}, blocks={total_feats}")
