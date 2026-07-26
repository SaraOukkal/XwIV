#!/usr/bin/env python3

# =========================================================
# GC and Depth distributions with BUSCO reference, per specimen
# - Exports SVG plot and TSV summary table
# - Adds 5-95% ranges and BUSCO medians for GC and Depth
# - BUSCO histograms are gray: edge alpha 1.0, fill alpha 0.5
# Comments in English. No emojis.
# =========================================================

import os
import gzip
import pandas as pd
import numpy as np
from Bio import SeqIO
import matplotlib.pyplot as plt
import seaborn as sns

# -----------------------
# Config
# -----------------------
output_dir = "/beegfs/data/soukkal/Thesis/Hymenoptera_Project/PDV/Results/Viral_genes/IVSPER_Analysis/"
output_plot = os.path.join(output_dir, "GC_and_Depth_Distributions_All_BUSCO.svg")
output_table = os.path.join(output_dir, "Test_depth_GC_all_BUSCO.tsv")

specimens = {
    "sp1": {
        "depth_bed": "/beegfs/data/soukkal/Thesis/Hymenoptera_Project/PDV/Results/Viral_genes/IVSPER_Analysis/xiphosomella_wirra_specimen1/Mapping/xiphosomella_wirra:DHJPAR0041296.per-base.bed",
        "busco_file": "/beegfs/project/horizon/data/stats/busco/specimens/DHJPAR0041296/run_insecta_odb10/full_table.tsv",
        "assembly_fasta": "/beegfs/project/horizon/data/assembly/specimens/DHJPAR0041296/redundans/scaffolds.reduced.fa",
        "candidates": [
            "scaffold524|size39189",
            "scaffold187|size52917",
            "scaffold33401|size845"
        ]
    },
    "sp2": {
        "depth_bed": "/beegfs/data/soukkal/Thesis/Hymenoptera_Project/PDV/Results/Viral_genes/IVSPER_Analysis/xiphosomella_wirra_specimen2/Mapping/bed/DHJPAR0036296.per-base.bed",
        "busco_file": "/beegfs/project/horizon/data/stats/busco/specimens/DHJPAR0036296/run_insecta_odb10/full_table.tsv",
        "assembly_fasta": "/beegfs/project/horizon/data/assembly/specimens/DHJPAR0036296/redundans/scaffolds.reduced.fa",
        "candidates": [
            "scaffold3|size100523"
        ]
    }
}

# -----------------------
# Helpers
# -----------------------
def load_busco_scaffolds(path):
    """Return the set of scaffolds that contain BUSCOs from BUSCO full_table.tsv."""
    scaffolds = set()
    with open(path) as f:
        for line in f:
            if line.startswith("#") or "Missing" in line:
                continue
            parts = line.strip().split("\t")
            if len(parts) >= 3:
                scaffolds.add(parts[2])
    return scaffolds

def read_depth_median(path, scaffolds, cache_path):
    """Compute or load per-scaffold median depth from a per-base BED. Only for specified scaffolds."""
    if os.path.exists(cache_path):
        print(f"  Loading cached depth values from {cache_path}")
        return pd.read_csv(cache_path, sep="\t", index_col=0).to_dict()["median_depth"]

    print(f"  Calculating depth values from {path}")
    open_fn = gzip.open if path.endswith(".gz") else open
    depth_dict = {}
    with open_fn(path, "rt") as f:
        for line in f:
            parts = line.strip().split("\t")
            if len(parts) < 3:
                continue
            scaf, _, depth = parts
            if scaf not in scaffolds:
                continue
            try:
                depth = float(depth)
            except ValueError:
                continue
            depth_dict.setdefault(scaf, []).append(depth)

    medians = {k: np.median(v) for k, v in depth_dict.items() if len(v) > 0}
    pd.DataFrame.from_dict(medians, orient="index", columns=["median_depth"]).to_csv(cache_path, sep="\t")
    return medians

def read_gc(path, cache_path):
    """Compute or load GC content per scaffold from FASTA."""
    if os.path.exists(cache_path):
        print(f"  Loading cached GC values from {cache_path}")
        return pd.read_csv(cache_path, sep="\t", index_col=0).to_dict()["gc"]

    print(f"  Calculating GC content from {path}")
    gc_dict = {}
    for record in SeqIO.parse(path, "fasta"):
        seq = record.seq.upper()
        length = len(seq)
        if length == 0:
            gc = 0.0
        else:
            gc = (seq.count("G") + seq.count("C")) / length * 100.0
        gc_dict[record.id] = gc

    pd.DataFrame.from_dict(gc_dict, orient="index", columns=["gc"]).to_csv(cache_path, sep="\t")
    return gc_dict

def empirical_pval(dist, val):
    """Two-sided empirical p-value around the median using absolute deviation."""
    if not dist or val is None:
        return "NA"
    median = np.median(dist)
    delta = abs(val - median)
    extreme = [x for x in dist if abs(x - median) >= delta]
    return round((len(extreme) + 1) / (len(dist) + 1), 4)

# -----------------------
# Plot setup
# -----------------------
sns.set_style("whitegrid")
fig, axs = plt.subplots(2, 2, figsize=(12, 8))
results = []

# Gray colors with explicit alpha
gray_fill = (0.5, 0.5, 0.5, 0.5)   # facecolor alpha 0.5
gray_edge = (0.5, 0.5, 0.5, 1.0)   # edgecolor alpha 1.0

# -----------------------
# Main loop
# -----------------------
for i, (sp, data) in enumerate(specimens.items()):
    print(f"\nProcessing {sp}...")

    busco_scaffs = load_busco_scaffolds(data["busco_file"])
    candidate_scaffs = list(data["candidates"])  # keep order
    all_needed = busco_scaffs | set(candidate_scaffs)

    gc_dict = read_gc(data["assembly_fasta"], os.path.join(output_dir, f"{sp}_gc.tsv"))
    depth_by_scaffold = read_depth_median(data["depth_bed"], all_needed, os.path.join(output_dir, f"{sp}_depth.tsv"))

    gc_busco = [gc_dict[s] for s in busco_scaffs if s in gc_dict]
    depth_busco = [depth_by_scaffold[s] for s in busco_scaffs if s in depth_by_scaffold]

    # Medians and quantiles; handle empty lists robustly
    median_gc_busco = np.median(gc_busco) if len(gc_busco) > 0 else None
    median_depth_busco = np.median(depth_busco) if len(depth_busco) > 0 else None
    q_gc = np.quantile(gc_busco, [0.05, 0.95]) if len(gc_busco) > 0 else [None, None]
    q_depth = np.quantile(depth_busco, [0.05, 0.95]) if len(depth_busco) > 0 else [None, None]

    # Palette for candidate vertical lines
    color_palette = sns.color_palette("husl", len(candidate_scaffs)) if len(candidate_scaffs) > 0 else []

    # --- GC histogram ---
    ax_gc = axs[i, 0]
    sns.histplot(gc_busco, bins=30, kde=False, ax=ax_gc,
                 color=gray_fill, edgecolor=gray_edge, linewidth=0.8)
    # Candidate lines
    for idx, c in enumerate(candidate_scaffs):
        val = gc_dict.get(c, None)
        if val is not None:
            ax_gc.axvline(val, linestyle="--", label=c, color=color_palette[idx])
    # 5-95% range lines
    if None not in q_gc:
        ax_gc.axvline(q_gc[0], color="black", linestyle=":")
        ax_gc.axvline(q_gc[1], color="black", linestyle=":")
    ax_gc.set_title(f"{sp} - GC")
    ax_gc.set_xlabel("GC content (%)")
    ax_gc.legend(fontsize="small")

    # --- Depth histogram ---
    ax_d = axs[i, 1]
    sns.histplot(depth_busco, bins=30, kde=False, ax=ax_d,
                 color=gray_fill, edgecolor=gray_edge, linewidth=0.8)
    for idx, c in enumerate(candidate_scaffs):
        val = depth_by_scaffold.get(c, None)
        if val is not None:
            ax_d.axvline(val, linestyle="--", label=c, color=color_palette[idx])
    if None not in q_depth:
        ax_d.axvline(q_depth[0], color="black", linestyle=":")
        ax_d.axvline(q_depth[1], color="black", linestyle=":")
    ax_d.set_title(f"{sp} - Depth")
    ax_d.set_xlabel("Median depth")
    ax_d.legend(fontsize="small")

    # --- Collect per-candidate stats ---
    for c in candidate_scaffs:
        gc_val = gc_dict.get(c, None)
        depth_val = depth_by_scaffold.get(c, None)

        p_gc = empirical_pval(gc_busco, gc_val)
        p_depth = empirical_pval(depth_busco, depth_val)

        test_stat_gc = (
            f"{round(q_gc[0], 2)}-{round(q_gc[1], 2)}" if None not in q_gc else "NA"
        )
        test_stat_depth = (
            f"{round(q_depth[0], 2)}-{round(q_depth[1], 2)}" if None not in q_depth else "NA"
        )

        nuclear_gc = (
            "Yes" if (gc_val is not None and None not in q_gc and q_gc[0] <= gc_val <= q_gc[1])
            else ("No" if gc_val is not None else "NA")
        )
        nuclear_depth = (
            "Yes" if (depth_val is not None and None not in q_depth and q_depth[0] <= depth_val <= q_depth[1])
            else ("No" if depth_val is not None else "NA")
        )

        results.append({
            "specimen": sp,
            "candidate": c,
            "GC_candidate": round(gc_val, 2) if gc_val is not None else "NA",
            "pval_gc": p_gc,
            "Test_Stat_gc": test_stat_gc,
            "Median_BUSCO_GC": round(median_gc_busco, 2) if median_gc_busco is not None else "NA",
            "Nuclear_GC": nuclear_gc,
            "Median_Depth": round(depth_val, 2) if depth_val is not None else "NA",
            "pval_depth": p_depth,
            "Test_Stat_depth": test_stat_depth,
            "Median_BUSCO_Depth": round(median_depth_busco, 2) if median_depth_busco is not None else "NA",
            "Nuclear_depth": nuclear_depth,
            "BUSCO_count": len(busco_scaffs)
        })

# -----------------------
# Save plot and table
# -----------------------
plt.tight_layout()
fig.savefig(output_plot)
print(f"\nSaved plot to: {output_plot}")

df = pd.DataFrame(results)
df.to_csv(output_table, sep="\t", index=False)
print(f"Saved summary table to: {output_table}")

