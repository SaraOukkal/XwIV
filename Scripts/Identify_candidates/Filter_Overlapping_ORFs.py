#!/usr/bin/env python3

import pandas as pd
from Bio import SeqIO
import os
import re

# Input paths
tsv_path = "/beegfs/data/soukkal/Thesis/Hymenoptera_Project/PDV/Results/Viral_genes/IVSPER_Analysis/xiphosomella_wirra_specimen2/ORF_prediction/scaffold3_ORFs_70AA.tsv"
fasta_path = "/beegfs/data/soukkal/Thesis/Hymenoptera_Project/PDV/Results/Viral_genes/IVSPER_Analysis/xiphosomella_wirra_specimen2/ORF_prediction/scaffold3_ORFs_70AA.faa"
annotated_path = "/beegfs/data/soukkal/Thesis/Hymenoptera_Project/PDV/Results/Viral_genes/IVSPER_Analysis/xiphosomella_wirra_specimen2/ORF_prediction/ORFs_annotated.txt"

# Output paths
output_dir = os.path.dirname(tsv_path)
filtered_tsv = os.path.join(output_dir, "scaffold3_ORFs_70AA_filtered.tsv")
filtered_faa = os.path.join(output_dir, "scaffold3_ORFs_70AA_filtered.faa")

# Load ORF table
orf_df = pd.read_csv(tsv_path, sep="\t")

# Compute strand and length
orf_df["strand"] = orf_df.apply(lambda row: "+" if row["end"] >= row["start"] else "-", axis=1)
orf_df["length"] = abs(orf_df["end"] - orf_df["start"]) + 1

# Coordinates normalisées pour les comparaisons seulement
orf_df["start_norm"] = orf_df[["start", "end"]].min(axis=1)
orf_df["end_norm"] = orf_df[["start", "end"]].max(axis=1)

# Liste des ORFs annotés
with open(annotated_path) as f:
    annotated = set(f.read().split())

# Fonction pour détecter si row_i est contenu (>50%) dans row_j
def is_contained(row_i, row_j):
    start_i, end_i = row_i["start_norm"], row_i["end_norm"]
    start_j, end_j = row_j["start_norm"], row_j["end_norm"]
    overlap = max(0, min(end_i, end_j) - max(start_i, start_j) + 1)
    return (overlap / row_i["length"]) > 0.5

# Application du filtre (tous brins confondus)
keep_orfs = set()
for idx_i, row_i in orf_df.iterrows():
    keep = True
    for idx_j, row_j in orf_df.iterrows():
        if idx_i == idx_j:
            continue
        if is_contained(row_i, row_j):
            if row_i["length"] < row_j["length"] and str(row_i["orf_number"]) not in annotated:
                keep = False
                break
    if keep:
        keep_orfs.add(row_i["orf_number"])

# Génération du TSV filtré (avec start/end d'origine + strand)
filtered_df = orf_df[orf_df["orf_number"].isin(keep_orfs)].copy()
filtered_df = filtered_df[["scaffold", "orf_number", "start", "end", "strand"]]
filtered_df.to_csv(filtered_tsv, sep="\t", index=False)

# Extraction du numéro d’ORF à partir du header FASTA
def extract_orf_number(header):
    match = re.search(r'ORF(\d+)_', header)
    return int(match.group(1)) if match else None

# FASTA filtré
records = list(SeqIO.parse(fasta_path, "fasta"))
filtered_records = [r for r in records if extract_orf_number(r.id) in keep_orfs]
SeqIO.write(filtered_records, filtered_faa, "fasta")

# Résumé
print(f"Initial ORFs: {len(orf_df)}")
print(f"Filtered ORFs: {len(filtered_df)}")
print(f"Removed ORFs: {len(orf_df) - len(filtered_df)}")
print(f"Filtered TSV written to: {filtered_tsv}")
print(f"Filtered FASTA written to: {filtered_faa}")

