#!/bin/bash
#SBATCH --job-name=orf_scaffold3
#SBATCH --output=orf_scaffold3.out
#SBATCH --error=orf_scaffold3.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --time=01:00:00
#SBATCH --mem=1G

set -euo pipefail

# === CONFIGURATION ===
OUTDIR="/beegfs/data/soukkal/Thesis/Hymenoptera_Project/PDV/Results/Viral_genes/IVSPER_Analysis/xiphosomella_wirra_specimen2/ORF_prediction"
INPUT_FASTA="/beegfs/data/soukkal/Thesis/Hymenoptera_Project/PDV/Results/Viral_genes/IVSPER_Analysis/xiphosomella_wirra_specimen2/scaffold3.fa"
CLEAN_FASTA="$OUTDIR/scaffold3_clean.fa"
ORF_FASTA="$OUTDIR/scaffold3_ORFs_70AA.faa"
ORF_TABLE="$OUTDIR/scaffold3_ORFs_70AA.tsv"
GENE_STATS="$OUTDIR/scaffold3_gene_stats.tsv"

# === 1. Nettoyage du header ===
sed 's/|/-/g' "$INPUT_FASTA" > "$CLEAN_FASTA"

# === 2. ORFfinder local (min 70 AA = 210 nt) ===
ORFfinder -in "$CLEAN_FASTA" -out "$ORF_FASTA" -ml 210 -outfmt 0

# === 3. Parsing FASTA -> TSV + calcul stats ===
python3 - <<EOF
import sys
from Bio import SeqIO
import statistics

fasta_path = "$ORF_FASTA"
scaffold_fasta = "$CLEAN_FASTA"
ts_table = "$ORF_TABLE"
out_stats = "$GENE_STATS"

# Parse ORFs
orfs = []
with open(ts_table, "w") as out:
    out.write("scaffold\torf_number\tstart\tend\n")
    for record in SeqIO.parse(fasta_path, "fasta"):
        header = record.description
        try:
            parts = header.split()
            main = parts[0]  # lcl|ORF2_scaffold3-size100523:1738:1974
            main = main.replace("lcl|", "")
            orf_and_scaffold, coords = main.split(":")[0], main.split(":")[1:]
            start, end = map(int, coords)
            orf_tag, scaffold = orf_and_scaffold.split("_", 1)
            orf_number = int(orf_tag.replace("ORF", ""))
            scaffold = scaffold.replace("-", "|")
            orfs.append((scaffold, orf_number, start, end))
            out.write(f"{scaffold}\t{orf_number}\t{start}\t{end}\n")
        except Exception as e:
            print(f"Erreur de parsing sur l'en-t\u00eate : {header} -> {e}")

# Load scaffold length
for record in SeqIO.parse(scaffold_fasta, "fasta"):
    scaffold_length = len(record.seq)
    break

# Générer le masque codant
mask = [0] * scaffold_length
for _, _, start, end in orfs:
    for i in range(min(start, end)-1, max(start, end)):
        if 0 <= i < scaffold_length:
            mask[i] = 1
coding_bp = sum(mask)
coding_pct = coding_bp / scaffold_length * 100

# Distances entre ORFs consécutifs
orfs_sorted = sorted(orfs, key=lambda x: min(x[2], x[3]))
distances = []
for i in range(1, len(orfs_sorted)):
    prev_end = max(orfs_sorted[i-1][2], orfs_sorted[i-1][3])
    curr_start = min(orfs_sorted[i][2], orfs_sorted[i][3])
    distances.append(max(0, curr_start - prev_end))

# Nombre d'ORFs > 100 AA (300 nt)
orfs_over_100aa = [1 for _, _, s, e in orfs if abs(e - s) >= 300]

# Sauvegarde des stats
with open(out_stats, "w") as out:
    out.write("scaffold\ttotal_bp\tcoding_bp\tcoding_pct\tn_ORFs\tn_ORFs_>100AA\tdist_min\tdist_mean\tdist_max\n")
    out.write(f"scaffold3\t{scaffold_length}\t{coding_bp}\t{coding_pct:.2f}\t{len(orfs)}\t{len(orfs_over_100aa)}\t")
    if distances:
        out.write(f"{min(distances)}\t{statistics.mean(distances):.2f}\t{max(distances)}\n")
    else:
        out.write("NA\tNA\tNA\n")
EOF

