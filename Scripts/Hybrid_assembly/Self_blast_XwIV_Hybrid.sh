#!/bin/bash

#SBATCH --job-name=XwIV_self_alignment
#SBATCH --output=/beegfs/data/soukkal/Thesis/Hymenoptera_Project/PDV/Results/Viral_genes/IVSPER_Analysis/Long_reads/07_Hybrid_assembly/Self_alignment/Logs/self_%j.out
#SBATCH --error=/beegfs/data/soukkal/Thesis/Hymenoptera_Project/PDV/Results/Viral_genes/IVSPER_Analysis/Long_reads/07_Hybrid_assembly/Self_alignment/Logs/self_%j.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --time=02:00:00
#SBATCH --mem=8G

set -euo pipefail


#################
# Paths         #
#################

BASE="/beegfs/data/soukkal/Thesis/Hymenoptera_Project/PDV/Results/Viral_genes/IVSPER_Analysis"

ASSEMBLY="${BASE}/Long_reads/07_Hybrid_assembly/MaSuRCA/CA.mr.99.17.15.0.02/primary.genome.scf.fasta"

OUTDIR="${BASE}/Long_reads/07_Hybrid_assembly/Self_alignment"

LOGDIR="${OUTDIR}/Logs"

MMSEQS_DIR="${OUTDIR}/MMseqs"

TMPDIR="${MMSEQS_DIR}/tmp"

VIRAL_SCAFFOLD="jcf7180000126179"


#################
# Software      #
#################

MMSEQS="/beegfs/data/soukkal/TOOLS/mmseqs/bin/mmseqs"

SAMTOOLS="/beegfs/data/soft/samtools-1.19/bin/samtools"

THREADS=8


#################
# Output files  #
#################

VIRAL_FASTA="${OUTDIR}/XwIV_hybrid_scaffold.fa"

ALL_HITS="${OUTDIR}/XwIV_self_alignment_ALL.tsv"

FILTERED_HITS="${OUTDIR}/XwIV_self_alignment_filtered.tsv"


#################
# Directories   #
#################

mkdir -p "${OUTDIR}"
mkdir -p "${LOGDIR}"
mkdir -p "${MMSEQS_DIR}"
mkdir -p "${TMPDIR}"


#########################################
# 1. Extract XwIV scaffold              #
#########################################

echo "Extracting XwIV scaffold..."

${SAMTOOLS} faidx "${ASSEMBLY}"

${SAMTOOLS} faidx \
    "${ASSEMBLY}" \
    "${VIRAL_SCAFFOLD}" \
    > "${VIRAL_FASTA}"


#########################################
# 2. Create MMseqs nucleotide database  #
#########################################

echo "Creating MMseqs database..."

${MMSEQS} createdb \
    "${VIRAL_FASTA}" \
    "${MMSEQS_DIR}/XwIV_DB"


#########################################
# 3. Self-search                        #
#########################################

echo "Running nucleotide self-search..."

${MMSEQS} search \
    "${MMSEQS_DIR}/XwIV_DB" \
    "${MMSEQS_DIR}/XwIV_DB" \
    "${MMSEQS_DIR}/XwIV_self_result" \
    "${TMPDIR}" \
    --search-type 3 \
    -e 0.001 \
    --threads ${THREADS}


#########################################
# 4. Export alignments                  #
#########################################

${MMSEQS} convertalis \
    "${MMSEQS_DIR}/XwIV_DB" \
    "${MMSEQS_DIR}/XwIV_DB" \
    "${MMSEQS_DIR}/XwIV_self_result" \
    "${ALL_HITS}" \
    --format-output "query,target,qstart,qend,tstart,tend,alnlen,pident,evalue,bits"


##############################################
# 5. Remove trivial self-alignment           #
##############################################

python3 - <<EOF

import pandas as pd

input_file = "${ALL_HITS}"
output_file = "${FILTERED_HITS}"

columns = [
    "query",
    "target",
    "query_start",
    "query_end",
    "target_start",
    "target_end",
    "alignment_length",
    "identity_percent",
    "evalue",
    "bitscore"
]

df = pd.read_csv(
    input_file,
    sep="\t",
    names=columns
)

print(f"All alignments: {len(df)}")

# Remove only alignments in which a region is aligned
# to exactly the same coordinates on the scaffold.
#
# Non-trivial alignments between different regions of
# the same scaffold are retained.

same_region = (
    (df["query"] == df["target"]) &
    (df["query_start"] == df["target_start"]) &
    (df["query_end"] == df["target_end"])
)

filtered = df[~same_region].copy()

# E-value threshold, retained explicitly here even
# though MMseqs search already used -e 0.001.

filtered = filtered[
    filtered["evalue"] <= 0.001
].copy()

# Sort strongest alignments first.

filtered = filtered.sort_values(
    ["evalue", "bitscore"],
    ascending=[True, False]
)

filtered.to_csv(
    output_file,
    sep="\t",
    index=False
)

print(f"Non-trivial homologous alignments: {len(filtered)}")
print(f"Output: {output_file}")

EOF


#################
# Done          #
#################

echo ""
echo "Self-alignment completed."
echo ""
echo "All hits:"
echo "${ALL_HITS}"
echo ""
echo "Non-trivial homologous regions:"
echo "${FILTERED_HITS}"
