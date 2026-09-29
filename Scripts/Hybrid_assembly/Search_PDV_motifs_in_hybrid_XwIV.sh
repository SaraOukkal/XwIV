#!/bin/bash

#SBATCH --job-name=XwIV_motifs_Glypta
#SBATCH --output=%x_%j.out
#SBATCH --error=%x_%j.err
#SBATCH --time=01:00:00
#SBATCH --cpus-per-task=2
#SBATCH --mem=10G
#SBATCH --exclude=pbil-deb[14-27]

set -euo pipefail


# ============================================================
# PATHS
# ============================================================

MMSEQS="/beegfs/data/soukkal/TOOLS/mmseqs/bin/mmseqs"
SAMTOOLS="/beegfs/data/soft/samtools-1.19/bin/samtools"

ASSEMBLY="/beegfs/data/soukkal/Thesis/Hymenoptera_Project/PDV/Results/Viral_genes/IVSPER_Analysis/Long_reads/07_Hybrid_assembly/MaSuRCA/CA.mr.99.17.15.0.02/primary.genome.scf.fasta"

CONTIG="jcf7180000126179"

MOTIFS="/beegfs/data/soukkal/Thesis/Hymenoptera_Project/PDV/Data/All_motifs.fa"

GLYPTA_SEGMENTS="/ceph/projets/horizon/PDV_Segments/Data/Segments/Glypta_fumiferanae_segments.fa"

OUTPUT_DIR="/beegfs/data/soukkal/Thesis/Hymenoptera_Project/PDV/Results/Viral_genes/IVSPER_Analysis/Long_reads/07_Hybrid_assembly/Motif_search"

VIRAL_GENOME="${OUTPUT_DIR}/XwIV_hybrid_genome.fa"

MOTIF_DIR="${OUTPUT_DIR}/PDV_motifs"
GLYPTA_DIR="${OUTPUT_DIR}/Glypta_segments"


# ============================================================
# CREATE OUTPUT DIRECTORIES
# ============================================================

mkdir -p "${OUTPUT_DIR}"
mkdir -p "${MOTIF_DIR}"
mkdir -p "${GLYPTA_DIR}"


# ============================================================
# CHECK INPUT FILES
# ============================================================

for FILE in "${ASSEMBLY}" "${MOTIFS}" "${GLYPTA_SEGMENTS}"; do
    if [[ ! -s "${FILE}" ]]; then
        echo "ERROR: input file not found or empty:"
        echo "${FILE}"
        exit 1
    fi
done


# ============================================================
# EXTRACT NEW XwIV GENOME
# ============================================================

echo
echo "============================================"
echo "Extracting new XwIV genome"
echo "============================================"

"${SAMTOOLS}" faidx "${ASSEMBLY}"

"${SAMTOOLS}" faidx \
    "${ASSEMBLY}" \
    "${CONTIG}" \
    > "${VIRAL_GENOME}"

echo "XwIV genome:"
echo "${VIRAL_GENOME}"


# ============================================================
# CREATE XwIV DATABASE
# ============================================================

echo
echo "============================================"
echo "Creating XwIV MMseqs database"
echo "============================================"

"${MMSEQS}" createdb \
    "${VIRAL_GENOME}" \
    "${OUTPUT_DIR}/XwIV_DB"


# ============================================================
# SEARCH 1
# XwIV vs PDV MOTIFS
# ============================================================

echo
echo "============================================"
echo "XwIV vs PDV motifs"
echo "============================================"

"${MMSEQS}" createdb \
    "${MOTIFS}" \
    "${MOTIF_DIR}/motifs_DB"


"${MMSEQS}" search \
    "${OUTPUT_DIR}/XwIV_DB" \
    "${MOTIF_DIR}/motifs_DB" \
    "${MOTIF_DIR}/result" \
    "${MOTIF_DIR}/tmp" \
    --search-type 2 \
    -e 0.01 \
    -s 7.5 \
    -a \
    --threads 2 \
    --remove-tmp-files


"${MMSEQS}" convertalis \
    "${OUTPUT_DIR}/XwIV_DB" \
    "${MOTIF_DIR}/motifs_DB" \
    "${MOTIF_DIR}/result" \
    "${MOTIF_DIR}/XwIV_vs_PDV_motifs.m8" \
    --format-output 'query,qlen,tlen,target,pident,alnlen,mismatch,gapopen,qstart,qend,tstart,tend,evalue,bits,qaln,tcov' \
    --search-type 2


# ============================================================
# SEARCH 2
# XwIV vs GLYPTA FUMIFERANAE SEGMENTS
# ============================================================

echo
echo "============================================"
echo "XwIV vs Glypta fumiferanae segments"
echo "============================================"

"${MMSEQS}" createdb \
    "${GLYPTA_SEGMENTS}" \
    "${GLYPTA_DIR}/glypta_DB"


"${MMSEQS}" search \
    "${OUTPUT_DIR}/XwIV_DB" \
    "${GLYPTA_DIR}/glypta_DB" \
    "${GLYPTA_DIR}/result" \
    "${GLYPTA_DIR}/tmp" \
    --search-type 2 \
    -e 0.01 \
    -s 7.5 \
    -a \
    --threads 2 \
    --remove-tmp-files


"${MMSEQS}" convertalis \
    "${OUTPUT_DIR}/XwIV_DB" \
    "${GLYPTA_DIR}/glypta_DB" \
    "${GLYPTA_DIR}/result" \
    "${GLYPTA_DIR}/XwIV_vs_Glypta_segments.m8" \
    --format-output 'query,qlen,tlen,target,pident,alnlen,mismatch,gapopen,qstart,qend,tstart,tend,evalue,bits,qaln,tcov' \
    --search-type 2


# ============================================================
# DONE
# ============================================================

echo
echo "============================================"
echo "ANALYSIS COMPLETE"
echo "============================================"

echo
echo "XwIV genome:"
echo "${VIRAL_GENOME}"

echo
echo "PDV motif results:"
echo "${MOTIF_DIR}/XwIV_vs_PDV_motifs.m8"

echo
echo "Glypta segment results:"
echo "${GLYPTA_DIR}/XwIV_vs_Glypta_segments.m8"

echo
echo "============================================"
