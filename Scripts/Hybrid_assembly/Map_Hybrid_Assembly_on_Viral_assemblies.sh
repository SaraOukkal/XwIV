#!/bin/bash

#SBATCH --job-name=XwIV_vs_Hybrid
#SBATCH --cpus-per-task=8
#SBATCH --mem=8G
#SBATCH --time=00:30:00

set -euo pipefail

# =========================
# Paths
# =========================

VIRUS="/beegfs/data/soukkal/Thesis/Hymenoptera_Project/PDV/Results/Viral_genes/IVSPER_Analysis/Viral_Genome/Viral_genome.fa"

HYBRID="/beegfs/data/soukkal/Thesis/Hymenoptera_Project/PDV/Results/Viral_genes/IVSPER_Analysis/Long_reads/07_Hybrid_assembly/MaSuRCA/CA.mr.99.17.15.0.02/primary.genome.scf.fasta"

OUTDIR="/beegfs/data/soukkal/Thesis/Hymenoptera_Project/PDV/Results/Viral_genes/IVSPER_Analysis/Long_reads/07_Hybrid_assembly/Viral_search/Mapping_Hybrid_on_XwIV"

mkdir -p "${OUTDIR}"

# Redirect SLURM logs to output directory
exec > "${OUTDIR}/XwIV_vs_Hybrid_${SLURM_JOB_ID}.out" \
     2> "${OUTDIR}/XwIV_vs_Hybrid_${SLURM_JOB_ID}.err"

OUT="${OUTDIR}/XwIV_vs_Hybrid_assembly.paf"

# =========================
# Alignment
# =========================

echo "Starting minimap2..."
echo "Virus reference:  ${VIRUS}"
echo "Hybrid assembly:  ${HYBRID}"
echo "Output:           ${OUT}"
echo

minimap2 \
    -t "${SLURM_CPUS_PER_TASK}" \
    -x asm20 \
    -c \
    --cs \
    "${VIRUS}" \
    "${HYBRID}" \
    > "${OUT}"

# =========================
# Summary
# =========================

echo
echo "Alignment finished."

echo "Number of alignments:"
wc -l "${OUT}"

echo
echo "20 longest alignments:"
sort -k11,11nr "${OUT}" | head -20

echo
echo "Done."
