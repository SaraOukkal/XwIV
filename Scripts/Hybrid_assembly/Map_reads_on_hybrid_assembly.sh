#!/bin/bash

#SBATCH --job-name=map_Xwirra_hybrid
#SBATCH --output=/beegfs/data/soukkal/Thesis/Hymenoptera_Project/PDV/Results/Viral_genes/IVSPER_Analysis/Long_reads/07_Hybrid_assembly/Mapping/Logs/map_%j.out
#SBATCH --error=/beegfs/data/soukkal/Thesis/Hymenoptera_Project/PDV/Results/Viral_genes/IVSPER_Analysis/Long_reads/07_Hybrid_assembly/Mapping/Logs/map_%j.err
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=30:00:00
#SBATCH --constraint='skylake|haswell|broadwell'
#SBATCH --exclude=pbil-deb27


set -euo pipefail


#################
# Paths         #
#################

BASE="/beegfs/data/soukkal/Thesis/Hymenoptera_Project/PDV/Results/Viral_genes/IVSPER_Analysis/Long_reads"

ASSEMBLY="${BASE}/07_Hybrid_assembly/MaSuRCA/CA.mr.99.17.15.0.02/primary.genome.scf.fasta"

ONT_READS="${BASE}/01_Reads/Xiphosomella_wirra_ONT.fastq.gz"

ILLUMINA_R1="${BASE}/07_Hybrid_assembly/Reads/Xwirra_combined_R1.fastq.gz"
ILLUMINA_R2="${BASE}/07_Hybrid_assembly/Reads/Xwirra_combined_R2.fastq.gz"

OUTDIR="${BASE}/07_Hybrid_assembly/Mapping"

ONT_OUTDIR="${OUTDIR}/ONT"
ILLUMINA_OUTDIR="${OUTDIR}/Illumina"
LOGDIR="${OUTDIR}/Logs"


#################
# Software      #
#################

minimap2="/beegfs/home/soukkal/miniconda3/bin/minimap2"

samtools="/beegfs/data/soft/samtools-1.19/bin/samtools"

THREADS=8


#################
# Output files  #
#################

# ONT

ONT_BAM="${ONT_OUTDIR}/Xwirra_ONT_vs_hybrid.sorted.bam"

ONT_COVERAGE="${ONT_OUTDIR}/Xwirra_ONT.samtools_coverage.tsv"


# Illumina

ILLUMINA_BAM="${ILLUMINA_OUTDIR}/Xwirra_Illumina_vs_hybrid.sorted.bam"

ILLUMINA_COVERAGE="${ILLUMINA_OUTDIR}/Xwirra_Illumina.samtools_coverage.tsv"


#################
# Directories   #
#################

mkdir -p "${ONT_OUTDIR}"
mkdir -p "${ILLUMINA_OUTDIR}"
mkdir -p "${LOGDIR}"


#####################################
# 1. Mapping ONT -> hybrid assembly #
#####################################

echo "Mapping ONT reads to hybrid assembly..."

${minimap2} \
    -t ${THREADS} \
    -ax map-ont \
    "${ASSEMBLY}" \
    "${ONT_READS}" \
    | ${samtools} view -@ ${THREADS} -b \
    | ${samtools} sort -@ ${THREADS} -o "${ONT_BAM}"

${samtools} index "${ONT_BAM}"


#######################################
# 2. samtools coverage - ONT          #
#######################################

echo "Calculating ONT coverage..."

${samtools} coverage "${ONT_BAM}" > "${ONT_COVERAGE}"


##########################################
# 3. Mapping Illumina -> hybrid assembly #
##########################################

echo "Mapping Illumina reads to hybrid assembly..."

${minimap2} \
    -t ${THREADS} \
    -ax sr \
    "${ASSEMBLY}" \
    "${ILLUMINA_R1}" \
    "${ILLUMINA_R2}" \
    | ${samtools} view -@ ${THREADS} -b \
    | ${samtools} sort -@ ${THREADS} -o "${ILLUMINA_BAM}"

${samtools} index "${ILLUMINA_BAM}"


########################################
# 4. samtools coverage - Illumina      #
########################################

echo "Calculating Illumina coverage..."

${samtools} coverage "${ILLUMINA_BAM}" > "${ILLUMINA_COVERAGE}"


###############################
# Done                        #
###############################

echo ""
echo "Mapping completed."
echo ""
echo "ONT:"
echo "BAM:      ${ONT_BAM}"
echo "Coverage: ${ONT_COVERAGE}"
echo ""
echo "Illumina:"
echo "BAM:      ${ILLUMINA_BAM}"
echo "Coverage: ${ILLUMINA_COVERAGE}"
