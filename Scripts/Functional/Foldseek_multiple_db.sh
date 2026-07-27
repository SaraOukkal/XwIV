#!/bin/bash
#SBATCH --job-name=foldseek_search
#SBATCH --cpus-per-task=8
#SBATCH --exclude=pbil-deb[19-27]
#SBATCH --output=/beegfs/data/soukkal/Thesis/Hymenoptera_Project/PDV/Results/Viral_genes/IVSPER_Analysis/Foldseek/foldseek_%j.out
#SBATCH --error=/beegfs/data/soukkal/Thesis/Hymenoptera_Project/PDV/Results/Viral_genes/IVSPER_Analysis/Foldseek/foldseek_%j.err

# Exit on error
set -euo pipefail

# Paths
QUERY="/beegfs/data/soukkal/Thesis/Hymenoptera_Project/PDV/Results/Viral_genes/IVSPER_Analysis/xiphosomella_wirra_specimen2/ORF_prediction/scaffold3_ORFs_70AA_filtered.faa"

DB_SWISSPROT="/beegfs/home/soukkal/Tools/swissprot_db/swissprot"
DB_UNIPROT50="/beegfs/data/soukkal/Thesis/Hymenoptera_Project/PDV/Results/Viral_genes/IVSPER_Analysis/Foldseek/uniprot50"
DB_BFVD="/beegfs/data/soukkal/Thesis/Hymenoptera_Project/PDV/Results/Viral_genes/IVSPER_Analysis/Foldseek/bfvd"

OUTDIR="/beegfs/data/soukkal/Thesis/Hymenoptera_Project/PDV/Results/Viral_genes/IVSPER_Analysis/Foldseek"

mkdir -p ${OUTDIR}

################################
# SwissProt search
################################

OUT_SWISSPROT="${OUTDIR}/swissprot"
TMP_SWISSPROT="${OUT_SWISSPROT}/tmp"
RESULT_SWISSPROT="${OUT_SWISSPROT}/scaffold3_foldseek_swissprot.m8"

mkdir -p ${OUT_SWISSPROT}
mkdir -p ${TMP_SWISSPROT}

foldseek easy-search \
${QUERY} \
${DB_SWISSPROT} \
${RESULT_SWISSPROT} \
${TMP_SWISSPROT} \
--threads 8 \
-e 0.01 \
--format-output query,target,alnlen,qstart,qend,tstart,tend,evalue,bits,tcov,qcov

################################
# UniProt50 search
################################

OUT_UNIPROT50="${OUTDIR}/uniprot50"
TMP_UNIPROT50="${OUT_UNIPROT50}/tmp"
RESULT_UNIPROT50="${OUT_UNIPROT50}/scaffold3_foldseek_uniprot50.m8"

mkdir -p ${OUT_UNIPROT50}
mkdir -p ${TMP_UNIPROT50}

foldseek easy-search \
${QUERY} \
${DB_UNIPROT50} \
${RESULT_UNIPROT50} \
${TMP_UNIPROT50} \
--threads 8 \
-e 0.01 \
--format-output query,target,alnlen,qstart,qend,tstart,tend,evalue,bits,tcov,qcov

################################
# BFVD search
################################

OUT_BFVD="${OUTDIR}/bfvd"
TMP_BFVD="${OUT_BFVD}/tmp"
RESULT_BFVD="${OUT_BFVD}/scaffold3_foldseek_bfvd.m8"

mkdir -p ${OUT_BFVD}
mkdir -p ${TMP_BFVD}

foldseek easy-search \
${QUERY} \
${DB_BFVD} \
${RESULT_BFVD} \
${TMP_BFVD} \
--threads 8 \
-e 0.01 \
--format-output query,target,alnlen,qstart,qend,tstart,tend,evalue,bits,tcov,qcov

echo "Foldseek searches finished"
echo "Results written to:"
echo "${RESULT_SWISSPROT}"
echo "${RESULT_UNIPROT50}"
echo "${RESULT_BFVD}"
