#!/bin/bash

#SBATCH --job-name=Xwirra_hybrid
#SBATCH --output=/beegfs/data/soukkal/Thesis/Hymenoptera_Project/PDV/Results/Viral_genes/IVSPER_Analysis/Long_reads/07_Hybrid_assembly/Logs/hybrid_%j.out
#SBATCH --error=/beegfs/data/soukkal/Thesis/Hymenoptera_Project/PDV/Results/Viral_genes/IVSPER_Analysis/Long_reads/07_Hybrid_assembly/Logs/hybrid_%j.err
#SBATCH --cpus-per-task=16
#SBATCH --mem=64G
#SBATCH --time=120:00:00
#SBATCH --constraint='skylake|haswell|broadwell'
#SBATCH --exclude=pbil-deb27

set -euo pipefail

export PATH="$HOME/local/numactl/usr/bin:$PATH"


########################
# General paths        #
########################

BASE="/beegfs/data/soukkal/Thesis/Hymenoptera_Project/PDV/Results/Viral_genes/IVSPER_Analysis"

LONG_BASE="${BASE}/Long_reads"

OUTDIR="${LONG_BASE}/07_Hybrid_assembly"

READDIR="${OUTDIR}/Reads"

ASSEMBLY_DIR="${OUTDIR}/MaSuRCA"

LOGDIR="${OUTDIR}/Logs"

THREADS=24


########################
# Software             #
########################

MASURCA="/beegfs/data/soft/MaSuRCA-4.0.8/bin/masurca"

SINGULARITY="/beegfs/data/soft/singularity3.11.0/bin/singularity"

BUSCO_SIF="/beegfs/home/soukkal/BUSCO/busco_v5.5.0_cv1.sif"

BUSCO_LINEAGE="/beegfs/data/soukkal/Thesis/Databases/busco_downloads/lineages/insecta_odb10"

BUSCO_DOWNLOAD="/beegfs/data/soukkal/Thesis/Databases/busco_downloads/"


########################
# Input reads          #
########################

ILLUMINA_DIR="${BASE}/Xiphosomella_wirra_assembly"

ILLUMINA_1_R1="${ILLUMINA_DIR}/DHJPAR0036296_1.fastq.gz"

ILLUMINA_1_R2="${ILLUMINA_DIR}/DHJPAR0036296_2.fastq.gz"

ILLUMINA_2_R1="${ILLUMINA_DIR}/DHJPAR0041296_1.fastq.gz"

ILLUMINA_2_R2="${ILLUMINA_DIR}/DHJPAR0041296_2.fastq.gz"

ONT="${LONG_BASE}/01_Reads/Xiphosomella_wirra_ONT.fastq.gz"


########################
# Combined reads       #
########################

R1="${READDIR}/Xwirra_combined_R1.fastq.gz"

R2="${READDIR}/Xwirra_combined_R2.fastq.gz"


########################
# Create directories   #
########################

mkdir -p "${READDIR}"
mkdir -p "${ASSEMBLY_DIR}"
mkdir -p "${LOGDIR}"


########################
# Check input files    #
########################

for FILE in "${ILLUMINA_1_R1}" "${ILLUMINA_1_R2}" "${ILLUMINA_2_R1}" "${ILLUMINA_2_R2}" "${ONT}"
do
	if [ ! -s "${FILE}" ]; then
		echo "ERROR: missing or empty input file: ${FILE}"
		exit 1
	fi
done


####################################
# Concatenate Illumina specimens   #
####################################

echo "Concatenating Illumina R1 files..."

cat "${ILLUMINA_1_R1}" "${ILLUMINA_2_R1}" > "${R1}"

echo "Concatenating Illumina R2 files..."

cat "${ILLUMINA_1_R2}" "${ILLUMINA_2_R2}" > "${R2}"


####################################
# MaSuRCA configuration            #
####################################

cd "${ASSEMBLY_DIR}"

CONFIG="${ASSEMBLY_DIR}/masurca_config.txt"

cat > "${CONFIG}" <<EOF
DATA
PE= pe 350 50 ${R1} ${R2}
NANOPORE=${ONT}
END

PARAMETERS
GRAPH_KMER_SIZE = auto
USE_LINKING_MATES = 1
LIMIT_JUMP_COVERAGE = 300
CA_PARAMETERS = cgwErrorRate=0.15
KMER_COUNT_THRESHOLD = 1
NUM_THREADS = ${THREADS}
JF_SIZE = 2000000000
SOAP_ASSEMBLY = 0
END
EOF


####################################
# Generate MaSuRCA assembly script #
####################################

echo "Generating MaSuRCA assembly pipeline..."

"${MASURCA}" "${CONFIG}"


if [ ! -f "${ASSEMBLY_DIR}/assemble.sh" ]; then
	echo "ERROR: MaSuRCA did not generate assemble.sh"
	exit 1
fi


####################################
# Run hybrid assembly              #
####################################

echo "Starting MaSuRCA hybrid assembly..."

bash "${ASSEMBLY_DIR}/assemble.sh"


####################################
# Locate final assembly            #
####################################

FINAL_ASSEMBLY="${ASSEMBLY_DIR}/CA/final.genome.scf.fasta"

if [ ! -s "${FINAL_ASSEMBLY}" ]; then
	echo "ERROR: final.genome.scf.fasta was not found."
	exit 1
fi

cp "${FINAL_ASSEMBLY}" "${OUTDIR}/Xiphosomella_wirra_hybrid_MaSuRCA.fa"

HYBRID_ASSEMBLY="${OUTDIR}/Xiphosomella_wirra_hybrid_MaSuRCA.fa"


####################################
# BUSCO                            #
####################################

BUSCO_PARENT="${OUTDIR}/BUSCO"

BUSCO_NAME="BUSCO_Xiphosomella_wirra_hybrid"

mkdir -p "${BUSCO_PARENT}"

cd "${BUSCO_PARENT}"

rm -rf "${BUSCO_NAME}"

echo "Running BUSCO..."

"${SINGULARITY}" exec --bind /beegfs/:/beegfs/ "${BUSCO_SIF}" busco -i "${HYBRID_ASSEMBLY}" -l "${BUSCO_LINEAGE}" -m genome -o "${BUSCO_NAME}" -c "${THREADS}" -f --offline --download_path "${BUSCO_DOWNLOAD}"


####################################
# Done                             #
####################################

echo ""
echo "Hybrid assembly completed."
echo ""
echo "Final assembly:"
echo "${HYBRID_ASSEMBLY}"
echo ""
echo "BUSCO:"
echo "${BUSCO_PARENT}/${BUSCO_NAME}/run_insecta_odb10/full_table.tsv"
```

