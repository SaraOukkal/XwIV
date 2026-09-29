#!/bin/bash

#SBATCH --job-name=XwIV_DUST_TRF
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

SAMTOOLS="/beegfs/data/soft/samtools-1.19/bin/samtools"
TRF="/beegfs/home/soukkal/Tools/trf"

ASSEMBLY="/beegfs/data/soukkal/Thesis/Hymenoptera_Project/PDV/Results/Viral_genes/IVSPER_Analysis/Long_reads/07_Hybrid_assembly/MaSuRCA/CA.mr.99.17.15.0.02/primary.genome.scf.fasta"

CONTIG="jcf7180000126179"

OUTPUT_DIR="/beegfs/data/soukkal/Thesis/Hymenoptera_Project/PDV/Results/Viral_genes/IVSPER_Analysis/Long_reads/07_Hybrid_assembly/Terminal_complexity"

VIRAL_GENOME="${OUTPUT_DIR}/XwIV_hybrid_genome.fa"

mkdir -p "${OUTPUT_DIR}"


# ============================================================
# CHECK PROGRAMS
# ============================================================

if [[ ! -x "${TRF}" ]]; then
    echo "ERROR: TRF not found or not executable:"
    echo "${TRF}"
    exit 1
fi

if ! python -c "import pydustmasker" 2>/dev/null; then
    echo "ERROR: pydustmasker cannot be imported with the current Python."
    echo "Python used:"
    which python
    exit 1
fi


# ============================================================
# EXTRACT XwIV GENOME
# ============================================================

echo
echo "============================================"
echo "Extracting XwIV genome"
echo "============================================"

"${SAMTOOLS}" faidx "${ASSEMBLY}"

"${SAMTOOLS}" faidx \
    "${ASSEMBLY}" \
    "${CONTIG}" \
    > "${VIRAL_GENOME}"

"${SAMTOOLS}" faidx "${VIRAL_GENOME}"

GENOME_LENGTH=$(cut -f2 "${VIRAL_GENOME}.fai")

echo "Contig: ${CONTIG}"
echo "Genome length: ${GENOME_LENGTH} bp"


# ============================================================
# PYDUSTMASKER
# ============================================================

echo
echo "============================================"
echo "Running pydustmasker"
echo "============================================"

cat > "${OUTPUT_DIR}/run_pydustmasker.py" <<'PYTHON'

from pathlib import Path
from pydustmasker import DustMasker


input_fasta = Path("XwIV_hybrid_genome.fa")
output_file = Path("XwIV_DUST_regions.tsv")


# ------------------------------------------------------------
# Read FASTA
# ------------------------------------------------------------

sequence = []

with open(input_fasta) as f:
    for line in f:
        if line.startswith(">"):
            continue
        sequence.append(line.strip())

sequence = "".join(sequence).upper()


# ------------------------------------------------------------
# Run DUST
# ------------------------------------------------------------

masker = DustMasker(sequence)

regions = masker.intervals


# ------------------------------------------------------------
# Write coordinates
#
# pydustmasker intervals use Python-style coordinates.
# Convert to 1-based inclusive genomic coordinates.
# ------------------------------------------------------------

with open(output_file, "w") as out:

    out.write("start\tend\tlength\n")

    for start, end in regions:

        genomic_start = start + 1
        genomic_end = end

        length = genomic_end - genomic_start + 1

        out.write(
            f"{genomic_start}\t"
            f"{genomic_end}\t"
            f"{length}\n"
        )


# ------------------------------------------------------------
# Summary
# ------------------------------------------------------------

covered = set()

for start, end in regions:
    covered.update(range(start, end))

masked_bp = len(covered)

percent = (
    100 * masked_bp / len(sequence)
    if len(sequence) > 0
    else 0
)

print()
print("DUST SUMMARY")
print("------------")
print(f"Genome length: {len(sequence)} bp")
print(f"Low-complexity regions: {len(regions)}")
print(f"Low-complexity bp: {masked_bp}")
print(f"Low-complexity percentage: {percent:.3f}%")

PYTHON


cd "${OUTPUT_DIR}"

python run_pydustmasker.py


# ============================================================
# TANDEM REPEATS FINDER
# ============================================================

echo
echo "============================================"
echo "Running Tandem Repeats Finder"
echo "============================================"

# Standard TRF parameters:
#
# Match       = 2
# Mismatch    = 7
# Delta       = 7
# PM          = 80
# PI          = 10
# MinScore    = 50
# MaxPeriod   = 500
#
# -d = data output
# -h = suppress HTML output

"${TRF}" \
    "${VIRAL_GENOME}" \
    2 7 7 80 10 50 500 \
    -d \
    -h


# ============================================================
# FIND TRF OUTPUT
# ============================================================

FASTA_NAME=$(basename "${VIRAL_GENOME}")

TRF_RAW="${FASTA_NAME}.2.7.7.80.10.50.500.dat"

if [[ ! -s "${TRF_RAW}" ]]; then
    echo "ERROR: TRF output not found:"
    echo "${OUTPUT_DIR}/${TRF_RAW}"
    exit 1
fi


# ============================================================
# CONVERT TRF OUTPUT TO CLEAN TSV
# ============================================================

TRF_TSV="${OUTPUT_DIR}/XwIV_TRF_regions.tsv"

echo -e "start\tend\tlength\tperiod_size\tcopy_number\tconsensus_size\tpercent_matches\tpercent_indels\tscore\tA_percent\tC_percent\tG_percent\tT_percent\tentropy\tconsensus_sequence" > "${TRF_TSV}"

awk '
    /^[0-9]/ {

        start=$1
        end=$2
        length=end-start+1

        print start "\t" \
              end "\t" \
              length "\t" \
              $3 "\t" \
              $4 "\t" \
              $5 "\t" \
              $6 "\t" \
              $7 "\t" \
              $8 "\t" \
              $9 "\t" \
              $10 "\t" \
              $11 "\t" \
              $12 "\t" \
              $13 "\t" \
              $14
    }
' "${TRF_RAW}" >> "${TRF_TSV}"


# ============================================================
# TRF SUMMARY
# ============================================================

echo
echo "TRF SUMMARY"
echo "-----------"

awk -F'\t' '
    NR > 1 {
        n++
        total += $3

        if (n == 1 || $3 > max) {
            max = $3
            max_start = $1
            max_end = $2
        }
    }

    END {
        print "Number of tandem repeats:", n+0
        print "Sum of repeat lengths:", total+0, "bp"

        if (n > 0) {
            print "Longest tandem-repeat region:", max, "bp"
            print "Coordinates:", max_start "-" max_end
        }
    }
' "${TRF_TSV}"


# ============================================================
# FINAL OUTPUT
# ============================================================

echo
echo "============================================"
echo "ANALYSIS COMPLETE"
echo "============================================"

echo
echo "Genome:"
echo "${VIRAL_GENOME}"

echo
echo "DUST regions:"
echo "${OUTPUT_DIR}/XwIV_DUST_regions.tsv"

echo
echo "TRF regions:"
echo "${TRF_TSV}"

echo
echo "Raw TRF output:"
echo "${OUTPUT_DIR}/${TRF_RAW}"

echo
echo "Genome length: ${GENOME_LENGTH} bp"

echo
echo "============================================"
