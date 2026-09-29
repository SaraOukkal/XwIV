#!/bin/bash

#SBATCH --job-name=XwIV_hybrid_ORFs
#SBATCH --output=/beegfs/data/soukkal/Thesis/Hymenoptera_Project/PDV/Results/Viral_genes/IVSPER_Analysis/Long_reads/07_Hybrid_assembly/ORF_prediction/Logs/orfs_%j.out
#SBATCH --error=/beegfs/data/soukkal/Thesis/Hymenoptera_Project/PDV/Results/Viral_genes/IVSPER_Analysis/Long_reads/07_Hybrid_assembly/ORF_prediction/Logs/orfs_%j.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --time=04:00:00
#SBATCH --mem=16G

set -euo pipefail


#################
# Paths         #
#################

BASE="/beegfs/data/soukkal/Thesis/Hymenoptera_Project/PDV/Results/Viral_genes/IVSPER_Analysis"

HYBRID_ASSEMBLY="${BASE}/Long_reads/07_Hybrid_assembly/MaSuRCA/CA.mr.99.17.15.0.02/primary.genome.scf.fasta"

OLD_ORFS="${BASE}/xiphosomella_wirra_specimen2/ORF_prediction/scaffold3_ORFs_70AA_filtered.faa"

OUTDIR="${BASE}/Long_reads/07_Hybrid_assembly/ORF_prediction"

LOGDIR="${OUTDIR}/Logs"

MMSEQS_DIR="${OUTDIR}/MMseqs"

TMPDIR="${MMSEQS_DIR}/tmp"

VIRAL_SCAFFOLD="jcf7180000126179"


#################
# Software      #
#################

ORFFINDER="ORFfinder"

MMSEQS="/beegfs/data/soukkal/TOOLS/mmseqs/bin/mmseqs"

SAMTOOLS="/beegfs/data/soft/samtools-1.19/bin/samtools"

THREADS=8


#################
# Output files  #
#################

VIRAL_FASTA="${OUTDIR}/XwIV_hybrid_scaffold.fa"

RAW_ORFS="${OUTDIR}/XwIV_hybrid_ORFfinder_70AA_raw.faa"

RAW_ORF_TABLE="${OUTDIR}/XwIV_hybrid_ORFfinder_70AA_raw.tsv"

MMSEQS_ALL="${MMSEQS_DIR}/old_ORFs_vs_hybrid_all_hits.tsv"

MMSEQS_BEST="${OUTDIR}/old_ORFs_vs_hybrid_best_hits.tsv"

NEW_ORFS="${OUTDIR}/XwIV_new_ORFs.faa"

NEW_ORF_TABLE="${OUTDIR}/XwIV_new_ORFs.tsv"

FINAL_ORFS="${OUTDIR}/XwIV_hybrid_final_ORFs.faa"

FINAL_TABLE="${OUTDIR}/XwIV_hybrid_final_ORFs.tsv"


#################
# Directories   #
#################

mkdir -p "${OUTDIR}"
mkdir -p "${LOGDIR}"
mkdir -p "${MMSEQS_DIR}"
mkdir -p "${TMPDIR}"


#########################################
# 1. Extract hybrid XwIV scaffold       #
#########################################

echo "Extracting hybrid XwIV scaffold..."

${SAMTOOLS} faidx "${HYBRID_ASSEMBLY}"

${SAMTOOLS} faidx \
    "${HYBRID_ASSEMBLY}" \
    "${VIRAL_SCAFFOLD}" \
    > "${VIRAL_FASTA}"


#################################################
# 2. MMseqs database: previous XwIV ORFs        #
#################################################

echo "Creating MMseqs database for previous XwIV ORFs..."

${MMSEQS} createdb \
    "${OLD_ORFS}" \
    "${MMSEQS_DIR}/old_ORFs_DB"


#################################################
# 3. MMseqs translated database: hybrid XwIV    #
#################################################

echo "Creating nucleotide database for hybrid XwIV..."

${MMSEQS} createdb \
    "${VIRAL_FASTA}" \
    "${MMSEQS_DIR}/hybrid_XwIV_DB"


#################################################
# 4. Search old proteins against hybrid genome  #
#################################################

echo "Searching previous XwIV ORFs against hybrid XwIV genome..."

${MMSEQS} search \
    "${MMSEQS_DIR}/old_ORFs_DB" \
    "${MMSEQS_DIR}/hybrid_XwIV_DB" \
    "${MMSEQS_DIR}/old_vs_hybrid_result" \
    "${TMPDIR}" \
    --search-type 2 \
    --threads ${THREADS}


#################################################
# 5. Export MMseqs hits                         #
#################################################

${MMSEQS} convertalis \
    "${MMSEQS_DIR}/old_ORFs_DB" \
    "${MMSEQS_DIR}/hybrid_XwIV_DB" \
    "${MMSEQS_DIR}/old_vs_hybrid_result" \
    "${MMSEQS_ALL}" \
    --format-output "query,target,pident,alnlen,qstart,qend,qlen,tstart,tend,tlen,evalue,bits"


#################################################
# 6. Select best hit for every previous ORF     #
#################################################

echo "Selecting best hit for each previous ORF..."

python3 - <<EOF

import pandas as pd

input_file = "${MMSEQS_ALL}"
output_file = "${MMSEQS_BEST}"

columns = [
    "Previous_ORF",
    "target",
    "pident",
    "alignment_length",
    "qstart",
    "qend",
    "qlen",
    "tstart",
    "tend",
    "tlen",
    "evalue",
    "bitscore"
]

df = pd.read_csv(
    input_file,
    sep="\t",
    names=columns
)

# Query coverage
df["query_coverage"] = (
    abs(df["qend"] - df["qstart"]) + 1
) / df["qlen"]

# Determine genomic strand from target coordinates.
df["strand"] = df.apply(
    lambda x: "+" if x["tend"] >= x["tstart"] else "-",
    axis=1
)

# Normalized genomic coordinates.
df["start"] = df[["tstart", "tend"]].min(axis=1)
df["end"] = df[["tstart", "tend"]].max(axis=1)

# Highest bitscore = best hit.
# E-value is used as secondary criterion.
df = df.sort_values(
    ["Previous_ORF", "bitscore", "evalue"],
    ascending=[True, False, True]
)

best = df.drop_duplicates(
    subset="Previous_ORF",
    keep="first"
).copy()

best = best[
    [
        "Previous_ORF",
        "start",
        "end",
        "strand",
        "pident",
        "alignment_length",
        "query_coverage",
        "evalue",
        "bitscore"
    ]
]

best.to_csv(
    output_file,
    sep="\t",
    index=False
)

print(f"Previous ORFs with a MMseqs hit: {len(best)}")

EOF


#########################################
# 7. ORFfinder                         #
#########################################

echo "Running ORFfinder..."

# Minimum length:
# 70 amino acids = 210 nucleotides

${ORFFINDER} \
    -in "${VIRAL_FASTA}" \
    -out "${RAW_ORFS}" \
    -ml 210 \
    -outfmt 0


#################################################
# 8. Parse ORFfinder + identify NEW ORFs        #
#################################################

echo "Identifying additional ORFs..."

python3 - <<EOF

import pandas as pd
from Bio import SeqIO
import re

orf_fasta = "${RAW_ORFS}"
best_hits_file = "${MMSEQS_BEST}"

raw_table = "${RAW_ORF_TABLE}"

new_fasta = "${NEW_ORFS}"
new_table = "${NEW_ORF_TABLE}"

final_fasta = "${FINAL_ORFS}"
final_table = "${FINAL_TABLE}"

scaffold = "${VIRAL_SCAFFOLD}"


########################################
# Load previous ORFs located by MMseqs #
########################################

old = pd.read_csv(
    best_hits_file,
    sep="\t"
)

old_regions = []

for _, row in old.iterrows():

    old_regions.append({
        "previous_orf": row["Previous_ORF"],
        "start": int(row["start"]),
        "end": int(row["end"]),
        "strand": row["strand"]
    })


########################################
# Parse ORFfinder predictions          #
########################################

predicted = []

for record in SeqIO.parse(orf_fasta, "fasta"):

    header = record.description

    match = re.search(
        r'ORF(\d+)_.*?:(\d+):(\d+)',
        header
    )

    if not match:
        raise ValueError(
            f"Cannot parse ORFfinder header: {header}"
        )

    original_number = int(match.group(1))

    coord1 = int(match.group(2))
    coord2 = int(match.group(3))

    strand = "+" if coord2 >= coord1 else "-"

    start = min(coord1, coord2)
    end = max(coord1, coord2)

    length_nt = end - start + 1

    predicted.append({
        "original_number": original_number,
        "start": start,
        "end": end,
        "strand": strand,
        "length_nt": length_nt,
        "record": record
    })


########################################
# Write complete ORFfinder table       #
########################################

with open(raw_table, "w") as out:

    out.write(
        "ORFfinder_ORF\tstart\tend\tstrand\tlength_nt\n"
    )

    for x in predicted:

        out.write(
            f"ORF{x['original_number']}\t"
            f"{x['start']}\t"
            f"{x['end']}\t"
            f"{x['strand']}\t"
            f"{x['length_nt']}\n"
        )


########################################
# Overlap function                     #
########################################

def overlap_bp(a_start, a_end, b_start, b_end):

    return max(
        0,
        min(a_end, b_end)
        - max(a_start, b_start)
        + 1
    )


##################################################
# Remove ORFfinder ORFs already represented     #
# by previous XwIV ORFs                         #
##################################################

candidate_new = []

for candidate in predicted:

    represented = False

    for old_orf in old_regions:

        overlap = overlap_bp(
            candidate["start"],
            candidate["end"],
            old_orf["start"],
            old_orf["end"]
        )

        fraction = overlap / candidate["length_nt"]

        if fraction > 0.50:
            represented = True
            break

    if not represented:
        candidate_new.append(candidate)


##################################################
# Filter overlap among NEW ORFs                  #
##################################################

# Process longest ORFs first.
# If >50% of the shorter ORF overlaps a longer
# retained ORF, the shorter ORF is removed.

candidate_new = sorted(
    candidate_new,
    key=lambda x: (-x["length_nt"], x["start"])
)

kept_new = []

for candidate in candidate_new:

    conflicting = False

    for existing in kept_new:

        overlap = overlap_bp(
            candidate["start"],
            candidate["end"],
            existing["start"],
            existing["end"]
        )

        shorter = min(
            candidate["length_nt"],
            existing["length_nt"]
        )

        if overlap / shorter > 0.50:
            conflicting = True
            break

    if not conflicting:
        kept_new.append(candidate)


########################################
# Sort new ORFs by genomic position    #
########################################

kept_new = sorted(
    kept_new,
    key=lambda x: x["start"]
)


########################################
# Write NEW ORFs only                  #
########################################

with open(new_table, "w") as out:

    out.write(
        "New_candidate\tstart\tend\tstrand\tlength_nt\n"
    )

    for i, x in enumerate(kept_new, start=1):

        out.write(
            f"New_candidate_{i}\t"
            f"{x['start']}\t"
            f"{x['end']}\t"
            f"{x['strand']}\t"
            f"{x['length_nt']}\n"
        )


new_records = []

for i, x in enumerate(kept_new, start=1):

    record = x["record"]

    record.id = f"New_candidate_{i}"
    record.name = record.id

    record.description = (
        f"{record.id} "
        f"start={x['start']} "
        f"end={x['end']} "
        f"strand={x['strand']}"
    )

    new_records.append(record)

SeqIO.write(
    new_records,
    new_fasta,
    "fasta"
)


##################################################
# Combine previous + newly detected ORFs         #
##################################################

combined = []

for old_orf in old_regions:

    combined.append({
        "previous_orf": old_orf["previous_orf"],
        "start": old_orf["start"],
        "end": old_orf["end"],
        "strand": old_orf["strand"],
        "record": None,
        "source": "previous"
    })


for x in kept_new:

    combined.append({
        "previous_orf": "NA",
        "start": x["start"],
        "end": x["end"],
        "strand": x["strand"],
        "record": x["record"],
        "source": "new"
    })


########################################
# Sort all ORFs by position            #
########################################

combined = sorted(
    combined,
    key=lambda x: x["start"]
)


########################################
# Final numbering                      #
########################################

for i, x in enumerate(combined, start=1):

    x["new_name"] = f"XwIV_ORF_{i}"


########################################
# Final table                          #
########################################

with open(final_table, "w") as out:

    out.write(
        "New_ORF\tPrevious_ORF\tstart\tend\tstrand\n"
    )

    for x in combined:

        out.write(
            f"{x['new_name']}\t"
            f"{x['previous_orf']}\t"
            f"{x['start']}\t"
            f"{x['end']}\t"
            f"{x['strand']}\n"
        )


########################################
# Final protein FASTA                  #
########################################

# For previous ORFs, retrieve the original
# protein sequence from the old FASTA.

old_sequences = {
    record.id: record
    for record in SeqIO.parse("${OLD_ORFS}", "fasta")
}

final_records = []

for x in combined:

    if x["source"] == "previous":

        record = old_sequences[x["previous_orf"]]

    else:

        record = x["record"]

    record.id = x["new_name"]
    record.name = x["new_name"]

    record.description = (
        f"{x['new_name']} "
        f"previous_ORF={x['previous_orf']} "
        f"start={x['start']} "
        f"end={x['end']} "
        f"strand={x['strand']}"
    )

    final_records.append(record)


SeqIO.write(
    final_records,
    final_fasta,
    "fasta"
)


########################################
# Summary                              #
########################################

print("")
print("XwIV ORF reconstruction")
print("-----------------------")

print(
    f"Previous ORFs recovered by MMseqs: "
    f"{len(old_regions)}"
)

print(
    f"Additional ORFs predicted by ORFfinder: "
    f"{len(kept_new)}"
)

print(
    f"Final number of XwIV ORFs: "
    f"{len(combined)}"
)

EOF


###############################
# Done                        #
###############################

echo ""
echo "Analysis completed."
echo ""
echo "Best hits of previous XwIV ORFs:"
echo "${MMSEQS_BEST}"
echo ""
echo "Additional ORFs:"
echo "${NEW_ORF_TABLE}"
echo ""
echo "Final ORF table:"
echo "${FINAL_TABLE}"
echo ""
echo "Final protein FASTA:"
echo "${FINAL_ORFS}"
