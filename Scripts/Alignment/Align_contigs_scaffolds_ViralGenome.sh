#!/bin/bash
#SBATCH --job-name=MMseqs_orient_and_align
#SBATCH --output=%x_%j.out
#SBATCH --error=%x_%j.err
#SBATCH --time=02:00:00
#SBATCH --cpus-per-task=4
#SBATCH --mem=10G
#SBATCH --exclude=pbil-deb[14-27]

set -euo pipefail
: "${SLURM_CPUS_PER_TASK:=1}"

# ================================
# Tools
# ================================
MMSEQS="/beegfs/data/soukkal/TOOLS/mmseqs/bin/mmseqs"
ALIGNER="/beegfs/home/soukkal/miniconda3/bin/clustalo"
ALIGNER_OPTS="--threads ${SLURM_CPUS_PER_TASK} --outfmt=fasta"
PYTHON="python3"

command -v "${MMSEQS}" >/dev/null 2>&1 || { echo "ERROR: mmseqs not found"; exit 1; }
command -v "${ALIGNER}" >/dev/null 2>&1 || { echo "ERROR: clustalo not found"; exit 1; }
command -v "${PYTHON}" >/dev/null 2>&1 || { echo "ERROR: python3 not found"; exit 1; }

# ================================
# Inputs
# ================================
VIRAL_GENOME="/beegfs/data/soukkal/Thesis/Hymenoptera_Project/PDV/Results/Viral_genes/IVSPER_Analysis/Viral_Genome/Viral_genome.fa"
SP1_CONTIGS="/beegfs/data/soukkal/Thesis/Hymenoptera_Project/PDV/Results/Viral_genes/IVSPER_Analysis/sp1_contigs.fasta"
SP2_CONTIGS="/beegfs/data/soukkal/Thesis/Hymenoptera_Project/PDV/Results/Viral_genes/IVSPER_Analysis/sp2_contigs.fasta"
SP1_SCAFFOLDS="/beegfs/data/soukkal/Thesis/Hymenoptera_Project/PDV/Results/Viral_genes/IVSPER_Analysis/Viral_candidate_scaffolds.fa"
SP2_SCAFFOLDS="/beegfs/data/soukkal/Thesis/Hymenoptera_Project/PDV/Results/Viral_genes/IVSPER_Analysis/xiphosomella_wirra_specimen2/scaffold3.fa"

OUTDIR="/beegfs/data/soukkal/Thesis/Hymenoptera_Project/PDV/Results/Viral_genes/IVSPER_Analysis/Viral_Genome/Alignements"
mkdir -p "${OUTDIR}"

for f in "$VIRAL_GENOME" "$SP1_CONTIGS" "$SP2_CONTIGS" "$SP1_SCAFFOLDS" "$SP2_SCAFFOLDS"; do
  [[ -s "$f" ]] || { echo "ERROR: Missing or empty $f"; exit 1; }
done

# ================================
# Step 1: Prefix contigs
# ================================
REF_FA="${OUTDIR}/Viral_genome.fa"
cp "${VIRAL_GENOME}" "${REF_FA}"

SP1_CONTIGS_PREF="${OUTDIR}/sp1_contigs.prefixed.fa"
SP2_CONTIGS_PREF="${OUTDIR}/sp2_contigs.prefixed.fa"
ALL_CONTIGS_PREF="${OUTDIR}/all_contigs.prefixed.fa"

${PYTHON} - "$SP1_CONTIGS" "$SP1_CONTIGS_PREF" << 'PY'
import sys
inp, outp = sys.argv[1], sys.argv[2]
with open(inp) as fin, open(outp, "w") as fout:
    for line in fin:
        if line.startswith(">"):
            fout.write(">sp1_contig|" + line.strip()[1:] + "\n")
        else:
            fout.write(line)
PY

${PYTHON} - "$SP2_CONTIGS" "$SP2_CONTIGS_PREF" << 'PY'
import sys
inp, outp = sys.argv[1], sys.argv[2]
with open(inp) as fin, open(outp, "w") as fout:
    for line in fin:
        if line.startswith(">"):
            fout.write(">sp2_contig|" + line.strip()[1:] + "\n")
        else:
            fout.write(line)
PY

cat "${SP1_CONTIGS_PREF}" "${SP2_CONTIGS_PREF}" > "${ALL_CONTIGS_PREF}"

# ================================
# Step 2: MMseqs DB creation
# ================================
DB_REF="${OUTDIR}/viral_db"
DB_CONTIGS="${OUTDIR}/contigs_db"
TMPDIR="${OUTDIR}/tmp"
mkdir -p "${TMPDIR}"

"${MMSEQS}" createdb "${REF_FA}" "${DB_REF}"
"${MMSEQS}" createdb "${ALL_CONTIGS_PREF}" "${DB_CONTIGS}"

# ================================
# Step 3: MMseqs search
# ================================
RES_BASE="${OUTDIR}/viral_vs_contigs"
"${MMSEQS}" search "${DB_REF}" "${DB_CONTIGS}" "${RES_BASE}" "${TMPDIR}" -a -s 7.5 -e 0.001 --threads ${SLURM_CPUS_PER_TASK} --remove-tmp-files --search-type 3

# ================================
# Step 4: convertalis to m8
# ================================
M8="${OUTDIR}/viral_vs_contigs.m8"
"${MMSEQS}" convertalis "${DB_REF}" "${DB_CONTIGS}" "${RES_BASE}" "${M8}" \
  --format-output 'query,qlen,tlen,target,pident,alnlen,mismatch,gapopen,qstart,qend,tstart,tend,evalue,bits,tcov' --search-type 3

# ================================
# Step 5: Orientation by qstart/qend
# ================================
ORIENTED_CONTIGS="${OUTDIR}/all_contigs.oriented.fa"
REVERSED_LOG="${OUTDIR}/contigs_reversed.tsv"

${PYTHON} - "$M8" "$ALL_CONTIGS_PREF" "$ORIENTED_CONTIGS" "$REVERSED_LOG" << 'PY'
import sys
m8_path, fasta_in, fasta_out, log_out = sys.argv[1:5]

best = {}
with open(m8_path) as f:
    for line in f:
        if not line.strip():
            continue
        cols = line.rstrip("\n").split("\t")
        if len(cols) < 15:
            continue
        target = cols[3]
        try:
            qstart = int(cols[8]); qend = int(cols[9]); bits = float(cols[13])
        except ValueError:
            continue
        if (target not in best) or (bits > best[target][0]):
            best[target] = (bits, qstart, qend)

def read_fasta(path):
    seqs = {}
    name = None
    buf = []
    with open(path) as f:
        for line in f:
            if line.startswith(">"):
                if name is not None:
                    seqs[name] = "".join(buf).upper()
                name = line.strip()[1:]
                buf = []
            else:
                buf.append(line.strip())
    if name is not None:
        seqs[name] = "".join(buf).upper()
    return seqs

COMP = str.maketrans("ACGTRYMKBDHVNacgtrymkbdhvn", "TGCAYRKMVHDBNtgcayrkmvhdbn")
def revcomp(s):
    return s.translate(COMP)[::-1]

seqs = read_fasta(fasta_in)
reversed_list = []

with open(fasta_out, "w") as out:
    for name, seq in seqs.items():
        strand = "+"
        if name in best:
            _, qs, qe = best[name]
            strand = "-" if qs > qe else "+"
        if strand == "-":
            new_name = name + "_RC"
            seq = revcomp(seq)
            reversed_list.append((name, new_name))
            out.write(">" + new_name + "\n")
        else:
            out.write(">" + name + "\n")
        for i in range(0, len(seq), 80):
            out.write(seq[i:i+80] + "\n")

with open(log_out, "w") as log:
    log.write("original_id\toriented_id\tstrand\n")
    for o, n in reversed_list:
        log.write(f"{o}\t{n}\t-\n")
PY

# ================================
# Step 6: Scaffold RC for scaffold33401|size845
# ================================
SP1_SCAFFOLDS_ORIENT="${OUTDIR}/sp1_scaffolds.oriented.fa"
${PYTHON} - "$SP1_SCAFFOLDS" "$SP1_SCAFFOLDS_ORIENT" << 'PY'
import sys
inp, outp = sys.argv[1], sys.argv[2]
target = "scaffold33401|size845"

COMP = str.maketrans("ACGTRYMKBDHVNacgtrymkbdhvn", "TGCAYRKMVHDBNtgcayrkmvhdbn")
def revcomp(s):
    return s.translate(COMP)[::-1]

def read_fasta(path):
    seqs = {}
    name = None
    buf = []
    with open(path) as f:
        for line in f:
            if line.startswith(">"):
                if name is not None:
                    seqs[name] = "".join(buf).upper()
                name = line.strip()[1:]
                buf = []
            else:
                buf.append(line.strip())
    if name is not None:
        seqs[name] = "".join(buf).upper()
    return seqs

seqs = read_fasta(inp)
with open(outp, "w") as out:
    for name, seq in seqs.items():
        if name == target:
            rc = revcomp(seq)
            out.write(">" + name + "_RC\n")
            for i in range(0, len(rc), 80):
                out.write(rc[i:i+80] + "\n")
        else:
            out.write(">" + name + "\n")
            for i in range(0, len(seq), 80):
                out.write(seq[i:i+80] + "\n")
PY

SP2_SCAFFOLDS_ORIENT="${OUTDIR}/sp2_scaffolds.oriented.fa"
cp "${SP2_SCAFFOLDS}" "${SP2_SCAFFOLDS_ORIENT}"

# ================================
# Step 7: Final MSA
# ================================
MSA_INPUT="${OUTDIR}/msa_input_on_viral.fa"
cat "${REF_FA}" "${SP1_SCAFFOLDS_ORIENT}" "${SP2_SCAFFOLDS_ORIENT}" "${ORIENTED_CONTIGS}" > "${MSA_INPUT}"

ALN="${OUTDIR}/alignment_on_viral_genome.aln"
"${ALIGNER}" -i "${MSA_INPUT}" -o "${ALN}" ${ALIGNER_OPTS}

echo "Done."
echo "  Oriented contigs FASTA: ${ORIENTED_CONTIGS}"
echo "  Reversed contigs log:   ${REVERSED_LOG}"
echo "  MSA written to:         ${ALN}"

