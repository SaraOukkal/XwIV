import sys
from collections import defaultdict
from Bio import SeqIO

orf_table = sys.argv[1]
fasta_file = sys.argv[2]
output_file = sys.argv[3]

# Lecture des longueurs de scaffolds avec restauration des "|"
scaffold_lengths = {}
for rec in SeqIO.parse(fasta_file, "fasta"):
    scaffold_id = rec.id.replace("-", "|")
    scaffold_lengths[scaffold_id] = len(rec.seq)

# Initialisation des masques par scaffold
scaffold_masks = {scaffold: [0] * length for scaffold, length in scaffold_lengths.items()}

# Lecture du TSV d'ORFs
with open(orf_table) as f:
    next(f)  # skip header
    for line in f:
        scaffold, orf_number, start, end = line.strip().split("\t")
        start, end = int(start), int(end)
        if scaffold in scaffold_masks:
            for i in range(min(start, end) - 1, max(start, end)):
                if 0 <= i < scaffold_lengths[scaffold]:
                    scaffold_masks[scaffold][i] = 1

# Écriture du fichier de sortie
with open(output_file, "w") as out:
    out.write("scaffold\tlength\tcoding_bp\tcoding_pct\n")
    for scaffold, mask in scaffold_masks.items():
        coding_bp = sum(mask)
        length = scaffold_lengths[scaffold]
        pct = coding_bp / length * 100 if length > 0 else 0
        out.write(f"{scaffold}\t{length}\t{coding_bp}\t{pct:.2f}\n")

