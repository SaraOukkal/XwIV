import sys
from Bio import SeqIO

fasta_path = sys.argv[1]
output_tsv = sys.argv[2]

with open(output_tsv, "w") as out:
    out.write("scaffold\torf_number\tstart\tend\n")
    for record in SeqIO.parse(fasta_path, "fasta"):
        header = record.description
        try:
            # Exemple d'en-tête (outfmt 0) :
            # >lcl|ORF2_scaffold21070|size2618:1738:1974 unnamed protein product
            parts = header.split()
            main = parts[0]  # lcl|ORF2_scaffold21070|size2618:1738:1974

            main = main.replace("lcl|", "")
            orf_and_scaffold, coords = main.split(":")[0], main.split(":")[1:]
            start, end = map(int, coords)

            orf_tag, scaffold = orf_and_scaffold.split("_", 1)
            orf_number = orf_tag.replace("ORF", "")

            # Restauration du séparateur pipe dans le scaffold
            scaffold = scaffold.replace("-", "|")

            out.write(f"{scaffold}\t{orf_number}\t{start}\t{end}\n")
        except Exception as e:
            print(f"Erreur de parsing sur l'en-tête : {header} -> {e}")

