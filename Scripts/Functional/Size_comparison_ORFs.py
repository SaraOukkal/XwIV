#!/usr/bin/env python3

import csv
from pathlib import Path

# Input root directory containing the homology search subdirectories
INPUT_ROOT = Path("/beegfs/data/soukkal/Thesis/Hymenoptera_Project/PDV/Results/Viral_genes/IVSPER_Analysis/xiphosomella_wirra_specimen2/Homology_search/")

# Output directory and file
OUTPUT_DIR = Path("/beegfs/data/soukkal/Thesis/Hymenoptera_Project/PDV/Results/Viral_genes/IVSPER_Analysis/xiphosomella_wirra_specimen2/Size_comparison/")
OUTPUT_FILE = OUTPUT_DIR / "Predicted_proteins_size_comparison.tsv"

# Expected subdirectories
SUBDIRS = [
    "ORFs_vs_IVgenes",
    "ORFs_vs_RefseqVirus",
    "ORFs_vs_UniRef50"
]

# M8 columns
M8_COLUMNS = [
    "query", "qlen", "tlen", "target", "pident", "alnlen", "mismatch", "gapopen",
    "qstart", "qend", "tstart", "tend", "evalue", "bits", "qaln", "tcov"
]


def safe_float(value):
    try:
        return float(value)
    except Exception:
        return None


def safe_int(value):
    try:
        return int(float(value))
    except Exception:
        return None


def parse_db_name(subdir_name):
    if subdir_name.startswith("ORFs_vs_"):
        return subdir_name.replace("ORFs_vs_", "", 1)
    return subdir_name


def is_better_hit(current, best):
    """
    Return True if current hit is better than best hit.
    Ranking:
    1. Higher tcov
    2. Lower evalue
    """
    if best is None:
        return True

    if current["tcov"] != best["tcov"]:
        return current["tcov"] > best["tcov"]

    return current["evalue"] < best["evalue"]


def process_m8_file(m8_file, db_name):
    best_hits = {}

    with open(m8_file, "r") as f:
        reader = csv.reader(f, delimiter="\t")
        for line_number, row in enumerate(reader, start=1):
            if not row:
                continue

            if len(row) != len(M8_COLUMNS):
                print(f"Warning: skipping malformed line {line_number} in {m8_file}")
                continue

            record = dict(zip(M8_COLUMNS, row))

            query = record["query"]
            target = record["target"]

            qlen = safe_int(record["qlen"])
            tlen = safe_int(record["tlen"])
            alnlen = safe_int(record["alnlen"])
            tcov = safe_float(record["tcov"])
            evalue = safe_float(record["evalue"])

            if None in (qlen, tlen, alnlen, tcov, evalue):
                print(f"Warning: skipping line {line_number} with invalid numeric values in {m8_file}")
                continue

            length_ratio = qlen / tlen if tlen > 0 else None
            qcov = alnlen / qlen if qlen > 0 else None

            hit = {
                "predicted_protein": query,
                "db_protein": target,
                "db": db_name,
                "predicted_length": qlen,
                "db_length": tlen,
                "length_ratio": length_ratio,
                "qcov": qcov,
                "tcov": tcov,
                "evalue": evalue
            }

            if query not in best_hits or is_better_hit(hit, best_hits[query]):
                best_hits[query] = hit

    return list(best_hits.values())


def main():
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

    all_results = []

    for subdir in SUBDIRS:
        m8_file = INPUT_ROOT / subdir / "result.m8"
        db_name = parse_db_name(subdir)

        if not m8_file.exists():
            print(f"Warning: file not found, skipping: {m8_file}")
            continue

        print(f"Processing {m8_file}")
        results = process_m8_file(m8_file, db_name)
        all_results.extend(results)

    all_results.sort(key=lambda x: (x["predicted_protein"], x["db"]))

    output_columns = [
        "predicted_protein",
        "db_protein",
        "db",
        "predicted_length",
        "db_length",
        "length_ratio",
        "qcov",
        "tcov"
    ]

    with open(OUTPUT_FILE, "w", newline="") as out:
        writer = csv.DictWriter(out, fieldnames=output_columns, delimiter="\t")
        writer.writeheader()
        for row in all_results:
            writer.writerow({
                "predicted_protein": row["predicted_protein"],
                "db_protein": row["db_protein"],
                "db": row["db"],
                "predicted_length": row["predicted_length"],
                "db_length": row["db_length"],
                "length_ratio": f"{row['length_ratio']:.6f}" if row["length_ratio"] is not None else "NA",
                "qcov": f"{row['qcov']:.6f}" if row["qcov"] is not None else "NA",
                "tcov": f"{row['tcov']:.6f}" if row["tcov"] is not None else "NA"
            })

    print(f"Done. Output written to: {OUTPUT_FILE}")


if __name__ == "__main__":
    main()
