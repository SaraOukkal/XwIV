#!/usr/bin/env python3

import os
from pathlib import Path
from collections import defaultdict
from Bio import SeqIO

# ==== INPUT FILES ====

refseq_tsv = "/beegfs/data/soukkal/Thesis/Hymenoptera_Project/PDV/Results/Viral_genes/IVSPER_Analysis/xiphosomella_wirra_specimen2/Homology_search/ORFs_vs_RefseqVirus/ORFs_Refseq_Virus_homologs.tsv"
uniref_tsv = "/beegfs/data/soukkal/Thesis/Hymenoptera_Project/PDV/Results/Viral_genes/IVSPER_Analysis/xiphosomella_wirra_specimen2/Homology_search/ORFs_vs_UniRef50/ORFs_UniRef50_homologs.tsv"
iv_tsv = "/beegfs/data/soukkal/Thesis/Hymenoptera_Project/PDV/Results/Viral_genes/IVSPER_Analysis/xiphosomella_wirra_specimen2/Homology_search/ORFs_vs_IVgenes/ORFs_IVSPER_homologs.tsv"

orfs_fasta = "/beegfs/data/soukkal/Thesis/Hymenoptera_Project/PDV/Results/Viral_genes/IVSPER_Analysis/xiphosomella_wirra_specimen2/ORF_prediction/scaffold3_ORFs_70AA.faa"
refseq_fasta = "/beegfs/data/soukkal/Thesis/Databases/Refseq_Virus_formatted.faa"
uniref_fasta = "/beegfs/data/soukkal/Thesis/Databases/uniref50.fasta"
iv_fasta = "/beegfs/data/soukkal/Thesis/Hymenoptera_Project/PDV/Data/IV_domesticated.faa"

output_dir = Path("/beegfs/data/soukkal/Thesis/Hymenoptera_Project/PDV/Results/Viral_genes/IVSPER_Analysis/Phylogeny/Data")
output_dir.mkdir(parents=True, exist_ok=True)

# ==== FUNCTIONS ====

def parse_mapping_file(path):
    orf_map = defaultdict(list)
    with open(path) as f:
        for line in f:
            parts = line.strip().split()
            if len(parts) != 2:
                continue
            orfname = parts[0].split("_")[0]
            orf_map[orfname].append(parts[1])
    return orf_map

def get_required_ids(mapping):
    return set(val for v in mapping.values() for val in v)

def extract_uniref_sequences(fasta_path, accession_list):
    wanted = set(accession_list)
    results = defaultdict(str)
    current_header = None
    current_lines = []

    with open(fasta_path) as f:
        for line in f:
            if line.startswith(">"):
                if current_header and current_header in wanted:
                    results[current_header] = "".join(current_lines)
                header_line = line.strip()
                acc = header_line[1:].split()[0]
                current_header = acc
                current_lines = [line] if acc in wanted else []
            else:
                if current_header in wanted:
                    current_lines.append(line)
        if current_header and current_header in wanted:
            results[current_header] = "".join(current_lines)
    return results

def write_or_append_fasta(orf_name, seq_records, raw_strings=[]):
    output_path = output_dir / f"{orf_name}.faa"
    mode = "a" if output_path.exists() else "w"
    with open(output_path, mode) as f:
        if seq_records:
            SeqIO.write(seq_records, f, "fasta")
        if raw_strings:
            for s in raw_strings:
                f.write(s)

# ==== MAIN ====

print("Indexing ORF sequences...")
orf_index = {}
for record in SeqIO.parse(orfs_fasta, "fasta"):
    name = record.id.split("|")[-1].split("_")[0]
    orf_index[name] = record

print("Parsing mapping files...")
refseq_map = parse_mapping_file(refseq_tsv)
uniref_map = parse_mapping_file(uniref_tsv)
iv_map = parse_mapping_file(iv_tsv)

print("Preparing accession lists...")
refseq_ids = get_required_ids(refseq_map)
uniref_ids = get_required_ids(uniref_map)
iv_ids = get_required_ids(iv_map)

print("Indexing RefSeq...")
refseq_index = {}
for record in SeqIO.parse(refseq_fasta, "fasta"):
    acc = record.id.split("-")[0]  # ex: YP_010776174.1
    if acc in refseq_ids:
        refseq_index[acc] = record

print("Indexing IVSPER...")
iv_index = {}
for record in SeqIO.parse(iv_fasta, "fasta"):
    raw_id = record.id
    id_clean = raw_id.replace("lcl|", "")
    if raw_id in iv_ids or id_clean in iv_ids:
        iv_index[raw_id] = record
        iv_index[id_clean] = record

print("Extracting UniRef entries...")
uniref_entries = extract_uniref_sequences(uniref_fasta, uniref_ids)

# On ne conserve que les ORFs avec au moins un homologue
all_orfs = set(refseq_map) | set(uniref_map) | set(iv_map)

print("Processing ORFs with homologues only...")
for orf in sorted(all_orfs):
    seqs = []
    raw_fasta = []

    if orf in orf_index:
        seqs.append(orf_index[orf])
    else:
        print(f"[WARN] ORF sequence not found: {orf}")

    for acc in refseq_map.get(orf, []):
        if acc in refseq_index:
            seqs.append(refseq_index[acc])
        else:
            print(f"[WARN] RefSeq not found: {acc} (ORF {orf})")

    for acc in iv_map.get(orf, []):
        if acc in iv_index:
            seqs.append(iv_index[acc])
        else:
            print(f"[WARN] IVSPER not found: {acc} (ORF {orf})")

    for acc in uniref_map.get(orf, []):
        if acc in uniref_entries:
            raw_fasta.append(uniref_entries[acc])
        else:
            print(f"[WARN] UniRef not found: {acc} (ORF {orf})")

    if seqs or raw_fasta:
        write_or_append_fasta(orf, seqs, raw_strings=raw_fasta)
    else:
        print(f"[INFO] No sequences written for ORF {orf}")

print("Done.")

