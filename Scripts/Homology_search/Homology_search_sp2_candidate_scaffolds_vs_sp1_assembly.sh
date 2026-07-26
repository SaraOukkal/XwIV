#!/bin/bash
#SBATCH --job-name=sp2_candidates_vs_sp1
#SBATCH --output=%x_%j.out
#SBATCH --error=%x_%j.err
#SBATCH --time=01:00:00
#SBATCH --cpus-per-task=2
#SBATCH --mem=10G
#SBATCH --exclude=pbil-deb[14-27]

MMSEQS="/beegfs/data/soukkal/TOOLS/mmseqs/bin/mmseqs"
QUERY="/beegfs/data/soukkal/Thesis/Hymenoptera_Project/PDV/Results/Viral_genes/IVSPER_Analysis/xiphosomella_wirra_specimen2/sp2_candidate_scaffolds.fa"
TARGET="/beegfs/project/horizon/data/assembly/specimens/DHJPAR0041296/redundans/scaffolds.reduced.fa"
OUTPUT_DIR="/beegfs/data/soukkal/Thesis/Hymenoptera_Project/PDV/Results/Viral_genes/IVSPER_Analysis/xiphosomella_wirra_specimen2/Homology_search/sp2_candidate_scaffolds_vs_sp1_assembly"

# Création du dossier de sortie
mkdir -p $OUTPUT_DIR

# Création des bases de données
$MMSEQS createdb $QUERY $OUTPUT_DIR/candidates_DB
$MMSEQS createdb $TARGET $OUTPUT_DIR/assembly_DB

# Recherche (nucléotide vs nucléotide)
$MMSEQS search $OUTPUT_DIR/candidates_DB $OUTPUT_DIR/assembly_DB $OUTPUT_DIR/result $OUTPUT_DIR/tmp \
  --search-type 2 -e 0.01 -s 7.5 -a --threads 2 --remove-tmp-files

# Conversion du résultat au format m8
$MMSEQS convertalis $OUTPUT_DIR/candidates_DB $OUTPUT_DIR/assembly_DB $OUTPUT_DIR/result $OUTPUT_DIR/result.m8 \
  --format-output 'query,qlen,tlen,target,pident,alnlen,mismatch,gapopen,qstart,qend,tstart,tend,evalue,bits,qaln,tcov' \
  --search-type 2

