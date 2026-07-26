#!/bin/bash
#SBATCH --job-name=scaffold3_vs_PDV_motifs
#SBATCH --output=%x_%j.out
#SBATCH --error=%x_%j.err
#SBATCH --time=01:00:00
#SBATCH --cpus-per-task=2
#SBATCH --mem=10G
#SBATCH --exclude=pbil-deb[14-27]

MMSEQS="/beegfs/data/soukkal/TOOLS/mmseqs/bin/mmseqs"
SCAFFOLD="/beegfs/data/soukkal/Thesis/Hymenoptera_Project/PDV/Results/Viral_genes/IVSPER_Analysis/xiphosomella_wirra_specimen2/scaffold3.fa"
MOTIFS="/beegfs/data/soukkal/Thesis/Hymenoptera_Project/PDV/Data/All_motifs.fa"
OUTPUT_DIR="/beegfs/data/soukkal/Thesis/Hymenoptera_Project/PDV/Results/Viral_genes/IVSPER_Analysis/xiphosomella_wirra_specimen2/Homology_search/scaffold3_vs_PDV_motifs"

# Création du dossier de sortie
mkdir -p $OUTPUT_DIR

# Création des bases de données
$MMSEQS createdb $SCAFFOLD $OUTPUT_DIR/scaffold3_DB
$MMSEQS createdb $MOTIFS $OUTPUT_DIR/motifs_DB

# Recherche (nucléotide vs nucléotide)
$MMSEQS search $OUTPUT_DIR/scaffold3_DB $OUTPUT_DIR/motifs_DB $OUTPUT_DIR/result $OUTPUT_DIR/tmp \
  --search-type 2 -e 0.01 -s 7.5 -a --threads 2 --remove-tmp-files

# Conversion du résultat au format m8
$MMSEQS convertalis $OUTPUT_DIR/scaffold3_DB $OUTPUT_DIR/motifs_DB $OUTPUT_DIR/result $OUTPUT_DIR/result.m8 \
  --format-output 'query,qlen,tlen,target,pident,alnlen,mismatch,gapopen,qstart,qend,tstart,tend,evalue,bits,qaln,tcov' \
  --search-type 2

