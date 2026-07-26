#!/bin/bash
#SBATCH --job-name=MergedGenome_vs_ORFs
#SBATCH --output=%x_%j.out
#SBATCH --error=%x_%j.err
#SBATCH --time=12:00:00
#SBATCH --cpus-per-task=4
#SBATCH --mem=4G
#SBATCH --exclude=pbil-deb[14-27]

# Chemins
MMSEQS="/beegfs/data/soukkal/TOOLS/mmseqs/bin/mmseqs"
QUERY_GENOME="/beegfs/data/soukkal/Thesis/Hymenoptera_Project/PDV/Results/Viral_genes/IVSPER_Analysis/Xiphosomella_wirra_assembly/Merged_xiphosomella_wirra/redundans/scaffolds.reduced.fa"
TARGET_ORFS="/beegfs/data/soukkal/Thesis/Hymenoptera_Project/PDV/Results/Viral_genes/IVSPER_Analysis/xiphosomella_wirra_specimen2/ORF_prediction/scaffold3_ORFs_70AA.faa"
OUTPUT_DIR="/beegfs/data/soukkal/Thesis/Hymenoptera_Project/PDV/Results/Viral_genes/IVSPER_Analysis/xiphosomella_wirra_specimen2/Homology_search/MergedGenome_vs_Scaffold3ORFs"

# Création du dossier de sortie
mkdir -p $OUTPUT_DIR

# Création des bases MMseqs2
$MMSEQS createdb $QUERY_GENOME $OUTPUT_DIR/genome_DB
$MMSEQS createdb $TARGET_ORFS $OUTPUT_DIR/orfs_DB

# Recherche d'homologie (blastx-like)
$MMSEQS search $OUTPUT_DIR/genome_DB $OUTPUT_DIR/orfs_DB $OUTPUT_DIR/result $OUTPUT_DIR/tmp \
  --search-type 1 --threads 4 -s 7.5 -e 0.1 -a --remove-tmp-files

# Conversion des alignements au format .m8
$MMSEQS convertalis $OUTPUT_DIR/genome_DB $OUTPUT_DIR/orfs_DB $OUTPUT_DIR/result $OUTPUT_DIR/result.m8 \
  --format-output 'query,qlen,tlen,target,pident,alnlen,mismatch,gapopen,qstart,qend,tstart,tend,evalue,bits,qaln,taln,tcov' --search-type 1

