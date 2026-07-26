#!/bin/bash
#SBATCH --job-name=ORFs_vs_UniRef50
#SBATCH --output=%x_%j.out
#SBATCH --error=%x_%j.err
#SBATCH --time=06:00:00
#SBATCH --cpus-per-task=8
#SBATCH --mem=8G
#SBATCH --exclude=pbil-deb[14-27]

MMSEQS="/beegfs/data/soukkal/TOOLS/mmseqs/bin/mmseqs"
ORF_PROTEINS="/beegfs/data/soukkal/Thesis/Hymenoptera_Project/PDV/Results/Viral_genes/IVSPER_Analysis/xiphosomella_wirra_specimen2/ORF_prediction/scaffold3_ORFs_70AA.faa"
UNIRE50="/beegfs/data/soukkal/Thesis/Databases/uniref50.fasta"

OUTPUT_DIR="/beegfs/data/soukkal/Thesis/Hymenoptera_Project/PDV/Results/Viral_genes/IVSPER_Analysis/xiphosomella_wirra_specimen2/Homology_search/ORFs_vs_UniRef50"
mkdir -p $OUTPUT_DIR

$MMSEQS createdb $ORF_PROTEINS $OUTPUT_DIR/scaffold3_ORFs_DB
$MMSEQS createdb $UNIRE50 $OUTPUT_DIR/uniref50_DB

$MMSEQS search $OUTPUT_DIR/scaffold3_ORFs_DB $OUTPUT_DIR/uniref50_DB $OUTPUT_DIR/result $OUTPUT_DIR/tmp \
  -a -s 7.5 -e 0.1 --threads 8 --remove-tmp-files --search-type 3

$MMSEQS convertalis $OUTPUT_DIR/scaffold3_ORFs_DB $OUTPUT_DIR/uniref50_DB $OUTPUT_DIR/result $OUTPUT_DIR/result.m8 \
  --format-output 'query,qlen,tlen,target,pident,alnlen,mismatch,gapopen,qstart,qend,tstart,tend,evalue,bits,qaln,tcov' --search-type 3

