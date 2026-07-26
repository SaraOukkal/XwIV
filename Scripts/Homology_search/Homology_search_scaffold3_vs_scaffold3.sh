#!/bin/bash
#SBATCH --job-name=scaffold3_self_nucl
#SBATCH --output=%x_%j.out
#SBATCH --error=%x_%j.err
#SBATCH --time=01:00:00
#SBATCH --cpus-per-task=2
#SBATCH --mem=10G
#SBATCH --exclude=pbil-deb[14-27]

MMSEQS="/beegfs/data/soukkal/TOOLS/mmseqs/bin/mmseqs"
SCAFFOLD="/beegfs/data/soukkal/Thesis/Hymenoptera_Project/PDV/Results/Viral_genes/IVSPER_Analysis/xiphosomella_wirra_specimen2/scaffold3.fa"
OUTPUT_DIR="/beegfs/data/soukkal/Thesis/Hymenoptera_Project/PDV/Results/Viral_genes/IVSPER_Analysis/xiphosomella_wirra_specimen2/Homology_search/scaffold3_self_nucl"
#mkdir -p $OUTPUT_DIR

# Create DB
#$MMSEQS createdb $SCAFFOLD $OUTPUT_DIR/scaffold3_DB

# Self search (nucleotide vs nucleotide)
#$MMSEQS search $OUTPUT_DIR/scaffold3_DB $OUTPUT_DIR/scaffold3_DB $OUTPUT_DIR/result $OUTPUT_DIR/tmp \
#  --search-type 2 -e 0.01 -s 7.5 -a --threads 2 --remove-tmp-files

# Convert to m8 format
#$MMSEQS convertalis $OUTPUT_DIR/scaffold3_DB $OUTPUT_DIR/scaffold3_DB $OUTPUT_DIR/result $OUTPUT_DIR/result.m8 \
#  --format-output 'query,qlen,tlen,target,pident,alnlen,mismatch,gapopen,qstart,qend,tstart,tend,evalue,bits,qaln,tcov' \
#  --search-type 2

# Step 4: Filter out exact self-alignments (based on coordinates)
awk '($9 != $11 || $10 != $12)' $OUTPUT_DIR/result.m8 > $OUTPUT_DIR/result_filtered.m8

