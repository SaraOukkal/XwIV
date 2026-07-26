#!/bin/bash
#SBATCH --job-name=MMseqs_search_xiphosomella
#SBATCH --output=MMseqs_xiphosomella_%j.out
#SBATCH --error=MMseqs_xiphosomella_%j.err
#SBATCH --cpus-per-task=8
#SBATCH --time=05:00:00
#SBATCH --mem=1G
#SBATCH --exclude=pbil-deb[14-27]
#SBATCH --qos=horizon

# Paths
MMSEQS="/beegfs/data/soukkal/TOOLS/mmseqs/bin/mmseqs"
GENOME="/beegfs/data/soukkal/Thesis/Hymenoptera_Project/Genomes/xiphosomella_wirra.fa"
IVSPER_LIST="/beegfs/data/soukkal/Thesis/Hymenoptera_Project/PDV/Data/IV_domesticated.faa"
OUTPUT_DIR="/beegfs/data/soukkal/Thesis/Hymenoptera_Project/PDV/Results/Viral_genes/IVSPER_Analysis/xiphosomella_wirra/Mmseqs_results"
TMP_DIR="$OUTPUT_DIR/tmp"

# Create output directories if they do not exist
mkdir -p "$OUTPUT_DIR"
mkdir -p "$TMP_DIR"

# Step 2: Create MMseqs databases
echo "Creating MMseqs databases..."
$MMSEQS createdb "$GENOME" "$OUTPUT_DIR/genome_db"
$MMSEQS createdb "$IVSPER_LIST" "$OUTPUT_DIR/IVSPER_db"

# Step 3: Run MMseqs search
echo "Running MMseqs search..."
$MMSEQS search "$OUTPUT_DIR/genome_db" "$OUTPUT_DIR/IVSPER_db" "$OUTPUT_DIR/results" "$TMP_DIR" --threads 8 -e 1e-5 --cov-mode 2 -c 0.5 -a 

# Step 4: Convert results to m8 format
echo "Converting results to m8 format..."
$MMSEQS convertalis "$OUTPUT_DIR/genome_db" "$OUTPUT_DIR/IVSPER_db" "$OUTPUT_DIR/results" "$OUTPUT_DIR/results.m8" --format-output 'query,qlen,target,tlen,pident,alnlen,qstart,qend,tstart,tend,evalue,qaln,tcov'

# Cleanup temporary files
echo "Cleaning up temporary files..."
rm -rf "$TMP_DIR"

echo "MMseqs search completed."

