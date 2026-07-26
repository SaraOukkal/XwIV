#!/bin/bash
#SBATCH --job-name=orf_analysis
#SBATCH --output=orf_analysis.out
#SBATCH --error=orf_analysis.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --time=02:00:00
#SBATCH --mem=8G

set -euo pipefail

OUTDIR="/beegfs/data/soukkal/Thesis/Hymenoptera_Project/PDV/Results/Viral_genes/IVSPER_Analysis/Gene_Density_analysis"
#mkdir -p "$OUTDIR"

#echo "1. Nettoyage du FASTA"
INPUT_FASTA="/beegfs/data/soukkal/Thesis/Hymenoptera_Project/Genomes/xiphosomella_wirra.fa"
CLEAN_FASTA="$OUTDIR/xiphosomella_wirra_clean.fa"
#sed 's/|/-/g' "$INPUT_FASTA" > "$CLEAN_FASTA"

#echo "2. ORFfinder (min 70 AA)"
ORF_FASTA="$OUTDIR/orfs_xiphosomella_70AA.faa"
#ORFfinder -in "$CLEAN_FASTA" -out "$ORF_FASTA" -ml 210 -outfmt 0

#echo "3. Parsing du fichier FASTA ORFfinder en TSV"
ORF_TABLE="$OUTDIR/orfs_xiphosomella_70AA.tsv"
#python3 parse_orffinder_fasta.py "$ORF_FASTA" "$ORF_TABLE"

#echo "4. Calcul de la densité codante"
#python3 compute_gene_density.py "$ORF_TABLE" "$CLEAN_FASTA" "$OUTDIR/gene_density_table.tsv"

#echo "5. Extraction des scaffolds BUSCO"
#grep "Complete" "/beegfs/data/soukkal/Thesis/Hymenoptera_Project/Results/Stats/BUSCO/xiphosomella_wirra/run_hymenoptera_odb10/full_table.tsv" | cut -f3 | sort | uniq > "$OUTDIR/list_BUSCO_scaffolds.txt"

echo "6. Analyse R de la distribution"
Rscript analyze_density.R "$OUTDIR/gene_density_table.tsv" "$OUTDIR/list_BUSCO_scaffolds.txt" "$OUTDIR/list_candidate_scaffolds.txt"

