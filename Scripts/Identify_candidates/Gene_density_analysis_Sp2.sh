#!/bin/bash
#SBATCH --job-name=gene_density_sp2
#SBATCH --output=gene_density_sp2.out
#SBATCH --error=gene_density_sp2.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --time=02:00:00
#SBATCH --mem=8G

set -euo pipefail

# === Dossier de sortie ===
OUTDIR="/beegfs/data/soukkal/Thesis/Hymenoptera_Project/PDV/Results/Viral_genes/IVSPER_Analysis/xiphosomella_wirra_specimen2/Gene_density"
mkdir -p "$OUTDIR"

# === Fichiers d’entrée ===
INPUT_FASTA="/beegfs/project/horizon/data/assembly/specimens/DHJPAR0036296/redundans/scaffolds.reduced.fa"
CLEAN_FASTA="$OUTDIR/genome_DHJPAR0036296_clean.fa"
ORF_FASTA="$OUTDIR/orfs_DHJPAR0036296_70AA.faa"
ORF_TABLE="$OUTDIR/orfs_DHJPAR0036296_70AA.tsv"
DENSITY_TABLE="$OUTDIR/gene_density_table.tsv"
BUSCO_FULL="/beegfs/project/horizon/data/stats/busco/specimens/DHJPAR0036296/run_insecta_odb10/full_table.tsv"
BUSCO_LIST="$OUTDIR/list_BUSCO_scaffolds.txt"
CANDIDATE_LIST="$OUTDIR/list_candidate_scaffolds.txt"  # Déjà existant

# Étape 1 : Nettoyage du FASTA
#echo "1. Nettoyage du FASTA"
#sed 's/|/-/g' "$INPUT_FASTA" > "$CLEAN_FASTA"

# Étape 2 : ORFfinder (min 70 AA)
#echo "2. ORFfinder (min 70 AA)"
#/beegfs/home/soukkal/Tools/ORFfinder -in "$CLEAN_FASTA" -out "$ORF_FASTA" -ml 210 -outfmt 0

# Étape 3 : Parsing FASTA ORFfinder
#echo "3. Parsing du fichier FASTA ORFfinder en TSV"
#python3 parse_orffinder_fasta.py "$ORF_FASTA" "$ORF_TABLE"

# Étape 4 : Calcul de la densité codante
#echo "4. Calcul de la densité codante"
#python3 compute_gene_density.py "$ORF_TABLE" "$CLEAN_FASTA" "$DENSITY_TABLE"

# Étape 5 : Extraction des scaffolds BUSCO
#echo "5. Extraction des scaffolds BUSCO"
#grep "Complete" "$BUSCO_FULL" | cut -f3 | sort | uniq > "$BUSCO_LIST"

# Étape 6 : Analyse R
echo "6. Analyse R de la distribution"
Rscript analyze_density.R "$DENSITY_TABLE" "$BUSCO_LIST" "$CANDIDATE_LIST"

