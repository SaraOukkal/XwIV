#!/bin/bash
#SBATCH --job-name=DHJPAR0036296_vs_candidates
#SBATCH --output=%x_%j.out
#SBATCH --error=%x_%j.err
#SBATCH --time=01:00:00
#SBATCH --cpus-per-task=4
#SBATCH --mem=2G
#SBATCH --exclude=pbil-deb[14-27]

# Définir les chemins vers MMseqs2 et les fichiers d'entrée
MMSEQS="/beegfs/data/soukkal/TOOLS/mmseqs/bin/mmseqs"
GENOME="/beegfs/project/horizon/data/assembly/specimens/DHJPAR0036296/redundans/scaffolds.reduced.fa"
CANDIDATES="/beegfs/data/soukkal/Thesis/Hymenoptera_Project/PDV/Results/Viral_genes/IVSPER_Analysis/Xiphosomella_wirra_candidates_scaffolds.fa"
OUTPUT_DIR="/beegfs/data/soukkal/Thesis/Hymenoptera_Project/PDV/Results/Viral_genes/IVSPER_Analysis/xiphosomella_wirra_specimen2"

# Créer le répertoire de sortie s'il n'existe pas
mkdir -p $OUTPUT_DIR

# Créer les bases de données MMseqs2 pour le génome et les scaffolds candidats
$MMSEQS createdb $GENOME $OUTPUT_DIR/DHJPAR0036296_DB
$MMSEQS createdb $CANDIDATES $OUTPUT_DIR/Xiphosomella_candidates_DB

# Effectuer la recherche d'homologie : genome (query) vs candidates (target)
$MMSEQS search $OUTPUT_DIR/DHJPAR0036296_DB $OUTPUT_DIR/Xiphosomella_candidates_DB $OUTPUT_DIR/DHJPAR0036296_vs_candidates_result $OUTPUT_DIR/tmp -a -s 7.5 -e 0.001 --threads 4 --remove-tmp-files --search-type 3

# Convertir les résultats dans un format lisible
$MMSEQS convertalis $OUTPUT_DIR/DHJPAR0036296_DB $OUTPUT_DIR/Xiphosomella_candidates_DB $OUTPUT_DIR/DHJPAR0036296_vs_candidates_result $OUTPUT_DIR/DHJPAR0036296_vs_candidates_result.m8 --format-output 'query,qlen,tlen,target,pident,alnlen,mismatch,gapopen,qstart,qend,tstart,tend,evalue,bits,qaln,tcov' --search-type 3

# Nettoyage temporaire
rm -rf $OUTPUT_DIR/tmp

