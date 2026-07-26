#!/bin/bash
#SBATCH --job-name=ORFs_vs_RefseqVirus_only
#SBATCH --output=%x_%j.out
#SBATCH --error=%x_%j.err
#SBATCH --time=03:00:00
#SBATCH --cpus-per-task=4
#SBATCH --mem=1G
#SBATCH --exclude=pbil-deb[14-27]

# Définir les chemins
MMSEQS="/beegfs/data/soukkal/TOOLS/mmseqs/bin/mmseqs"
ORF_PROTEINS="/beegfs/data/soukkal/Thesis/Hymenoptera_Project/PDV/Results/Viral_genes/IVSPER_Analysis/xiphosomella_wirra_specimen2/ORF_prediction/scaffold3_ORFs_70AA.faa"
REFSEQ_VIRUS="/beegfs/data/soukkal/Thesis/Databases/Refseq_viral_db_nophages_noPDV.faa"

OUTPUT_DIR="/beegfs/data/soukkal/Thesis/Hymenoptera_Project/PDV/Results/Viral_genes/IVSPER_Analysis/xiphosomella_wirra_specimen2/Homology_search/ORFs_vs_RefseqVirus"
mkdir -p $OUTPUT_DIR

# Créer les bases MMseqs2 (dans dossier temporaire)
$MMSEQS createdb $ORF_PROTEINS $OUTPUT_DIR/scaffold3_ORFs_DB
$MMSEQS createdb $REFSEQ_VIRUS $OUTPUT_DIR/refseq_virus_DB

# Effectuer le search
$MMSEQS search $OUTPUT_DIR/scaffold3_ORFs_DB $OUTPUT_DIR/refseq_virus_DB $OUTPUT_DIR/result $OUTPUT_DIR/tmp -a -s 7.5 -e 0.1 --threads 8 --remove-tmp-files --search-type 3

# Convertir les résultats
$MMSEQS convertalis $OUTPUT_DIR/scaffold3_ORFs_DB $OUTPUT_DIR/refseq_virus_DB $OUTPUT_DIR/result $OUTPUT_DIR/result.m8 \
  --format-output 'query,qlen,tlen,target,pident,alnlen,mismatch,gapopen,qstart,qend,tstart,tend,evalue,bits,qaln,tcov' --search-type 3

# Nettoyage temporaire
rm -rf $OUTPUT_DIR/tmp

