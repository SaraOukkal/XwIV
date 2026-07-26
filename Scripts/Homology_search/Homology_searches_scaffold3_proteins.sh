#!/bin/bash
#SBATCH --job-name=scaffold3_mmseqs
#SBATCH --output=%x_%j.out
#SBATCH --error=%x_%j.err
#SBATCH --time=06:00:00
#SBATCH --cpus-per-task=4
#SBATCH --mem=1G
#SBATCH --exclude=pbil-deb[14-27]

# Définir les chemins
MMSEQS="/beegfs/data/soukkal/TOOLS/mmseqs/bin/mmseqs"

ORF_PROTEINS="/beegfs/data/soukkal/Thesis/Hymenoptera_Project/PDV/Results/Viral_genes/IVSPER_Analysis/xiphosomella_wirra_specimen2/ORF_prediction/scaffold3_ORFs_70AA.faa"
XIPHOSOMELLA_SP1="/beegfs/project/horizon/data/assembly/species/xiphosomella_wirra/gnm.fna"
XIPHOSOMELLA_SP2="/beegfs/project/horizon/data/assembly/specimens/DHJPAR0036296/redundans/scaffolds.reduced.fa"

BASE_OUTPUT="/beegfs/data/soukkal/Thesis/Hymenoptera_Project/PDV/Results/Viral_genes/IVSPER_Analysis/xiphosomella_wirra_specimen2/Homology_search"
mkdir -p $BASE_OUTPUT

# Créer les bases de données (répertoire commun)
DB_DIR="$BASE_OUTPUT/MMseqs_DBs"
mkdir -p $DB_DIR
$MMSEQS createdb $ORF_PROTEINS $DB_DIR/scaffold3_ORFs_DB
$MMSEQS createdb $XIPHOSOMELLA_SP1 $DB_DIR/xiphosomella_sp1_DB
$MMSEQS createdb $XIPHOSOMELLA_SP2 $DB_DIR/xiphosomella_sp2_DB
$MMSEQS createdb $REFSEQ_VIRUS $DB_DIR/refseq_virus_DB

###########################################
# 1. ORFs vs ORFs (paralogues internes)
###########################################
OUT1="$BASE_OUTPUT/ORFs_vs_ORFs"
mkdir -p $OUT1
$MMSEQS search $DB_DIR/scaffold3_ORFs_DB $DB_DIR/scaffold3_ORFs_DB $OUT1/result $OUT1/tmp -a -s 7.5 -e 0.001 --threads 8 --remove-tmp-files --search-type 3
$MMSEQS convertalis $DB_DIR/scaffold3_ORFs_DB $DB_DIR/scaffold3_ORFs_DB $OUT1/result $OUT1/result.m8 --format-output 'query,qlen,tlen,target,pident,alnlen,mismatch,gapopen,qstart,qend,tstart,tend,evalue,bits,qaln,tcov' --search-type 3

###########################################
# 2. Genome sp1 vs ORFs
###########################################
OUT2="$BASE_OUTPUT/GenomeSp1_vs_ORFs"
mkdir -p $OUT2
$MMSEQS search $DB_DIR/xiphosomella_sp1_DB $DB_DIR/scaffold3_ORFs_DB $OUT2/result $OUT2/tmp -a -s 7.5 -e 0.001 --threads 8 --remove-tmp-files --search-type 3
$MMSEQS convertalis $DB_DIR/xiphosomella_sp1_DB $DB_DIR/scaffold3_ORFs_DB $OUT2/result $OUT2/result.m8 --format-output 'query,qlen,tlen,target,pident,alnlen,mismatch,gapopen,qstart,qend,tstart,tend,evalue,bits,qaln,tcov' --search-type 3

###########################################
# 3. Genome sp2 vs ORFs
###########################################
OUT3="$BASE_OUTPUT/GenomeSp2_vs_ORFs"
mkdir -p $OUT3
$MMSEQS search $DB_DIR/xiphosomella_sp2_DB $DB_DIR/scaffold3_ORFs_DB $OUT3/result $OUT3/tmp -a -s 7.5 -e 0.001 --threads 8 --remove-tmp-files --search-type 3
$MMSEQS convertalis $DB_DIR/xiphosomella_sp2_DB $DB_DIR/scaffold3_ORFs_DB $OUT3/result $OUT3/result.m8 --format-output 'query,qlen,tlen,target,pident,alnlen,mismatch,gapopen,qstart,qend,tstart,tend,evalue,bits,qaln,tcov' --search-type 3


# Nettoyage temporaire
rm -rf $OUT1/tmp $OUT2/tmp $OUT3/tmp

