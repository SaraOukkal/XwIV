#!/bin/bash
#SBATCH --job-name=ClustalO_with_RC
#SBATCH --output=%x_%j.out
#SBATCH --error=%x_%j.err
#SBATCH --time=01:00:00
#SBATCH --cpus-per-task=4
#SBATCH --mem=4G
#SBATCH --exclude=pbil-deb[14-27]

# Outil Clustal Omega
CLUSTAL_O=/beegfs/home/soukkal/miniconda3/bin/clustalo

# Fichiers
SCAF3=/beegfs/data/soukkal/Thesis/Hymenoptera_Project/PDV/Results/Viral_genes/IVSPER_Analysis/xiphosomella_wirra_specimen2/scaffold3.fa
CANDIDATES=/beegfs/data/soukkal/Thesis/Hymenoptera_Project/PDV/Results/Viral_genes/IVSPER_Analysis/Viral_candidate_scaffolds.fa
OUTDIR=/beegfs/data/soukkal/Thesis/Hymenoptera_Project/PDV/Results/Viral_genes/IVSPER_Analysis/xiphosomella_wirra_specimen2/clustalo_alignment
mkdir -p $OUTDIR

# Extraire les deux scaffolds normaux
awk '/^>scaffold524\|size39189|^>scaffold187\|size52917/{f=1; print; next} /^>/{f=0} f' $CANDIDATES > $OUTDIR/scaffolds_normal.fa

# Extraire et reverse-complement le scaffold33401
awk '/^>scaffold33401\|size845/{f=1; h=$0; next} /^>/{f=0} f {seq=seq $0} END {
  rc=""; for(i=length(seq);i>0;i--) {
    base=substr(seq,i,1);
    rc=rc (base=="A"?"T":base=="T"?"A":base=="G"?"C":base=="C"?"G":base)
  }
  print h "_RC"; print rc
}' $CANDIDATES > $OUTDIR/scaffold33401_RC.fa

# Concaténer tous les fichiers
INPUT=$OUTDIR/combined_with_rc.fa
cat $SCAF3 $OUTDIR/scaffolds_normal.fa $OUTDIR/scaffold33401_RC.fa > $INPUT

# Lancer Clustal Omega
$CLUSTAL_O -i $INPUT -o $OUTDIR/alignment_with_rc.aln --threads 4

echo "✅ Alignement terminé : $OUTDIR/alignment_with_rc.aln"

