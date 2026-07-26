#!/bin/bash
#SBATCH --job-name=CircularityTest_faidx
#SBATCH --output=%x_%j.out
#SBATCH --error=%x_%j.err
#SBATCH --time=00:30:00
#SBATCH --cpus-per-task=2
#SBATCH --mem=2G
#SBATCH --exclude=pbil-deb[14-27]

# Outils
CLUSTALO=/beegfs/home/soukkal/miniconda3/bin/clustalo
SAMTOOLS=/beegfs/data/soft/samtools-1.19/bin/samtools

# Répertoire de sortie
OUTDIR=/beegfs/data/soukkal/Thesis/Hymenoptera_Project/PDV/Results/Viral_genes/IVSPER_Analysis/xiphosomella_wirra_specimen2/circularity_test_faidx
mkdir -p $OUTDIR

# Fichiers FASTA d'origine
SCAF3_FA=/beegfs/data/soukkal/Thesis/Hymenoptera_Project/PDV/Results/Viral_genes/IVSPER_Analysis/xiphosomella_wirra_specimen2/scaffold3.fa
SCAF524_FA=/beegfs/data/soukkal/Thesis/Hymenoptera_Project/PDV/Results/Viral_genes/IVSPER_Analysis/Viral_candidate_scaffolds.fa

# Indexer si nécessaire
$SAMTOOLS faidx $SCAF3_FA
$SAMTOOLS faidx $SCAF524_FA

# Extraire les séquences
$SAMTOOLS faidx $SCAF3_FA "scaffold3|size100523:1-9336" > $OUTDIR/scaffold3_1_9336.fa
$SAMTOOLS faidx $SCAF524_FA "scaffold524|size39189:35095-39189" > $OUTDIR/scaffold524_35095_39189.fa

# Générer le reverse-complement de la région de scaffold524
rc_seq=$($SAMTOOLS faidx $SCAF524_FA "scaffold524|size39189:35095-39189" | \
  awk '!/^>/ {seq = seq $0} END {
    rc = ""; for (i=length(seq); i>0; i--) {
      base = substr(seq, i, 1);
      rc = rc (base=="A"?"T":base=="T"?"A":base=="G"?"C":base=="C"?"G":base)
    }
    print rc
  }')

echo -e ">scaffold524_35095_39189_RC\n$rc_seq" > $OUTDIR/scaffold524_35095_39189_RC.fa

# Alignement direct
$CLUSTALO --infmt fasta -i <(cat $OUTDIR/scaffold3_1_9336.fa $OUTDIR/scaffold524_35095_39189.fa) \
  -o $OUTDIR/alignment_sc3_1_9336_vs_sc524_direct.aln --threads 2

# Alignement reverse-complement
$CLUSTALO --infmt fasta -i <(cat $OUTDIR/scaffold3_1_9336.fa $OUTDIR/scaffold524_35095_39189_RC.fa) \
  -o $OUTDIR/alignment_sc3_1_9336_vs_sc524_RC.aln --threads 2

echo "✅ Alignements terminés :"
echo "- Direct : $OUTDIR/alignment_sc3_1_9336_vs_sc524_direct.aln"
echo "- RC     : $OUTDIR/alignment_sc3_1_9336_vs_sc524_RC.aln"

