# Unraveling the origins of ichnoviruses: metagenomic discovery of a free-living relative

This repository contains the scripts and supplementary data used for the analyses presented in the manuscript. The scripts do not constitute a single automated pipeline, but rather independent analyses performed throughout the project.

---

# Scripts

## Assembly

- `Snakemake_Trimming_Assembly`: Trim short-reads (Fastp), assemble reads into contigs (Megahit) and perform scaffolding (Redundans python pipeline).

- `Snakemake_BUSCO_Quast_Specimens`: Assess assembly completeness and quality (BUSCO and QUAST).

- `Snakemake_Mapping_Cov`: Map sequencing reads of specimens 1/2 on specimen 1 assembly (reference genome) (Bwa-mem2) and compute genome depth/coverage (Samtools depth and coverage).

---

## Identify candidates

Identification and ORFs detection of candidate viral scaffolds containing IVSPER genes.

### Identification in XwIV specimen 1 assembly

- `Search_IVSPER_in_xiphosomella_wirra.sh`
Search IVSPER genes in the *Xiphosomella wirra* reference genome assembly (specimen 1) (Mmseqs 2 search)

### Identification in XwIV specimen 2 assembly

- `Search_candidate_scaffolds_vs_sp2_genome.sh`: Search candidate scaffolds from specimen 1 against the specimen 2 genome assembly (Mmseqs 2 search)
- `Plot_scaffold3_hits.R`: Visualize homologous regions detected on scaffold 3.
- `Plot_synteny_sp2_sp1_candidate_scaffolds.R`: Generate synteny plots between candidate scaffolds.
- `Align_scaffolds_sp1_sp2.sh`: Align homologous scaffolds between the two specimens (Clustal Omega)
- `Align_test_circularity.sh`: Test whether the reconstructed viral genome is circular.

### Annotate viral candidate of XwIV specimen 2
- `ORFfinder_scaffol3.sh`: Predict ORFs on scaffold 3.
- `Filter_Overlapping_ORFs.py`: Remove overlapping or nested ORFs.
- `parse_orffinder_fasta.py`: Format ORFfinder output files.
- `Plot_filtered_ORFs.py`: Visualize ORF annotations.

---

## Homology search

Homology searches and scaffold3 annotation.

### Searches

- `Homology_search_scaffold3_vs_scaffold3.sh`: Run nt scaffold3 against itself to detect repeated regions (Mmseqs 2 search)
- `Homology_searches_scaffold3_proteins.sh`: Run scaffold3 ORFs vs scaffold3 ORFs search (Mmseqs 2 search)
- `Homology_search_scaffold3_proteins_Refseq_virus.sh`: Search scaffold 3 proteins against RefSeq Virus (Mmseqs 2 search)
- `Homology_search_scaffold3_proteins_IVSPER.sh`: Search scaffold 3 proteins against IVSPER proteins (Mmseqs 2 search)
- `Homology_search_scaffold3_proteins_UniRef50.sh`: Search scaffold 3 proteins against UniRef50 (Mmseqs 2 search)
- `Homology_search_scaffold3_vs_Glypta_segments.sh`: Compare scaffold 3 with Glypta ichnovirus segments (Mmseqs 2 search)
- `Homology_search_scaffold3_vs_PDV_motifs.sh`: Search for known PDV sequence motifs (Mmseqs 2 search)
- `Homology_search_sp2_candidate_scaffolds_vs_sp1_assembly.sh`: Search specimen 2 candidate scaffolds against the specimen 1 assembly.

### Plots

- `Scaffold3_ORFs_plot_annotated.py`: Plot annotated ORFs and associated homology information.
- `plot_alignment.py`: Visualize sequence alignments.
- `Circos_plot_repetitions_annotated.R`: Generate Circos plots showing repeats and annotated genomic features.

---

## Alignment

Alignment of contigs and scaffolds of specimens 1 and 2

- `Align_contigs_scaffolds_ViralGenome.sh`: Align contigs, scaffolds of specimens 1 and 2 against the putative viral genome.
- `Align_contigs_scaffolds_ViralGenome_V2.sh`: Align contigs, scaffolds of specimens 1 and 2 against the putative viral genome V2.
- `Plot_contigs_scaffolds_alignement.py`: Visualize contigs/scaffolds/genome alignments.

---

## Mapping

- `Snakemake_mapping_sp1_sp2_ViralGenome`: Map sequencing reads from both specimens to the putative viral genome (Bwa-mem2 / Samtools depth)

---

## Test_endogenous_exogenous

### Depth and GC analysis

- `Test_depth_sp1_sp2_candidates.py`: For each specimen, the script compares the sequencing depth and GC content of candidate viral scaffolds with BUSCO-containing scaffolds (considered wasp scaffolds). It computes the median sequencing depth and GC content of each candidate scaffold, derives genome-wide distributions from BUSCO scaffolds, and assesses whether candidate values fall within the expected host range based on empirical 2.5-97.5% quantiles. An empirical two-sided p-value is also calculated to quantify how extreme each candidate is relative to the host genome distribution.
- `GC_Cov.R`: generated GC/Cov plots 

### Gene density analysis

- `Gene_density_analysis.sh`: Run coding-density analysis for specimen 1 (ORFfinder)
- `Gene_density_analysis_Sp2.sh`: Run coding-density analysis for specimen 2 (ORFfinder)
- `compute_gene_density.py`: Compute coding density for each scaffold, run by Gene_density_analysis script 
- `analyze_density.R`: Compare candidate scaffold gene density with genome-wide distributions, run by Gene_density_analysis script

---

## Functional

- `Foldseek_multiple_db.sh`: Perform protein structure prediction and homology searches against 3D protein structure databases: Swissprot, UniProt50, BFVD (Foldseek)
- `Size_comparison_ORFs.py`: Compares the lengths of predicted ORFs with those of their homologous proteins identified in reference databases (IVSPER, RefSeq Virus, and UniRef50) to assess the functional status of predicted viral genes.

---

## Phylogeny

- `Make_ORF_clusters.py`: Build homologous protein clusters for phylogenetic analyses.
- `Snakemake_align_phylogeny`: Perform multiple sequence alignment, trimming, and phylogenetic inference (Clustal Omega, Trimal, Iqtree)


# Supplementary data

- `Viral_genome.fa`: putative viral genome

## Predicted ORFs

- `scaffold3_ORFs_filtered.faa`: fasta file
- `scaffold3_ORFs_filtered.tsv`: table
- `Search_Results_IVSPER_vs_Xiphosomella_wirra_sp2.m8`

## Phylogenies

- `Alignments.zip`
- `Phylogenies_pdf/`
