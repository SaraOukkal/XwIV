# === Chargement des packages ===
packages_needed <- c("circlize", "RColorBrewer")
packages_to_install <- packages_needed[!(packages_needed %in% installed.packages()[,"Package"])]
if(length(packages_to_install)) {
  install.packages(packages_to_install, repos = "https://cloud.r-project.org/")
}
suppressPackageStartupMessages(library(circlize))
library(RColorBrewer)

# === Chemins ===
input_dir <- "/beegfs/data/soukkal/Thesis/Hymenoptera_Project/PDV/Results/Viral_genes/IVSPER_Analysis/xiphosomella_wirra_specimen2/"
blast_file <- file.path(input_dir, "Homology_search/scaffold3_self_nucl/result_filtered.m8")
orf_file <- file.path(input_dir, "ORF_prediction/scaffold3_ORFs_70AA_filtered.tsv")
output_file <- file.path(input_dir, "ORF_prediction/Circos_plot_annotated.svg")

# === Constantes ===
scaffold_length <- 100523
scaffold_name <- "scaffold3"

# === Données DRJ et motifs ===
drj_data <- data.frame(
  annotation = c("HdIV_DRJ1L", "HdIV_DRJ1R"),
  start = c(60683, 90615),
  end = c(60754, 90562),
  color = rep("#8B0000", 2)
)

motif_data <- data.frame(
  annotation = rep("Glypta_segment", 4),
  start = c(22549, 31043, 88703, 93788),
  end = c(22475, 31090, 88771, 93744),
  color = rep("#1E90FF", 4)
)

# === Lecture des ORFs filtrés ===
orfs <- read.table(orf_file, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
orf_tracks <- data.frame(
  start = pmin(orfs$start, orfs$end),
  end = pmax(orfs$start, orfs$end),
  strand = orfs$strand,
  color = "#808080" # gris
)

# === Lecture fichier m8 ===
blast <- read.table(blast_file, sep = "\t", header = FALSE, stringsAsFactors = FALSE)
colnames(blast) <- c("query", "qlen", "tlen", "target", "pident", "alnlen", "mismatch", "gapopen",
                     "qstart", "qend", "tstart", "tend", "evalue", "bitscore", "qaln", "tcov")
blast <- blast[!(blast$qstart == blast$tstart & blast$qend == blast$tend), ]

min_q <- pmin(blast$qstart, blast$qend)
max_q <- pmax(blast$qstart, blast$qend)
min_t <- pmin(blast$tstart, blast$tend)
max_t <- pmax(blast$tstart, blast$tend)
id1 <- paste(min_q, max_q, min_t, max_t, sep = "_")
id2 <- paste(min_t, max_t, min_q, max_q, sep = "_")
pair_id <- pmin(id1, id2)
blast$pair_id <- pair_id
blast_unique <- blast[!duplicated(blast$pair_id), ]

links <- data.frame(
  start1 = pmin(blast_unique$qstart, blast_unique$qend),
  end1   = pmax(blast_unique$qstart, blast_unique$qend),
  start2 = pmin(blast_unique$tstart, blast_unique$tend),
  end2   = pmax(blast_unique$tstart, blast_unique$tend)
)

# === Plot SVG ===
svg(output_file, width = 6, height = 6)
circos.clear()
circos.par(start.degree = 90, gap.after = 0)
circos.initialize(factors = scaffold_name, xlim = c(0, scaffold_length))

# Track principal avec axe
circos.track(ylim = c(0, 1), panel.fun = function(x, y) {
  circos.axis(h = "bottom", major.at = seq(0, scaffold_length, by = 10000), labels.cex = 0.5, labels.niceFacing = TRUE)
}, bg.border = NA)

# Liens de répétitions
palette <- colorRampPalette(brewer.pal(8, "Set3"))(nrow(links))
for (i in seq_len(nrow(links))) {
  circos.link(scaffold_name, c(links$start1[i], links$end1[i]),
              scaffold_name, c(links$start2[i], links$end2[i]),
              col = palette[i], border = NA)
}

# === Ajout des ORFs ===
circos.trackPlotRegion(factors = scaffold_name, ylim = c(0, 1), track.height = 0.05, bg.border = NA,
                       panel.fun = function(region, value, ...) {
                         for (i in 1:nrow(orf_tracks)) {
                           circos.rect(
                             xleft = orf_tracks$start[i],
                             xright = orf_tracks$end[i],
                             ybottom = 0,
                             ytop = 1,
                             col = orf_tracks$color[i],
                             border = NA
                           )
                         }
                       })

# === Track DRJ et Glypta motifs ===
circos.trackPlotRegion(factors = scaffold_name, ylim = c(0, 1), track.height = 0.06, bg.border = NA,
                       panel.fun = function(region, value, ...) {
                         for (i in 1:nrow(drj_data)) {
                           circos.rect(
                             xleft = min(drj_data$start[i], drj_data$end[i]),
                             xright = max(drj_data$start[i], drj_data$end[i]),
                             ybottom = 0.55,
                             ytop = 0.85,
                             col = drj_data$color[i],
                             border = NA
                           )
                         }
                         for (i in 1:nrow(motif_data)) {
                           circos.rect(
                             xleft = min(motif_data$start[i], motif_data$end[i]),
                             xright = max(motif_data$start[i], motif_data$end[i]),
                             ybottom = 0.15,
                             ytop = 0.45,
                             col = motif_data$color[i],
                             border = NA
                           )
                         }
                       })

# === Légende ===
legend("bottomleft",
       legend = c("Self-hits", "ORFs (filtered)", "DRJ", "Glypta motifs"),
       fill = c("grey", "#808080", "#8B0000", "#1E90FF"),
       border = NA, cex = 0.7)

dev.off()
