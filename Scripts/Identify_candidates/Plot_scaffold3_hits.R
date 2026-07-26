# Packages
library(ggplot2)
library(dplyr)
library(readr)

cat("📥 Lecture du fichier m8...\n")
output_dir <- "/beegfs/data/soukkal/Thesis/Hymenoptera_Project/PDV/Results/Viral_genes/IVSPER_Analysis/xiphosomella_wirra_specimen2/"
m8_path <- file.path(output_dir, "DHJPAR0036296_vs_candidates_result.m8")

# Lire le fichier
m8 <- read_tsv(m8_path,
               col_names = c("query", "qlen", "tlen", "target", "pident", "alnlen",
                             "mismatch", "gapopen", "qstart", "qend", "tstart", "tend",
                             "evalue", "bits", "qaln", "tcov"),
               col_types = cols(.default = "c"))

### === PARTIE 1 : alignements sur scaffold3|size100523 ===

cat("🔍 Recherche des lignes où 'query' est scaffold3|size100523...\n")
scaffold_query <- "scaffold3|size100523"

matching_queries <- m8 %>%
  filter(query == scaffold_query)

cat("✅ Nombre de lignes trouvées :", nrow(matching_queries), "\n")

if (nrow(matching_queries) == 0) {
  cat("⚠️ Aucun alignement trouvé pour scaffold3|size100523 en tant que query.\n")
  quit(save = "no")
}

# Conversion et positions
hits <- matching_queries %>%
  mutate(
    qstart = as.integer(qstart),
    qend = as.integer(qend),
    start = pmin(qstart, qend),
    end = pmax(qstart, qend)
  )

cat("📌 Exemple de positions d'alignement :\n")
print(head(hits[, c("target", "start", "end")]))

# Taille totale du scaffold
scaffold_length <- 100523

# Plot des positions sur scaffold3
cat("🖼️ Création du plot d'alignements...\n")
plot1 <- ggplot(hits, aes(xmin = start, xmax = end, y = 1, fill = target)) +
  geom_rect(aes(ymin = 0.6, ymax = 1.4)) +
  scale_x_continuous(limits = c(0, scaffold_length), expand = c(0, 0)) +
  labs(
    title = paste("Alignements sur", scaffold_query),
    x = "Position sur scaffold3 (query)",
    y = NULL,
    fill = "Target (scaffold)"
  ) +
  theme_minimal() +
  theme(axis.text.y = element_blank(), axis.ticks.y = element_blank())

ggsave(file.path(output_dir, "alignements_query_scaffold3.png"), plot1, width = 12, height = 2, dpi = 300)
cat("✅ Plot 1 sauvegardé.\n")

### === PARTIE 2 : couverture des scaffolds viraux par scaffold3 ===

cat("📊 Calcul de la couverture des scaffolds viraux par scaffold3...\n")

# Définir les scaffolds candidats et leurs tailles
scaffold_sizes <- c(
  "scaffold524|size39189" = 39189,
  "scaffold12966|size5621" = 5621,
  "scaffold9801|size7733" = 7733,
  "scaffold33401|size845" = 845,
  "scaffold187|size52917" = 52917,
  "scaffold21569|size2487" = 2487
)

# Extraire les hits où target est un des scaffolds candidats et query == scaffold3
coverage_hits <- m8 %>%
  filter(query == scaffold_query & target %in% names(scaffold_sizes)) %>%
  mutate(
    tstart = as.integer(tstart),
    tend = as.integer(tend),
    start = pmin(tstart, tend),
    end = pmax(tstart, tend)
  )

# Calculer les % de couverture (en fusionnant les intervalles)
library(IRanges)

coverage_summary <- coverage_hits %>%
  group_by(target) %>%
  summarise(
    total_covered = sum(reduce(IRanges(start, end))@width),
    scaffold_size = scaffold_sizes[unique(target)],
    coverage_pct = round(100 * total_covered / scaffold_size, 2)
  ) %>%
  arrange(desc(coverage_pct))

cat("📈 Couverture estimée (%):\n")
print(coverage_summary)

# Plot couverture
plot2 <- ggplot(coverage_summary, aes(x = reorder(target, -coverage_pct), y = coverage_pct)) +
  geom_col(fill = "#377EB8") +
  geom_text(aes(label = paste0(coverage_pct, "%")), vjust = -0.5, size = 3.5) +
  ylim(0, 110) +
  labs(
    title = "Couverture des scaffolds candidats par scaffold3",
    x = "Scaffold candidat (target)",
    y = "Pourcentage de couverture"
  ) +
  theme_minimal()

ggsave(file.path(output_dir, "couverture_par_scaffold3.png"), plot2, width = 8, height = 4, dpi = 300)
cat("✅ Plot 2 sauvegardé.\n")

