#!/usr/bin/env Rscript

# --- Libraries ---
library(ggplot2)
library(dplyr)
library(readr)
library(tools)
library(stringr)

# --- PARAMETERS ---
m8_path <- "/beegfs/data/soukkal/Thesis/Hymenoptera_Project/PDV/Results/Viral_genes/IVSPER_Analysis/xiphosomella_wirra_specimen2/Homology_search/sp2_candidate_scaffolds_vs_sp1_assembly/result.m8"
query_name <- "scaffold3|size100523"

# Correct order: 187 -> 33401 -> 524
ordered_targets <- c("scaffold187|size52917",
                     "scaffold33401|size845",
                     "scaffold524|size39189")

# --- LOAD AND FILTER ---
m8 <- read_tsv(
  m8_path,
  col_names = c("query", "qlen", "tlen", "target", "pident", "alnlen", "mismatch",
                "gapopen", "qstart", "qend", "tstart", "tend", "evalue", "bits",
                "qaln", "tcov"),
  show_col_types = FALSE
)

filtered <- m8 %>%
  filter(query == query_name, target %in% ordered_targets)

qlen <- unique(filtered$qlen)[1]

# --- TARGET LENGTHS + OFFSETS ---
target_lengths <- filtered %>%
  select(target, tlen) %>%
  distinct() %>%
  mutate(target = factor(target, levels = ordered_targets)) %>%
  arrange(target)

target_offsets <- target_lengths %>%
  mutate(offset = lag(cumsum(as.numeric(tlen)), default = 0)) %>%
  mutate(target = as.character(target))

combined_tlen <- sum(as.numeric(target_lengths$tlen))

# --- CONNECTOR SEGMENTS ---
segments <- filtered %>%
  mutate(
    qx1 = pmin(qstart, qend),
    qx2 = pmax(qstart, qend),
    tx1 = pmin(tstart, tend),
    tx2 = pmax(tstart, tend)
  ) %>%
  left_join(target_offsets %>% select(target, offset), by = "target") %>%
  mutate(
    tx1_shift = tx1 + offset,
    tx2_shift = tx2 + offset,
    x_query   = qx1,
    xend_target = (tx1_shift + tx2_shift) / 2,
    y_query   = 0,
    y_target  = 1,
    col_target = case_when(
      target == "scaffold187|size52917" ~ "red",
      target == "scaffold33401|size845" ~ "green",
      target == "scaffold524|size39189" ~ "blue",
      TRUE ~ "grey50"
    ),
    alpha_val = pmin(pmax(pident / 100, 0.3), 1.0)
  )

# --- GEOMETRIES ---
query_bar <- tibble(x = 0, xend = qlen, y = 0)
concat_bar <- tibble(x = 0, xend = combined_tlen, y = 1)

target_blocks <- target_offsets %>%
  mutate(
    x = offset,
    xend = offset + as.numeric(target_lengths$tlen[match(target, target_lengths$target)])
  ) %>%
  mutate(ymin = 1 - 0.06, ymax = 1 + 0.06)

target_labels <- target_offsets %>%
  mutate(
    x = offset + as.numeric(target_lengths$tlen[match(target, target_lengths$target)]) / 2,
    y = 1 + 0.15,
    label = str_c(target)
  )

# --- PLOT ---
p <- ggplot() +
  geom_segment(data = concat_bar,
               aes(x = x, xend = xend, y = y, yend = y),
               linewidth = 2, color = "black") +
  geom_rect(data = target_blocks,
            aes(xmin = x, xmax = xend, ymin = ymin, ymax = ymax),
            fill = "grey80", alpha = 0.4, color = NA) +
  geom_segment(data = query_bar,
               aes(x = x, xend = xend, y = y, yend = y),
               linewidth = 2, color = "black") +
  geom_segment(data = segments,
               aes(x = x_query, xend = xend_target, y = y_query, yend = y_target,
                   color = target, alpha = alpha_val),
               linewidth = 1.1) +
  geom_text(data = target_labels,
            aes(x = x, y = y, label = label), size = 3.5, vjust = 0) +
  annotate("text", x = qlen/2, y = -0.18, label = query_name, size = 4) +
  annotate("text", x = combined_tlen/2, y = 1.22,
           label = "Concatenated targets: 187 | 33401 | 524", size = 4) +
  scale_color_manual(values = c(
    "scaffold187|size52917" = "red",
    "scaffold33401|size845" = "green",
    "scaffold524|size39189" = "blue"
  )) +
  scale_alpha(range = c(0.35, 1.0)) +
  coord_cartesian(xlim = c(0, max(qlen, combined_tlen)), ylim = c(-0.4, 1.4)) +
  theme_minimal(base_size = 12) +
  theme(panel.grid = element_blank(),
        axis.text = element_blank(),
        axis.title = element_blank(),
        axis.ticks = element_blank(),
        legend.position = "none")

# --- SAVE ---
output_path <- file.path(dirname(m8_path), "synteny_concat_plot.svg")
ggsave(output_path, p, width = 12, height = 5)
cat("Plot saved to:", output_path, "\n")
