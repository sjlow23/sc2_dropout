#!/usr/bin/env Rscript

library(data.table)
library(dplyr)
library(ggplot2)
library(tidyr)
library(wesanderson)



args <- commandArgs(trailingOnly=TRUE)
load(args[1])


#Function for saving pheatmap
save_pheatmap_pdf <- function(x, filename, width=10, height=12) {
   stopifnot(!missing(x))
   stopifnot(!missing(filename))
   pdf(filename, width=width, height=height)
   grid::grid.newpage()
   grid::grid.draw(x$gtable)
   dev.off()
}

# Amplicon dropout (Ns) heatmap
heatmap_amp.df <- amplicons_missing %>% 
  filter(!sample %in% samples_perfect) %>% 
  #filter(run_id==241023002 & completeness>=90) %>%
  arrange(lineage) %>%
  select(sample_clean, primer_clean, amplicon_perc_N, consensus_perc_N, lineage) %>%
  as.data.frame()

heatmap_amp.mat <- heatmap_amp.df %>%
  select(-consensus_perc_N, -lineage) %>%
  pivot_wider(names_from=primer_clean, values_from=amplicon_perc_N) %>%
  as.data.frame()

rownames(heatmap_amp.mat) <- heatmap_amp.mat[, 1]
heatmap_amp.mat <- heatmap_amp.mat[, -1]

#Annotations
heatmap_amp_annotation <- amplicons_missing %>% 
  filter(sample_clean %in% rownames(heatmap_amp.mat)) %>% 
  select(sample_clean, lineage, qc_status, completeness, coverage_depth) %>%
  distinct() %>%
  as.data.frame()
rownames(heatmap_amp_annotation) <- heatmap_amp_annotation[, 1]
heatmap_amp_annotation <- heatmap_amp_annotation[, -1]



# Primer mismatches heatmap
heatmap_primer.df <- amplicons_missing %>% 
  filter(!sample %in% samples_perfect) %>% 
  #filter(run_id==241023002 & completeness>=90) %>%
  arrange(lineage) %>%
  select(sample_clean, primer_clean, fwd_mismatch, rev_mismatch, consensus_perc_N, lineage) %>%
  as.data.frame()

heatmap_primer.mat <- heatmap_primer.df %>%
  select(-consensus_perc_N, -lineage) %>%
  pivot_longer(-c(sample_clean, primer_clean), names_to="position", values_to="mismatch") %>%
  mutate(primer_final=case_when(position=="fwd_mismatch" ~ paste0(primer_clean, "F"), 
                          TRUE ~ paste0(primer_clean, "R"))) %>%
  select(-primer_clean, -position) %>%
  pivot_wider(names_from=primer_final, values_from=mismatch) %>%
  as.data.frame()

rownames(heatmap_primer.mat) <- heatmap_primer.mat[, 1]
heatmap_primer.mat <- heatmap_primer.mat[, -1]

#Annotations
heatmap_primer_annotation <- amplicons_missing %>% 
  filter(sample_clean %in% rownames(heatmap_primer.mat)) %>% 
  select(sample_clean, lineage, qc_status, completeness, coverage_depth) %>%
  distinct() %>%
  as.data.frame()
rownames(heatmap_primer_annotation) <- heatmap_primer_annotation[, 1]
heatmap_primer_annotation <- heatmap_primer_annotation[, -1]

save(heatmap_amp.mat, heatmap_primer.mat, file = args[2])
