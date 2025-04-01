#!/usr/bin/env Rscript

library(dplyr)
library(data.table)
library(pheatmap)
library(wesanderson)
library(RColorBrewer)

args <- commandArgs(trailingOnly=TRUE)
load(args[1])

mycolors1 <- c("#FFD700", "#FF0000","#2BE4E9", "#F19CBB", "#9C7C38", "#FF7540", "#107090",  "#D8BFD8", "#D5006A",
              "#B1A2CA","#EFBE7D","#61F4DE","#FFCACA", "#C0C0C0", "#676767","#FCF5C7","#000000",  "#87CEFA", 
              "#E9EC6B", "#BA55D3", "#FF4500","#FF1493", "#4D4DFF", "#FF00FF", 
              "#2F4F4F","#008000","#008080","#6B8E23","#00FF00","#00FF7F","#ADFF2F", "#acff80")

mycolors2 <- c(wes_palette("Darjeeling1"), wes_palette("FantasticFox1"), wes_palette("Chevalier1"),
               wes_palette("Royal2"), wes_palette("Moonrise3"),
               wes_palette("Zissou1"))

mycolorpalette <- c(mycolors1, mycolors2)


# Amplicon heatmap
heatmap_amp_annotation <- fread(args[4], header = T, sep = "\t") %>%
	select(sample_id, lineage, completeness, coverage_depth, mean_depth, qc_status) %>%
	filter(sample_id %in% rownames(heatmap_amp.mat)) 
heatmap_amp_annotation <- tibble::column_to_rownames(heatmap_amp_annotation, var = "sample_id")

# Generate color palette for lineages
lineages <- unique(heatmap_amp_annotation$lineage)
lineages_colors <- mycolorpalette[1:length(lineages)]
named_colors_lineages <- setNames(lineages_colors, lineages)

ann_colors <- list(
  lineage = named_colors_lineages,
  qc_status = c("pass" = "green", "fail" = "red"),
  mean_depth = brewer.pal(8, 'PuBu'),
  completeness = brewer.pal(8, 'BuGn'),
  coverage_depth = brewer.pal(8, 'PuRd')
)


pdf(file=args[2], width=15, height=15)
#Heatmap of amplicon dropout
pheatmap(heatmap_amp.mat, 
		cluster_cols = F, 
		cluster_rows = T, 
		scale = "none", 
		angle_col = 90, 
		na_col = "grey90",
		fontsize = 6, 
		display_numbers = F, 
		cellwidth = 11, 
		cellheight = 10,
		annotation_row = heatmap_amp_annotation,
		annotation_colors = ann_colors,
		main = "Proportion of Ns in amplicon (based on consensus sequence)")
dev.off()


# Primer mismatch heatmap
pdf(file = args[3], width = 18, height = 16)
pheatmap(heatmap_primer.mat, 
		cluster_cols = F, 
		cluster_rows = F, 
		scale = "none", 
		angle_col = 90, 
		na_col = "grey90",
		fontsize = 6, 
		display_numbers = F, 
		cellwidth = 12, 
		cellheight = 10.5,
		annotation_row = heatmap_amp_annotation,
		annotation_colors = ann_colors,
		main = "Mismatches in SARS-CoV-2 primer binding sites (based on consensus sequence)")
dev.off()

