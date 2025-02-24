#!/usr/bin/env Rscript

library(dplyr)
library(data.table)
library(pheatmap)

args <- commandArgs(trailingOnly=TRUE)
load(args[1])


# Amplicon heatmap
heatmap_amp_annotation <- fread(args[4], header=T, sep="\t") %>%
	select(sample_id, lineage, completeness, coverage_depth, mean_depth) %>%
	filter(sample_id %in% rownames(heatmap_amp.mat)) 
heatmap_amp_annotation <- tibble::column_to_rownames(heatmap_amp_annotation, var="sample_id")

pdf(file=args[2], width=12, height=15)
#Heatmap of amplicon dropout
pheatmap(heatmap_amp.mat, cluster_cols=F, cluster_rows=T, 
		 scale = "none", angle_col=90, na_col="grey90",
		 fontsize=6, display_numbers=F, 
		 cellwidth=11, cellheight=10,
		 annotation_row=heatmap_amp_annotation,
		 main="Proportion of Ns in amplicon")
dev.off()



# Primer mismatch heatmap
pdf(file = args[3], width=12, height=15)
pheatmap(heatmap_primer.mat, cluster_cols=F, cluster_rows=F, 
		 scale = "none", angle_col=90, na_col="grey90",
		 fontsize=6, display_numbers=F, 
		 cellwidth=11, cellheight=10,
		 annotation_row=heatmap_amp_annotation,
		 main="Mismatches in SARS-CoV-2 primer binding sites")
dev.off()

