#!/usr/bin/env Rscript

library(data.table)
library(dplyr)
library(ggplot2)
library(tidyr)
library(wesanderson)


args <- commandArgs(trailingOnly=TRUE)
load(args[1])
lineages <- fread(args[2], header=T, sep="\t")


amplicons_missing <- amplicons_missing %>%
	mutate(sample = as.character(sample))

lineage_summary <- lineages %>% 
	group_by(lineage) %>% 
	summarize(count_sequences_in_lineage = n_distinct(sample_id)) %>% 
	ungroup() 

barplot.df <- amplicons_missing %>%
	left_join(lineage_summary, by = "lineage") %>%
	#filter(coverage_depth>=80) %>%
	group_by(lineage, count_sequences_in_lineage, primer) %>% 
	summarize(`Amplicon present` = n_distinct(sample[!is.na(fwd_mismatch)])) %>%
	mutate(`Amplicon absent` = count_sequences_in_lineage - `Amplicon present`) %>% 
	pivot_longer(-c(lineage, count_sequences_in_lineage, primer), 
							 names_to = "category", values_to = "count") %>% 
	mutate(primer_clean = as.numeric(gsub("SARS-CoV-2_1200_", "", primer))) %>%
	mutate(primername = as.character(primer_clean)) %>% 
	distinct() 


#Rowcount for plot
rowcount <- ceiling(length(lineages_min)/8)

if (rowcount == 1 ) {
	myheight = 5
} else {
	myheight = 10
}

if (length(lineages_min) >= 5) {
	mywidth = 17
} else {
	mywidth = 12
}

myplot <- barplot.df %>% 
	filter(lineage %in% lineages_min) %>%
	ggplot(aes(x = reorder(primername, primer_clean), y = count)) +
	geom_bar(stat = "identity", position = "fill", aes(fill = category)) + 
	geom_text(aes(label = count), position = position_fill(vjust = 0.5), size = 2) +
	facet_wrap(.~lineage, scales = "free_x", nrow = rowcount) +
	scale_y_continuous(labels = scales::percent) +
	xlab("Primer/Amplicon") +
	ylab("Percentage of genomes") +
	ggtitle("Proportion of missing amplicons by lineage") +
	theme_bw() +
	scale_y_continuous(expand = c(0,0)) +
	theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 7), 
		plot.title = element_text(hjust = 0.5))


ggsave(args[3], myplot, width = mywidth, height = myheight)
