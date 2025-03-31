#!/usr/bin/env Rscript

library(data.table)
library(dplyr)
library(ggplot2)
library(tidyr)
library(wesanderson)
library(ggrepel)



args <- commandArgs(trailingOnly=TRUE)
load(args[1])

# N ratio plot faceted by lineage
nratio_lineagefacet <- amplicons_missing %>%
	#Limit to samples with less than 5% Ns in consensus sequence (like in Borcard et al. paper)
	#https://www.frontiersin.org/journals/virology/articles/10.3389/fviro.2022.840952/full
	filter(consensus_perc_N < 5) %>% 
	
	#Remove amplicons completely missing
	filter(!is.na(fwd_mismatch)) %>%
	filter(lineage %in% lineages_min) %>%
	
	#Calculate N ratio
	mutate(`N ratio` = count_N/consensus_Ncount) %>%
	
	#Label sample if ratio of amplicon_N to consensus_N >=20%
	mutate(mylabel_sample = case_when(`N ratio`>=0.2 ~ sample_clean, TRUE ~ NA)) %>% 
	ggplot(aes(x = primer_clean, y = `N ratio`)) +
	geom_point(aes(color = lineage, shape = qc_status), size = 2, alpha = 0.8) +
	scale_color_manual(values = c(wes_palette("Darjeeling1"),
															wes_palette("Moonrise2"),
															wes_palette("GrandBudapest1"),
															wes_palette("GrandBudapest2"),
															wes_palette("Darjeeling2"),
															wes_palette("Moonrise3"))) +
	scale_shape_manual(values = c(4,16)) +
	ggrepel::geom_text_repel(aes(label = mylabel_sample), size = 2.5) +
	geom_hline(yintercept = 0.2, linetype = "dashed", color = "navy") +
	facet_wrap(.~lineage, scales = "free_x") +
	ylab("Ratio of N (N_amplicon/N_consensus)") +
	xlab("Amplicon") +
	ggtitle("Ratio of N by lineage") +
	theme_bw() +
	theme(axis.text.x = element_text(angle=45, hjust=1, size=7),
				strip.text.x = element_text(size = 8), plot.title=element_text(hjust=0.5))



# N ratio plot faceted by primer
# Allows us to see if dropouts in amplicons (partial N in amplicons) are specific to lineages 
nratio_primerfacet <- amplicons_missing %>% 
	#Limit to samples with less than 5% Ns in consensus sequence (like in Borcard et al. paper)
	filter(consensus_perc_N<5) %>% 
	
	#filter(run_id==241023002) %>%
	#Remove amplicons completely missing
	filter(!is.na(fwd_mismatch)) %>%
	filter(lineage %in% lineages_min) %>%
	
	#Calculate N ratio
	mutate(`N ratio` = count_N/consensus_Ncount) %>%
	#Label sample if ratio of amplicon_N to consensus_N >=20%
	mutate(mylabel_lineage = case_when(`N ratio`>=0.2 ~ lineage, TRUE ~ NA)) %>% 
	ggplot(aes(x = coverage_depth, y = `N ratio`)) +
	geom_point(aes(color = lineage, shape = qc_status), size = 2, alpha = 0.8) +
	scale_color_manual(values = c(wes_palette("Darjeeling1"), 
															wes_palette("BottleRocket2"), 
															wes_palette("Moonrise2"), 
															wes_palette("GrandBudapest1"), 
															wes_palette("GrandBudapest2"), 
															wes_palette("Darjeeling2"), 
															wes_palette("Moonrise3"))) +
	scale_shape_manual(values = c(4,16)) +
	ggrepel::geom_text_repel(aes(label = mylabel_lineage), size = 2.5) +
	geom_hline(yintercept = 0.2, linetype = "dashed", color = "navy") +
	xlim(0, 100) +
	facet_wrap(.~primer, scales = "free_x") +
	xlab("Genome coverage ≥20x (%)") +
	ylab("N_amplicon/N_consensus") +
	ggtitle("Ratio of N by amplicon") +
	theme_bw() +
	theme(axis.text.x = element_text(angle = 0, size = 8),
				strip.text.x = element_text(size = 8), plot.title = element_text(hjust = 0.5))



ggsave(args[2], nratio_lineagefacet, height = 12, width = 14)
ggsave(args[3], nratio_primerfacet, height = 12, width = 14)

