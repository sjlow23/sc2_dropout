#!/usr/bin/env Rscript

library(data.table)
library(dplyr)
library(tidyr)


args <- commandArgs(trailingOnly=TRUE)

amplicons <- fread(args[1], header=T, sep="\t")
lineages <- fread(args[2], header=T, sep="\t")
consensus_nprop <- fread(args[3], header=T, sep="\t")

# Reformat input
lineages <- lineages %>% 
  mutate(sample=paste(sample_id, run_id, sep="|")) %>%
  mutate(sample=as.character(sample))

consensus_nprop <- consensus_nprop %>% 
  mutate(sample=as.character(sample))


# Generate intermediate dfs
#Average amplicons sizes for estimating `N ratio` for missing amplicons

amplicons <- amplicons %>%
  mutate(sample=as.character(sample))

mean_amplicon_size <- amplicons %>% 
  group_by(primer) %>% 
  summarize(mean_amplicon_length=mean(amplicon_size)) %>%
  ungroup()

#Samples with at least one missing amplicon
samples_nonperfect <- amplicons %>%
  left_join(lineages, by="sample") %>%
  group_by(sample, lineage) %>%
  summarize(count=n_distinct(primer)) %>% 
  ungroup() %>%
  filter(count!=27) %>% 
  pull(sample)

#Lineages to plot
lineages_min5 <- amplicons %>% 
  left_join(lineages, by="sample") %>%
  group_by(lineage) %>% 
  summarize(count=n_distinct(sample)) %>%
  filter(count>=5) %>%
  pull(lineage)

#Fill in missing sample/primer combinations
#Lineage and other information
amplicons_missing <- amplicons %>% 
  complete(sample, primer,
           fill=list(amplicon_perc_N=100)) %>% 
  mutate(primer_clean=as.numeric(gsub("SARS-CoV-2_1200_", "", primer))) %>%
  arrange(primer_clean) %>%
  mutate(primer_clean=as.factor(primer_clean)) %>% 
  separate(sample, into=c("sample_clean", "runid"), sep ="\\|", remove=F) %>%
  select(-runid) %>%
  left_join(lineages, by="sample") %>%
  left_join(consensus_nprop, by="sample")

#Samples with no missing or partial amplicons
samples_perfect <- amplicons_missing %>% 
  group_by(sample) %>%
  summarize(count_amp=n_distinct(primer),
            partial_amps=n_distinct(primer[amplicon_perc_N!=0])) %>%
  filter(partial_amps==0 & count_amp==27) %>%
  ungroup()

# Save objects
save(samples_nonperfect, samples_perfect, lineages_min5, amplicons_missing, file = args[4])
