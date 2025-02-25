INPUT_SEQ = "input/sequences_250217041.fasta"
INPUT_PRIMERS = "input/primers_sc2.bed"
INPUT_LINEAGES = "input/lineage_250217041.tsv"
OUTDIR = "output_250217041/"
RUNID = 250217041

rule all:
	input:
		OUTDIR + "status/process.txt",
		OUTDIR + "status/process_dfs.txt",
		OUTDIR + "status/plot_bargraph.txt",
		OUTDIR + "status/prepare_heatmaps.txt",
		OUTDIR + "status/plot_nratio.txt",
		OUTDIR + "status/plotheatmap.txt"
		

include: "rules/stats.smk"
include: "rules/processdf.smk"
include: "rules/plotgraphs.smk"
include: "rules/plotheatmap.smk"



