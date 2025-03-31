INPUT_SEQ = "input/sequences_250311053.fasta"
INPUT_PRIMERS = "input/primers_sc2.bed"
INPUT_LINEAGES = "input/lineage_250311053.tsv"
OUTDIR = "output_250311053/"
RUNID = 250311053
MIN_COUNT_LINEAGE = 3
MAX_PRIMER_MISMATCH = 7		#max mismatches when matching primers



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



