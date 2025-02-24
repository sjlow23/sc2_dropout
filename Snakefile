OUTDIR = "my_output/"
RUNID = 241217059

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



