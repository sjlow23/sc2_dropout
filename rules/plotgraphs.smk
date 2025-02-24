rule plot_bargraph:
	input:
		rules.process_dfs.output.status,
		dfs = rules.process_dfs.output.dfs,
		lineages = rules.process_dfs.output.lineagefiltered       
	output:
		barplot = OUTDIR + "plots/barplot.pdf",
		status = OUTDIR + "status/plot_bargraph.txt"
	conda: "../envs/R.yaml"
	shell:
		"""
		Rscript scripts/plot_barplot.R {input.dfs} {input.lineages} {output.barplot}
		touch {output.status}
		"""


rule prepare_heatmaps:
	input:
		rules.process_dfs.output.status,
		dfs = rules.process_dfs.output.dfs,
	output:
		heatmapdfs = OUTDIR + "plots/heatmap_dfs.RData",
		status = OUTDIR + "status/prepare_heatmaps.txt"
	conda: "../envs/R.yaml"
	shell:
		"""
		Rscript scripts/prepare_heatmaps.R {input.dfs} {output.heatmapdfs}
		touch {output.status}
		"""


rule plot_nratio:
	input:
		rules.process_dfs.output.status,
		dfs = rules.process_dfs.output.dfs,
	output:
		nratio_lineage = OUTDIR + "plots/nratio_lineagefacet.pdf",
		nratio_primer = OUTDIR + "plots/nratio_primerfacet.pdf",
		status = OUTDIR + "status/plot_nratio.txt"
	conda: "../envs/R.yaml"
	shell:
		"""
		Rscript scripts/plot_nratio.R {input.dfs} {output.nratio_lineage} {output.nratio_primer}
		touch {output.status}
		"""