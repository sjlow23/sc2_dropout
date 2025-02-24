rule plot_heatmap:
	input:
		rules.prepare_heatmaps.output.status,
		rules.prepare_heatmaps.output.heatmapdfs,
		dfs = rules.prepare_heatmaps.output.heatmapdfs,
		annot = rules.process_dfs.output.lineagefiltered,
	output:
		heatmap_amp = OUTDIR + "plots/heatmap_amplicons.pdf",
		heatmap_primer = OUTDIR + "plots/heatmap_primers_mismatch.pdf",
		status = OUTDIR + "status/plotheatmap.txt"
	conda: "../envs/pheatmap.yaml"
	shell:
		"""
		
		Rscript scripts/plot_heatmap.R {input.dfs} {output.heatmap_amp} {output.heatmap_primer} {input.annot}
		
		touch {output.status}
		"""