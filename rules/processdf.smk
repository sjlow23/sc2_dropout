rule process_dfs:
	input:
		rules.process_seqs.output.status,
		stats = rules.process_seqs.output.stats,
		prop = rules.process_seqs.output.prop,
		lineages = INPUT_LINEAGES
		
	output:
		lineagefiltered = OUTDIR + "samples_lineage_filtered.tsv",
		dfs = OUTDIR + "plots/plot_dfs.RData",
		status = OUTDIR + "status/process_dfs.txt"
	params:
		runid = RUNID,
		minlineagecount = MIN_COUNT_LINEAGE,
	conda: "../envs/process_df.yaml"
	shell:
		"""
		csvtk filter2 -t -f '$run_id=={params.runid}' {input.lineages} > {output.lineagefiltered}
		Rscript scripts/process_df.R {input.stats} {output.lineagefiltered} {input.prop} {params.minlineagecount} {output.dfs}

		touch {output.status}
		"""
