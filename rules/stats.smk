rule process_seqs:
	input:
		seq = INPUT_SEQ,
		primers = INPUT_PRIMERS,
		lineages = INPUT_LINEAGES
	output:
		filteredseq = OUTDIR + "all_consensus_filtered.fasta",
		prop = OUTDIR + "consensus_Nprop.tsv",
		ispcr = OUTDIR + "primers_ispcr.tsv",
		amplicon = OUTDIR + "all_consensus.amplicon.tsv",
		stats = OUTDIR + "all_consensus.stats.tsv",
		status = OUTDIR + "status/process.txt"
	conda: "../envs/process.yaml"
	params:
		runid = RUNID,
		mismatch = 5
	threads: 32
	shell:
		"""
		csvtk filter2 -t -f '$run_id=={params.runid}' {input.lineages} | awk -F "\\t" '{{ print $2"|"$1 }}' > tmp
		
		seqkit grep -n -f tmp {input.seq} -o {output.filteredseq}
		#rm tmp

		seqtk comp {output.filteredseq} | \
		cut -f1,2,9 | \
		awk '{{ print $0, $3/$2*100 }}' OFS="\\t" > {output.prop}

		echo -e "sample\\tlength\\tconsensus_Ncount\\tconsensus_perc_N\\n$(cat {output.prop})" > {output.prop}

		paste - - < {input.primers} | \
		cut -f1,3,6 > {output.ispcr}

		seqkit amplicon -p {output.ispcr} \
		--max-mismatch {params.mismatch} \
		-j {threads} \
		--output-mismatches \
		--bed \
		{output.filteredseq} > {output.amplicon}

		awk 'BEGIN {{ FS=OFS="\\t" }} {{ print $0, gsub(/N/,"",$7) }}' OFS="\\t" {output.amplicon} | \
		awk '{{ print $0, length($7) }}' OFS="\t" | \
		awk -F "\\t" '{{ print $1, $4, $9, $10, $11, $12, $11/$12*100 }}' OFS="\\t" > {output.stats}
		
		echo -e "sample\\tprimer\\tfwd_mismatch\\trev_mismatch\\tcount_N\\tamplicon_size\\tamplicon_perc_N\\n$(cat {output.stats})" > {output.stats}

		touch {output.status}
		"""
