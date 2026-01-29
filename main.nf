#!/usr/bin/env nextflow
nextflow.enable.dsl=2

process FILTER_SNV {
	module 'StdEnv/2023'
	module 'vcftools'
	module 'bcftools'
	input:
		tuple val(sample_id), path(vcf)
	output:
		tuple val(sample_id), path("SNV/*.recode.vcf"), emit: snv_vcf
		tuple val(sample_id), path("INDEL/*.recode.vcf"), emit: indel_vcf
	script:
		"""
		mkdir -p SNV INDEL
		vcftools --vcf ${vcf} --remove-filtered-all --remove-indels --recode --recode-INFO-all --out SNV/${vcf.baseName}.SNV
		vcftools --vcf ${vcf} --remove-filtered-all --keep-only-indels --recode --recode-INFO-all --out INDEL/${vcf.baseName}.INDEL
		"""
}

process ANNOT_VCF {
	module 'StdEnv/2023'
	module 'java/21.0.1'
	input:
		tuple val(sample_id), path(vcf), val(type)
	output:
		tuple val(sample_id),path("${vcf.baseName}.annotMotifPhast.vcf"), val(type)
	publishDir { "${params.output_dir}/${sample_id}/${type}" }, mode: 'copy'

	script:
		"""
		echo "Type: ${type}"
		java -Xmx8g -jar ${params.snpeff_jar}/snpEff.jar -v ${params.snpeff_genome} ${vcf} > ${vcf.baseName}.annotMotif.vcf
		java -Xmx8g -jar ${params.snpeff_jar}/SnpSift.jar phastCons ${params.snpeff_phast} ${vcf.baseName}.annotMotif.vcf > ${vcf.baseName}.annotMotifPhast.vcf
		rm ${vcf.baseName}.annotMotif.vcf
		rm -f ${vcf}
		rm snpEff_summary.html
		rm snpEff_genes.txt
		"""
}

process SIGPROFILER {
	tag "$sample_id:$subdir"
	input:
		tuple val(sample_id), path(annot_dir), val(subdir)
	output:
		tuple val(sample_id), path("sigprofiler"), val(subdir), path("input"), emit: sigprofiler_vcfs
	script:
		"""
		module load StdEnv/2023 python scipy-stack
		source ~/virtenv/bin/activate
		mkdir -p sigprofiler input
		mv ${annot_dir}/*.annotMotifPhast.vcf input
		python ${projectDir}/process/sigprofiler.py input sigprofiler ${params.snpeff_genome} ${subdir}
		"""
}

process SNV_CONTEXT {
	input:
		tuple val(sample_id), path(snv_file), path(decomposed_probs_files), val(subdir)
	output:
		path "*bpSNVSeqcontext.tsv", emit: snv_context
		path "*.bed", emit: snv_bed
	//publishDir { "${params.output_dir}/${sample_id}/${subdir}" }, mode: 'copy'
	script:
		"""
		module load StdEnv/2023 python scipy-stack
		source ~/ENV/bin/activate
		python ${projectDir}/process/mutcontextSNV.py ${snv_file} ${decomposed_probs_files} ${params.kmer_file} ${params.ref_fasta} ${params.flank_size} ${params.chrlist} ${params.flank_size}bpSNVSeqcontext.tsv ${snv_file}.bed
		"""
}

process INDEL_CONTEXT {
	input:
		tuple val(sample_id), path(indel_file), path(decomposed_probs_files), val(subdir)
	output:
		//path "*bpINDELSeqContext.tsv", emit: indel_context
		path "DEL*bpINDELSeqContext.tsv", emit: del_context
		path "INS*bpINDELSeqContext.tsv", emit: ins_context

		path "DEL*.bed", emit: del_bed
		path "INS*.bed", emit: ins_bed
	//publishDir { "${params.output_dir}/${sample_id}/${subdir}" }, mode: 'copy'
	script:
		"""
		module load StdEnv/2023 python scipy-stack
		source ~/ENV/bin/activate
		python ${projectDir}/process/mutcontextINDEL.py ${indel_file} ${decomposed_probs_files} ${params.kmer_file} ${params.ref_fasta} ${params.flank_size} ${params.chrlist} ${params.flank_size}bpINDELSeqContext.tsv ${indel_file}.bed
		"""
}

process MERGE_BED {
	input:
		tuple val(category), path(bed_files)
	output:
		path "merged_5bp${category}Seqcontext.bed"
	publishDir { "${params.output_dir}" }, mode: 'copy'
	script:
		"""
		cat ${bed_files.join(' ')} | sort -k1,1 -k2,2n > merged_${params.flank_size}bp${category}Seqcontext.bed
		"""
}

process MERGE_CONTEXT {
	input:
		tuple val(category), path(context_files)
	output:
		path "*.tsv"
	publishDir { "${params.output_dir}" }, mode: 'copy'
	script:
		"""
		module load StdEnv/2023 python scipy-stack
		source ~/ENV/bin/activate
		python ${projectDir}/process/mergingContext.py merged_${params.flank_size}bp${category}Seqcontext.tsv ${context_files}
		"""
}

workflow {
	if (!params.samplesheet) {error "Please provide a samplesheet using --samplesheet"}

	samples_ch = Channel.fromPath(params.samplesheet).splitCsv(header: true, sep: '\t').map { row ->tuple(row.sample_id,file(row.vcf))}

	filtered = FILTER_SNV(samples_ch)

	snv_vcf = filtered.snv_vcf.map {sid, vcf -> tuple(sid,vcf, 'SNV') }
	indel_vcf = filtered.indel_vcf.map {sid, vcf -> tuple(sid,vcf, 'INDEL') }
	combined_vcfs = snv_vcf.concat(indel_vcf)

	annotated = ANNOT_VCF(combined_vcfs)
	annotated_snv = annotated.filter { it[2] == 'SNV' }
	annotated_indel = annotated.filter { it[2] == 'INDEL' }

	all_annotated = annotated_snv.concat(annotated_indel)
	sigprofiler_ch = SIGPROFILER(all_annotated.map {sid, vcf, type -> tuple(sid, vcf.getParent(), type) })
	snv_sig = sigprofiler_ch.filter { it[2] == 'SNV' }
	indel_sig = sigprofiler_ch.filter { it[2] == 'INDEL' }

	snv_pairs =snv_sig.map {sid, sig_dir, subdir, vcf_dir ->
		tuple(sid, file("${vcf_dir}/*.vcf"),file("${sig_dir}/Assignment_Solution/Activities/Decomposed_Mutation_Probabilities/Decomposed_Mutation_Probabilities_*.txt"), 'SNV')}
	snv_context_result = SNV_CONTEXT(snv_pairs)

	indel_pairs =indel_sig.map {sid, sig_dir, subdir, vcf_dir ->
		tuple(sid, file("${vcf_dir}/*.vcf"),file("${sig_dir}/Assignment_Solution/Activities/Decomposed_Mutation_Probabilities/Decomposed_Mutation_Probabilities_*.txt"), 'INDEL')}
	indel_context_result = INDEL_CONTEXT(indel_pairs)

	all_beds_ch = indel_context_result.del_bed.map { bed -> tuple('DEL', bed) }
		.mix(indel_context_result.ins_bed.map { bed -> tuple('INS', bed) })
		.mix(snv_context_result.snv_bed.map { bed -> tuple('SNV', bed) })
		.groupTuple()
		.map { type, beds -> tuple(type, beds) }
	MERGE_BED(all_beds_ch)

	all_contexts_ch = indel_context_result.del_context.map { context -> tuple('DEL', context) }
		.mix(indel_context_result.ins_context.map { context -> tuple('INS', context) })
		.mix(snv_context_result.snv_bed.map { context -> tuple('SNV', context) })
		.groupTuple()
		.map { type, contexts -> tuple(type, contexts) }
	MERGE_CONTEXT(all_beds_ch)

}


