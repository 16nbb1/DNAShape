#!/usr/bin/env nextflow
nextflow.enable.dsl=2

process FILTER_SNV {
	module 'StdEnv/2023'
	module 'vcftools'
	module 'bcftools'
	input:
		path vcf
	output:
		path "SNV/*.recode.vcf", emit: snv_vcf
		path "INDEL/*.recode.vcf", emit: indel_vcf
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
		tuple path(vcf), val(type)
	output:
		tuple path("${vcf.baseName}.annotMotifPhast.vcf"), val(type)
	publishDir { "${params.output_dir}/${type}" }, mode: 'copy'

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
	tag "$subdir"
	input:
		tuple val(annot_dir), val(subdir)
	output:
		tuple path("sigprofiler"), val(subdir), path("input"), emit: sigprofiler_vcfs

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
		tuple path(snv_file), path(decomposed_probs_files), val(subdir)
	output:
		path "*.tsv"
		path "*.bed"
	publishDir { "${params.output_dir}/${subdir}" }, mode: 'copy'
	script:
		"""
		module load StdEnv/2023 python scipy-stack
		source ~/ENV/bin/activate
		python ${projectDir}/process/mutcontextSNV.py ${snv_file} ${decomposed_probs_files} ${params.kmer_file} ${params.ref_fasta} ${params.flank_size} ${params.chrlist} ${params.flank_size}bpSNVSeqcontext.tsv ${snv_file}.bed
		"""
}

process INDEL_CONTEXT {
	input:
		tuple path(indel_file), path(decomposed_probs_files), val(subdir)
	output:
		path "*.tsv"
		path "*.bed"

	publishDir { "${params.output_dir}/${subdir}" }, mode: 'copy'
	script:
		"""
		module load StdEnv/2023 python scipy-stack
		source ~/ENV/bin/activate
		python ${projectDir}/process/mutcontextINDEL.py ${indel_file} ${decomposed_probs_files} ${params.kmer_file} ${params.ref_fasta} ${params.flank_size} ${params.chrlist} ${params.flank_size}bpINDELContext.tsv ${indel_file}.bed
		"""
}


workflow {
	if (!params.input) {error "Please provide an input file using --input"}

	input_ch = Channel.fromPath(params.input)
	filtered = FILTER_SNV(input_ch)

	snv_vcf = filtered.snv_vcf.map { vcf -> tuple(vcf, 'SNV') }
	indel_vcf = filtered.indel_vcf.map { vcf -> tuple(vcf, 'INDEL') }
	combined_vcfs = snv_vcf.concat(indel_vcf)

	annotated = ANNOT_VCF(combined_vcfs)
	annotated_snv = annotated.filter { it[1] == 'SNV' }
	annotated_indel = annotated.filter { it[1] == 'INDEL' }

	all_annotated = annotated_snv.concat(annotated_indel)
	sigprofiler_ch = SIGPROFILER(all_annotated.map { vcf, type -> tuple(vcf.getParent(), type) })
	snv_sig = sigprofiler_ch.filter { it[1] == 'SNV' }
	indel_sig = sigprofiler_ch.filter { it[1] == 'INDEL' }

	snv_pairs =snv_sig.map { sig_dir, subdir, vcf_dir ->
		tuple(file("${vcf_dir}/*.vcf"),file("${sig_dir}/Assignment_Solution/Activities/Decomposed_Mutation_Probabilities/Decomposed_Mutation_Probabilities_*.txt"), 'SNV')}
	SNV_CONTEXT(snv_pairs)

	indel_pairs =indel_sig.map { sig_dir, subdir, vcf_dir ->
		tuple(file("${vcf_dir}/*.vcf"),file("${sig_dir}/Assignment_Solution/Activities/Decomposed_Mutation_Probabilities/Decomposed_Mutation_Probabilities_*.txt"), 'INDEL')}
	INDEL_CONTEXT(indel_pairs)
}

