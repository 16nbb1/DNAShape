#!/usr/bin/env nextflow
nextflow.enable.dsl=2

params.input = null
params.output_dir = 'results'
params.snpeff_jar = '/home/nboev/projects/def-sushant/nboev/tools/snpEff/snpEff'
params.snpeff_genome = 'GRCh38.99'
params.snpeff_phast = '/home/nboev/projects/def-sushant/nboev/tools/snpEff/snpEff/db/phastCons'
params.kmer_file = '/home/nboev/projects/def-sushant/nboev/data/refgenome/Jellyfish/11mer_hg38/DNAShape/DNAShape_11.hg38.kmer.tsv'
params.ref_fasta = '/home/nboev/projects/def-sushant/nboev/data/refgenome/hg38.fa'
params.flank_size = '5'

process FILTER_SNV {
	module 'StdEnv/2023'
	module 'vcftools'
	input:
		path vcf
	output:
		path "SNV/*.recode.vcf", emit: snv_vcf
		path "INDEL/*.recode.vcf", emit: indel_vcf
		publishDir "${params.output_dir}", mode: 'copy'
	script:
		"""
		mkdir -p SNV INDEL

		vcftools --vcf ${vcf} --remove-filtered-all --not-chr chrY --remove-indels --recode --recode-INFO-all --out SNV/${vcf.baseName}.SNV
		vcftools --vcf ${vcf} --remove-filtered-all --not-chr chrY --keep-only-indels --recode --recode-INFO-all --out INDEL/${vcf.baseName}.INDEL
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
		java -Xmx8g -jar ${params.snpeff_jar}/snpEff.jar -v -motif ${params.snpeff_genome} ${vcf} > ${vcf.baseName}.annotMotif.vcf
		java -Xmx8g -jar ${params.snpeff_jar}/SnpSift.jar phastCons ${params.snpeff_phast} ${vcf.baseName}.annotMotif.vcf > ${vcf.baseName}.annotMotifPhast.vcf
		rm ${vcf.baseName}.annotMotif.vcf
		rm -f ${vcf}
		"""
}

process SIGPROFILER {
	tag "$subdir"
	input:
		tuple val(annot_dir), val(subdir)
	output:
		path("sigprofiler"), emit: sigprofiler_out
		path("input/*.vcf"),  emit: sigprofiler_vcfs
	publishDir { "${params.output_dir}/${subdir}/sigprofiler" }, mode: 'copy'
	script:
		"""
		module load StdEnv/2023 python scipy-stack
		source ~/virtenv/bin/activate
		mkdir -p sigprofiler input
		mv ${annot_dir}/*.annotMotifPhast.vcf input
		python ${projectDir}/process/sigprofiler.py input sigprofiler ${params.snpeff_genome}
		"""
}

process SNV_CONTEXT {
	input:
		path snv_file
		path decomposed_probs_file
	output:
		path "*.txt", emit: snv_results
	script:
		"""
		module load StdEnv/2023 python scipy-stack
		source ~/ENV/bin/activate
		python ${projectDir}/process/mutcontextSNV.py ${snv_file} ${decomposed_probs_file} ${params.kmer_file} ${params.ref_fasta} ${params.flank_size} >> mutcontextSNV.py.txt
		"""
}

workflow {
	if (!params.input) {
		error "Please provide an input file using --input"
	}

	input_ch = Channel.fromPath(params.input)
	filtered = FILTER_SNV(input_ch)

	snv_vcf = filtered.snv_vcf.map { vcf -> tuple(vcf, 'SNV') }
	indel_vcf = filtered.indel_vcf.map { vcf -> tuple(vcf, 'INDEL') }

	combined_vcfs = snv_vcf.concat(indel_vcf)
	annotated = ANNOT_VCF(combined_vcfs)
	annotated_snv = annotated.filter { it[1] == 'SNV' }
	annotated_indel = annotated.filter { it[1] == 'INDEL' }

	all_annotated = annotated_snv.concat(annotated_indel)
	sigprofiler_out = SIGPROFILER(all_annotated.map { vcf, type -> tuple(vcf.getParent(), type) })

	sig_vcfs = sigprofiler_out.sigprofiler_vcfs
	snv_from_sig = sig_vcfs.filter { it.name.contains('SNV') }

	decomposed_probs = sigprofiler_out.sigprofiler_out.map { folder ->file("${folder}/Assignment_Solution/Activities/Decomposed_Mutation_Probabilities/Decomposed_Mutation_Probabilities_*.txt")}.flatten()
	SNV_CONTEXT(snv_from_sig, decomposed_probs)
}
