#!/usr/bin/env nextflow

/*
* Variation on VEP pipeline
* AUTHOR: CERC Genomic Medicine,  Daniel Taliun, PhD <daniel.taliun@mcgill.ca>
* VERSION: 3.0
* YEAR: 2022
*/

process annotate {
	label "VEP"

	cache "lenient"
        scratch true
	cpus 1

	containerOptions "-B ${params.vep_cache}:/opt/vep/.vep"

	input:
	tuple path(vcf), path(vcf_index)

	output:
	path("${vcf.getSimpleName()}.vep")
	path "*.vep.log", emit: logs
	publishDir "results", mode: "copy"


	script:
	if (params.assembly == "GRCh38")
		"""
		export PERL5LIB=/opt/vep/.vep/Plugins/:\$PERL5LIB
		loftee_args=human_ancestor_fa:/opt/vep/.vep/loftee_db_${params.assembly}/human_ancestor.fa.gz,gerp_bigwig:/opt/vep/.vep/loftee_db_${params.assembly}/gerp_conservation_scores.homo_sapiens.GRCh38.bw,conservation_file:/opt/vep/.vep/loftee_db_${params.assembly}/loftee.sql${params.loftee_flags}
		bcftools view ${params.drop_genotypes} ${vcf}| vep --cache --offline --assembly ${params.assembly} --format vcf --force_overwrite --dir_cache /opt/vep/.vep/ --plugin LoF,loftee_path:/opt/vep/.vep/loftee_${params.assembly},\${loftee_args} --dir_plugins /opt/vep/.vep/loftee_${params.assembly} --plugin CADD,/opt/vep/.vep/CADD_${params.assembly}/whole_genome_SNVs.tsv.gz,/opt/vep/.vep/CADD_${params.assembly}/InDels.tsv.gz --plugin CONTEXT ${params.vep_flags} --warning_file STDERR --output_file STDOUT > ${vcf.getSimpleName()}.vep 2> ${vcf.getSimpleName()}.vep.log
	 	"""
	else if (params.assembly == "GRCh37")
		"""
		export PERL5LIB=/opt/vep/.vep/Plugins/:\$PERL5LIB
		loftee_args=human_ancestor_fa:/opt/vep/.vep/loftee_db_${params.assembly}/human_ancestor.fa.gz,conservation_file:/opt/vep/.vep/loftee_db_${params.assembly}/phylocsf_gerp.sql${params.loftee_flags}
		bcftools view ${params.drop_genotypes} ${vcf} | vep --cache --offline --assembly ${params.assembly} --format vcf --force_overwrite --dir_cache /opt/vep/.vep/ --plugin LoF,loftee_path:/opt/vep/.vep/loftee_${params.assembly},\${loftee_args} --dir_plugins /opt/vep/.vep/loftee_${params.assembly} --plugin CADD,/opt/vep/.vep/CADD_${params.assembly}/whole_genome_SNVs.tsv.gz,/opt/vep/.vep/CADD_${params.assembly}/InDels.tsv.gz --plugin CONTEXT ${params.vep_flags} --warning_file STDERR --output_file STDOUT > ${vcf.getSimpleName()}.vep 2> ${vcf.getSimpleName()}.vep.log
		"""
	else
		error "Invalid assembly name: ${params.assembly}"
}




workflow {
	vcfs = Channel.fromPath(params.vcfs).map{ vcf -> [ vcf, vcf + ".tbi" ] }
	vcfs_chunks_annotated = annotate(vcfs)
}
