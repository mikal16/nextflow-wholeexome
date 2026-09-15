// Annotation step, ported from ngstk's ANNOTATE_VCF stage.
// Uses ngstk's curated vcfanno resource set (assets/vcfanno_grch3{7,8}.toml:
// ClinVar, COSMIC, dbSNP, gnomAD exomes+genomes, CADD, REVEL, dbscSNV,
// SpliceAI, PrimateAI, UCSC segdup/repeat-masker/refseq-select) instead of
// sarek's own annotation resources.
//
// Container tags are pinned to real, published biocontainers images
// (verified against quay.io at the time this pipeline was written); nothing
// here has been run end-to-end, so treat the pins as a starting point and
// re-verify before production use.

process VCFANNO {
    tag "$meta.id"
    label 'process_medium'
    container 'quay.io/biocontainers/vcfanno:0.3.9--h1079eea_0'
    conda 'bioconda::vcfanno=0.3.9'

    input:
    tuple val(meta), path(vcf), path(tbi)
    path toml
    val  resources_dir

    output:
    tuple val(meta), path("${meta.id}.vcfanno.vcf"), emit: vcf

    script:
    """
    RESOURCES=${resources_dir} vcfanno -p ${task.cpus} ${toml} ${vcf} > ${meta.id}.vcfanno.vcf
    """
}

process BGZIP_TABIX {
    tag "$meta.id"
    label 'process_low'
    publishDir "${params.outdir}/annotated_vcf/${meta.id}", mode: 'copy'
    container 'quay.io/biocontainers/htslib:1.24--ha79157c_0'
    conda 'bioconda::htslib=1.24'

    input:
    tuple val(meta), path(vcf)

    output:
    tuple val(meta), path("*.vcf.gz"), path("*.vcf.gz.tbi"), emit: vcf_tbi

    script:
    """
    bgzip -c -@ ${task.cpus} ${vcf} > ${vcf}.gz
    tabix -p vcf ${vcf}.gz
    """
}
