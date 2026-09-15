// Annotation step, ported from ngstk's ANNOTATE_VCF stage.
// Uses ngstk's curated vcfanno resource set (assets/vcfanno_grch3{7,8}.toml:
// ClinVar, COSMIC, dbSNP, gnomAD exomes+genomes, CADD, REVEL, dbscSNV,
// SpliceAI, PrimateAI, UCSC segdup/repeat-masker/refseq-select) instead of
// sarek's own annotation resources.
//
// Container tags are pinned to real, published biocontainers images
// (verified against quay.io) and this whole module has been run end-to-end
// against a synthetic VCF + a full set of dummy files standing in for every
// path in assets/vcfanno_grch3{7,8}.toml (see repo commit history for the
// smoke test). One thing that testing caught: vcfanno does NOT expand
// $RESOURCES itself -- setting the RESOURCES env var and leaving
// "$RESOURCES/..." literal in the toml just makes vcfanno try to open a
// path containing the literal string "$RESOURCES" and fail. The toml has
// to be pre-resolved before vcfanno ever sees it, which is what the sed
// substitution below does.

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
    sed 's#\\\$RESOURCES#${resources_dir}#g' ${toml} > resolved.toml
    vcfanno -p ${task.cpus} resolved.toml ${vcf} > ${meta.id}.vcfanno.vcf
    """
}

process BGZIP_TABIX {
    tag "$meta.id"
    label 'process_low'
    publishDir path: { "${params.outdir}/annotated_vcf/${meta.id}" }, mode: 'copy'
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
