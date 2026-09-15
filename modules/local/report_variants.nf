// Wraps ngstk's reportVariants.py / run_report_variants.py (vendored
// unmodified into bin/, aside from the config .ini files it reads from --
// see conf/README_ngstk_configs.md for what was patched and why) to turn
// the annotated VCF into ngstk's flattened per-sample TSV variant report,
// enriched with HGNC gene symbols and HPO gene-disease terms.
//
// Nextflow automatically puts this pipeline's bin/ on $PATH, so
// run_report_variants.py, reportVariants.py and common_tools.py are found
// without any install step -- exactly as they'd import in the original
// ngstk checkout.

process REPORT_VARIANTS {
    tag "$meta.id"
    label 'process_low'
    publishDir "${params.outdir}/report_variants/${meta.id}", mode: 'copy'
    conda "${moduleDir}/environment.yml"

    input:
    tuple val(meta), path(vcf), path(tbi)
    path config_ini
    path reporting_file
    val  resources_dir

    output:
    path "report/*.tsv", emit: reports

    script:
    def rf_arg = reporting_file.name != 'NO_FILE' ? "-rf ${reporting_file}" : ''
    """
    export RESOURCES=${resources_dir}
    mkdir -p report
    run_report_variants.py \\
        -vcf ${vcf} \\
        -cf ${config_ini} \\
        -o report/ \\
        -vc ${params.variant_caller} \\
        -at ${params.annotation_tool} \\
        ${rf_arg}
    """
}
