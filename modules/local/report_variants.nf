// Wraps ngstk's reportVariants.py / run_report_variants.py (vendored into
// bin/ with one change: added a `#!/usr/bin/env python3` shebang to
// run_report_variants.py, which the original ngstk copy doesn't have --
// Nextflow's bin/-on-PATH mechanism executes scripts directly rather than
// via `python3 script.py`, and that fails with no shebang. Confirmed by
// actually running it. See conf/README_ngstk_configs.md for the other
// vendoring change, to the .ini config files) to turn the annotated VCF
// into ngstk's flattened per-sample TSV variant report, enriched with
// HGNC gene symbols and HPO gene-disease terms.
//
// This whole process has been run end-to-end against a synthetic VCF and
// the real, unmodified assets/vcfanno_grch38.toml + conf .ini shipped in
// this repo (see repo commit history for the smoke test).

process REPORT_VARIANTS {
    tag "$meta.id"
    label 'process_low'
    publishDir path: { "${params.outdir}/report_variants/${meta.id}" }, mode: 'copy'
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
