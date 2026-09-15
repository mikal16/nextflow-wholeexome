#!/usr/bin/env nextflow
/*
 * nextflow-wholeexome
 *
 * Stage 2 of a two-stage whole-exome workflow:
 *   1. nf-core/sarek (run separately, see README.md) does trimming,
 *      alignment, duplicate marking, BQSR and germline variant calling,
 *      producing a SnpEff-annotated (or annotate here) VCF per sample.
 *   2. THIS pipeline re-annotates that VCF with the ngstk exome resource
 *      set (assets/vcfanno_*.toml) and generates ngstk's per-sample TSV
 *      variant report (bin/reportVariants.py), matching the reporting
 *      ngstk's own dnaseq_exome module has always produced.
 *
 * nf-core/exoseq was deliberately NOT used as the base: it has been
 * archived since 2019, is DSL1, and nf-core's own repo now points users to
 * nf-core/sarek instead.
 */

nextflow.enable.dsl = 2

include { VCFANNO; BGZIP_TABIX } from './modules/local/vcfanno'
include { REPORT_VARIANTS }      from './modules/local/report_variants'

def exomeConfigs() {
    [
        'GRCh38': [
            'default'        : params.grch38_config,
        ],
        'GRCh37': [
            'default'        : params.grch37_config,
            'ucsc_all_exons' : params.grch37_config_ucsc_all_exons,
            'ucsc_coding'    : params.grch37_config_ucsc_coding,
        ],
    ]
}

workflow {

    if (!params.input) {
        error "Please provide --input, a glob to sarek's germline VCF(s), " +
              "e.g. --input 'results/sarek/variant_calling/haplotypecaller/*/*.filtered.vcf.gz'"
    }
    if (!(params.genome in ['GRCh37', 'GRCh38'])) {
        error "--genome must be 'GRCh37' or 'GRCh38' (got '${params.genome}')"
    }
    if (!params.resources_dir) {
        error "Please provide --resources_dir, the root of the ngstk resource " +
              "tree (gnomad/clinvar/cosmic/dbsnp/cadd/revel/dbscsnv/segdup/hgnc/hpo). " +
              "See docs/usage.md."
    }

    def configs = exomeConfigs()
    def configPath = configs[params.genome][params.exome_target_set]
    if (!configPath) {
        error "Unknown --exome_target_set '${params.exome_target_set}' for genome " +
              "${params.genome}. Valid options: ${configs[params.genome].keySet().join(', ')}"
    }
    config_ini = file(configPath, checkIfExists: true)

    toml = file(
        params.genome == 'GRCh37' ? params.vcfanno_toml_grch37 : params.vcfanno_toml_grch38,
        checkIfExists: true
    )

    reporting_file = params.reporting_file ? file(params.reporting_file, checkIfExists: true) : file('NO_FILE')

    vcf_ch = Channel.fromPath(params.input, checkIfExists: true)
        .map { vcf ->
            def tbiCandidate = file("${vcf}.tbi")
            if (!tbiCandidate.exists()) {
                error "Missing tabix index for ${vcf} (expected ${tbiCandidate}). " +
                      "Run 'tabix -p vcf ${vcf}' first."
            }
            def meta = [ id: vcf.simpleName ]
            tuple(meta, vcf, tbiCandidate)
        }

    VCFANNO(vcf_ch, toml, params.resources_dir)
    BGZIP_TABIX(VCFANNO.out.vcf)
    REPORT_VARIANTS(BGZIP_TABIX.out.vcf_tbi, config_ini, reporting_file, params.resources_dir)

    // A custom workflow.onComplete{} handler was tried here and dropped: on
    // Nextflow's newer strict-syntax parser (confirmed against 26.04.6) it
    // can no longer live at the top level ("Statements cannot be mixed with
    // script declarations"), and nesting it inside this workflow block
    // instead makes `workflow` resolve to null within the closure
    // (NullPointerException on workflow.complete), confirmed by actually
    // running it. Nextflow's own default completion summary covers this.
}
