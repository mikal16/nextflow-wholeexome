import configparser
import argparse
import common_tools
import reportVariants


def main(vcf_file, config_file, output_directory, variant_caller, annotation_tool, reporting_file):
    """
    Get parameters from config file and call function to generate report
    Args:
        vcf_file: annotated vcf file
        config_file: config file
        output_directory: output directory
        variant_caller: variant caller
        annotation_tool: annotation tool
        reporting_file: the path to the CSV file describing the reports to generate.
    Returns: -
    """
    config = configparser.ConfigParser()
    config.read(config_file)

    # validate output directory
    output_directory = common_tools.validate_directory(output_directory)

    # create the output directory
    common_tools.create_directory(output_directory)

    # values to get from config in section FILTER_VCF
    minimum_depth = config.get("FILTER_VCF", "DEPTH_MIN")
    minimum_depth_alternate = config.get("FILTER_VCF", "DEPTH_ALT_MIN")
    minimum_alternate_frequency = config.get("FILTER_VCF", "ALT_FREQ_MIN")

    # values to get from config in section REPORT_VARIANTS_COLUMNS
    report_parameters = {
        "SAMPLE_ID": config.getboolean("REPORT_VARIANTS_COLUMNS", "SAMPLE_ID"),
        "POSITION": config.getboolean("REPORT_VARIANTS_COLUMNS", "POSITION"),
        "ALLELES": config.getboolean("REPORT_VARIANTS_COLUMNS", "ALLELES"),
        "SAMPLE_GT": config.getboolean("REPORT_VARIANTS_COLUMNS", "SAMPLE_GT"),
        "SAMPLE_GT_QUAL": config.getboolean("REPORT_VARIANTS_COLUMNS", "SAMPLE_GT_QUAL"),
        "FILTER": config.getboolean("REPORT_VARIANTS_COLUMNS", "FILTER"),
        "SAMPLE_DEPTH": config.getboolean("REPORT_VARIANTS_COLUMNS", "SAMPLE_DEPTH"),
        "SAMPLE_DEPTH_REF": config.getboolean("REPORT_VARIANTS_COLUMNS", "SAMPLE_DEPTH_REF"),
        "SAMPLE_DEPTH_ALT": config.getboolean("REPORT_VARIANTS_COLUMNS", "SAMPLE_DEPTH_ALT"),
        "SAMPLE_ALT_AF": config.getboolean("REPORT_VARIANTS_COLUMNS", "SAMPLE_ALT_AF"),
        "VARIANT_GENE_NAME": config.getboolean("REPORT_VARIANTS_COLUMNS", "VARIANT_GENE_NAME"),
        "GENE_LIST": config.getboolean("REPORT_VARIANTS_COLUMNS", "GENE_LIST"),
        "VARIANT_ANNOTATION": config.getboolean("REPORT_VARIANTS_COLUMNS", "VARIANT_ANNOTATION"),
        "VARIANT_HGVS_C": config.getboolean("REPORT_VARIANTS_COLUMNS", "VARIANT_HGVS_C"),
        "VARIANT_HGVS_P": config.getboolean("REPORT_VARIANTS_COLUMNS", "VARIANT_HGVS_P"),
        "VARIANT_EXON_INTRON": config.getboolean("REPORT_VARIANTS_COLUMNS", "VARIANT_EXON_INTRON"),
        "CANONICAL_HGVS_C": config.getboolean("REPORT_VARIANTS_COLUMNS", "CANONICAL_HGVS_C"),
        "CANONICAL_HGVS_P": config.getboolean("REPORT_VARIANTS_COLUMNS", "CANONICAL_HGVS_P"),
        "SELECT_HGVS_C": config.getboolean("REPORT_VARIANTS_COLUMNS", "SELECT_HGVS_C"),
        "SELECT_HGVS_P": config.getboolean("REPORT_VARIANTS_COLUMNS", "SELECT_HGVS_P"),
        "VARIANT_CODING": config.getboolean("REPORT_VARIANTS_COLUMNS", "VARIANT_CODING"),
        "VARIANT_LOF": config.getboolean("REPORT_VARIANTS_COLUMNS", "VARIANT_LOF"),
        "DBSNP_ID": config.getboolean("REPORT_VARIANTS_COLUMNS", "DBSNP_ID"),
        "CLINVAR_ID": config.getboolean("REPORT_VARIANTS_COLUMNS", "CLINVAR_ID"),
        "CLINVAR_SIG": config.getboolean("REPORT_VARIANTS_COLUMNS", "CLINVAR_SIG"),
        "CLINVAR_SIG_CONF": config.getboolean("REPORT_VARIANTS_COLUMNS", "CLINVAR_SIG_CONF"),
        "CLINVAR_DN": config.getboolean("REPORT_VARIANTS_COLUMNS", "CLINVAR_DN"),
        "CLINVAR_REVSTAT": config.getboolean("REPORT_VARIANTS_COLUMNS", "CLINVAR_REVSTAT"),
        "COSMIC_ID": config.getboolean("REPORT_VARIANTS_COLUMNS", "COSMIC_ID"),
        "COSMIC_GENE": config.getboolean("REPORT_VARIANTS_COLUMNS", "COSMIC_GENE"),
        "COSMIC_CDS": config.getboolean("REPORT_VARIANTS_COLUMNS", "COSMIC_CDS"),
        "COSMIC_AA": config.getboolean("REPORT_VARIANTS_COLUMNS", "COSMIC_AA"),
        "COSMIC_CNT": config.getboolean("REPORT_VARIANTS_COLUMNS", "COSMIC_CNT"),
        "GNOMAD_POS": config.getboolean("REPORT_VARIANTS_COLUMNS", "GNOMAD_POS"),
        "GNOMAD_AC": config.getboolean("REPORT_VARIANTS_COLUMNS", "GNOMAD_AC"),
        "GNOMAD_AN": config.getboolean("REPORT_VARIANTS_COLUMNS", "GNOMAD_AN"),
        "GNOMAD_AF": config.getboolean("REPORT_VARIANTS_COLUMNS", "GNOMAD_AF"),
        "GNOMAD_HOM": config.getboolean("REPORT_VARIANTS_COLUMNS", "GNOMAD_HOM"),
        "GNOMAD_FILT": config.getboolean("REPORT_VARIANTS_COLUMNS", "GNOMAD_FILT"),
        "GNOMAD_EX_AC": config.getboolean("REPORT_VARIANTS_COLUMNS", "GNOMAD_EX_AC"),
        "GNOMAD_EX_AN": config.getboolean("REPORT_VARIANTS_COLUMNS", "GNOMAD_EX_AN"),
        "GNOMAD_EX_AF": config.getboolean("REPORT_VARIANTS_COLUMNS", "GNOMAD_EX_AF"),
        "GNOMAD_EX_HOM": config.getboolean("REPORT_VARIANTS_COLUMNS", "GNOMAD_EX_HOM"),
        "GNOMAD_EX_FILT": config.getboolean("REPORT_VARIANTS_COLUMNS", "GNOMAD_EX_FILT"),
        "GNOMAD_WG_AC": config.getboolean("REPORT_VARIANTS_COLUMNS", "GNOMAD_WG_AC"),
        "GNOMAD_WG_AN": config.getboolean("REPORT_VARIANTS_COLUMNS", "GNOMAD_WG_AN"),
        "GNOMAD_WG_AF": config.getboolean("REPORT_VARIANTS_COLUMNS", "GNOMAD_WG_AF"),
        "GNOMAD_WG_HOM": config.getboolean("REPORT_VARIANTS_COLUMNS", "GNOMAD_WG_HOM"),
        "GNOMAD_WG_FILT": config.getboolean("REPORT_VARIANTS_COLUMNS", "GNOMAD_WG_FILT"),
        "GNOMAD_WG_POPMAX": config.getboolean("REPORT_VARIANTS_COLUMNS", "GNOMAD_WG_POPMAX"),
        "GNOMAD_WG_AN_POPMAX": config.getboolean("REPORT_VARIANTS_COLUMNS", "GNOMAD_WG_AN_POPMAX"),
        "GNOMAD_WG_AC_POPMAX": config.getboolean("REPORT_VARIANTS_COLUMNS", "GNOMAD_WG_AC_POPMAX"),
        "GNOMAD_WG_AF_POPMAX": config.getboolean("REPORT_VARIANTS_COLUMNS", "GNOMAD_WG_AF_POPMAX"),
        "GNOMAD_WG_NHOM_POPMAX": config.getboolean("REPORT_VARIANTS_COLUMNS", "GNOMAD_WG_NHOM_POPMAX"),
        "GNOMAD_EX_POPMAX": config.getboolean("REPORT_VARIANTS_COLUMNS", "GNOMAD_EX_POPMAX"),
        "GNOMAD_EX_AN_POPMAX": config.getboolean("REPORT_VARIANTS_COLUMNS", "GNOMAD_EX_AN_POPMAX"),
        "GNOMAD_EX_AC_POPMAX": config.getboolean("REPORT_VARIANTS_COLUMNS", "GNOMAD_EX_AC_POPMAX"),
        "GNOMAD_EX_AF_POPMAX": config.getboolean("REPORT_VARIANTS_COLUMNS", "GNOMAD_EX_AF_POPMAX"),
        "GNOMAD_EX_NHOM_POPMAX": config.getboolean("REPORT_VARIANTS_COLUMNS", "GNOMAD_EX_NHOM_POPMAX"),
        "GNOMAD_EX_NON_CANCER_AF": config.getboolean("REPORT_VARIANTS_COLUMNS", "GNOMAD_EX_NON_CANCER_AF"),
        "SEGDUP": config.getboolean("REPORT_VARIANTS_COLUMNS", "SEGDUP"),
        "REP_MASK": config.getboolean("REPORT_VARIANTS_COLUMNS", "REP_MASK"),
        "REP_SIMPLE": config.getboolean("REPORT_VARIANTS_COLUMNS", "REP_SIMPLE"),
        "CADD_RAW": config.getboolean("REPORT_VARIANTS_COLUMNS", "CADD_RAW"),
        "CADD_PHRED": config.getboolean("REPORT_VARIANTS_COLUMNS", "CADD_PHRED"),
        "CADD_PHRED_SNV": config.getboolean("REPORT_VARIANTS_COLUMNS", "CADD_PHRED_SNV"),
        "CADD_PHRED_INDEL": config.getboolean("REPORT_VARIANTS_COLUMNS", "CADD_PHRED_INDEL"),
        "SPLICEAI_GENE": config.getboolean("REPORT_VARIANTS_COLUMNS", "SPLICEAI_GENE"),
        "SPLICEAI_DS_MAX": config.getboolean("REPORT_VARIANTS_COLUMNS", "SPLICEAI_DS_MAX"),
        "SPLICEAI_TYPE": config.getboolean("REPORT_VARIANTS_COLUMNS", "SPLICEAI_TYPE"),
        "SPLICEAI_POS": config.getboolean("REPORT_VARIANTS_COLUMNS", "SPLICEAI_POS"),
        "PRIMATEAI": config.getboolean("REPORT_VARIANTS_COLUMNS", "PRIMATEAI"),
        "REMM": config.getboolean("REPORT_VARIANTS_COLUMNS", "REMM"),
        "REVEL": config.getboolean("REPORT_VARIANTS_COLUMNS", "REVEL"),
        "BAYESDEL_NO_AF_SCORE": config.getboolean("REPORT_VARIANTS_COLUMNS", "BAYESDEL_NO_AF_SCORE"),
        "BAYESDEL_ADD_AF_SCORE": config.getboolean("REPORT_VARIANTS_COLUMNS", "BAYESDEL_ADD_AF_SCORE"),
        "METASVM_SCORE": config.getboolean("REPORT_VARIANTS_COLUMNS", "METASVM_SCORE"),
        "METALR_SCORE": config.getboolean("REPORT_VARIANTS_COLUMNS", "METALR_SCORE"),
        "METARNN_SCORE": config.getboolean("REPORT_VARIANTS_COLUMNS", "METARNN_SCORE"),
        "DBSCSNV_ADA": config.getboolean("REPORT_VARIANTS_COLUMNS", "DBSCSNV_ADA"),
        "DBSCSNV_RF": config.getboolean("REPORT_VARIANTS_COLUMNS", "DBSCSNV_RF"),
        "GNOMAD_SYN_Z": config.getboolean("REPORT_VARIANTS_COLUMNS", "GNOMAD_SYN_Z"),
        "GNOMAD_MIS_Z": config.getboolean("REPORT_VARIANTS_COLUMNS", "GNOMAD_MIS_Z"),
        "GNOMAD_LOF_Z": config.getboolean("REPORT_VARIANTS_COLUMNS", "GNOMAD_LOF_Z"),
        "HPO_DISEASE_OLD": config.getboolean("REPORT_VARIANTS_COLUMNS", "HPO_DISEASE"),
        "HPO_DISEASE": config.getboolean("REPORT_VARIANTS_COLUMNS", "HPO_DISEASE"),
        "VARIANT_NUM_CALLED": config.getboolean("REPORT_VARIANTS_COLUMNS", "VARIANT_NUM_CALLED"),
        "VARIANT_NUM_HOM_REF": config.getboolean("REPORT_VARIANTS_COLUMNS", "VARIANT_NUM_HOM_REF"),
        "VARIANT_NUM_HET": config.getboolean("REPORT_VARIANTS_COLUMNS", "VARIANT_NUM_HET"),
        "VARIANT_NUM_HOM_ALT": config.getboolean("REPORT_VARIANTS_COLUMNS", "VARIANT_NUM_HOM_ALT"),
        "VARIANT_NUM_UNKNOWN": config.getboolean("REPORT_VARIANTS_COLUMNS", "VARIANT_NUM_UNKNOWN"),
        "DENOVO": config.getboolean("REPORT_VARIANTS_COLUMNS", "DENOVO"),
        "RECESSIVE": config.getboolean("REPORT_VARIANTS_COLUMNS", "RECESSIVE"),
        "PHASING": config.getboolean("REPORT_VARIANTS_COLUMNS", "PHASING"),
        "GENE_COUNT": config.getboolean("REPORT_VARIANTS_COLUMNS", "GENE_COUNT"),
        "COMMENTS": config.getboolean("REPORT_VARIANTS_COLUMNS", "COMMENTS"),
        "AM_SCORE": config.getboolean("REPORT_VARIANTS_COLUMNS", "AM_SCORE"),
        "AM_CLASS": config.getboolean("REPORT_VARIANTS_COLUMNS", "AM_CLASS"),
    }

    command_generator = reportVariants.ReportVariants(config_file, vcf_file, variant_caller, output_directory, annotation_tool, reporting_file, report_parameters)

    command_generator.generate_report(minimum_depth, minimum_depth_alternate, minimum_alternate_frequency)


if __name__ == "__main__":
    # enable command line arguments parsing. Easy to use and to add arguments. See : https://docs.python.org/fr/2/howto/argparse.html
    parser = argparse.ArgumentParser()
    parser.add_argument("-vcf", help="File path to a vcf file", required=True)
    parser.add_argument("-cf", help="Config file", required=True)
    parser.add_argument("-o", help="Output directory where all files will be created", required=True)
    parser.add_argument("-vc", help="variant caller", required=True)
    parser.add_argument("-at", help="annotation tool", required=True)
    parser.add_argument("-rf", help="CSV file to organize the report generation. "
                                    "Comma separated list of sample names to include in a report with no header. "
                                    "One line per report. If no csv is provided, all samples will be in one report.", default="")
    args = parser.parse_args()
    main(args.vcf, args.cf, args.o, args.vc, args.at, args.rf)
