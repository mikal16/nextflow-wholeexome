# About these config files

`dnaseq_exome_grch37_config.ini`, `dnaseq_exome_grch38_config.ini`,
`dnaseq_exome_grch37_config_ucsc_all_exons.ini` and
`dnaseq_exome_grch37_config_ucsc_coding.ini` are copied from
[RI-MUHC-Bioinformatics/ngstk](https://github.com/RI-MUHC-Bioinformatics/ngstk)
(`dnaseq_exome/`, `develop` branch) essentially unchanged, with one fix:

`bin/run_report_variants.py` (also vendored from ngstk, `develop` branch)
reads more `[REPORT_VARIANTS_COLUMNS]` keys than these four `.ini` files
defined -- `VARIANT_EXON_INTRON`, `CADD_PHRED_SNV`, `CADD_PHRED_INDEL`,
`BAYESDEL_NO_AF_SCORE`, `BAYESDEL_ADD_AF_SCORE`, `METASVM_SCORE`,
`METALR_SCORE`, `METARNN_SCORE`, `GNOMAD_EX_NON_CANCER_AF`, `AM_SCORE`,
`AM_CLASS` -- so as vendored they would fail with a `configparser`
`NoOptionError` before producing a report. This means the exome-module
configs in ngstk's own repo are currently out of sync with its
`run_report_variants.py`; worth flagging back upstream.

All eleven were appended set to `False` (i.e. omit that column from the
report), since none of the corresponding annotations are present in
`assets/vcfanno_grch37.toml` / `vcfanno_grch38.toml` as vendored here.
If you add the underlying resources (e.g. AlphaMissense for
`AM_SCORE`/`AM_CLASS`, split SNV/indel CADD, dbNSFP for the BayesDel/MetaSVM/
MetaLR/MetaRNN meta-scores) to the vcfanno TOML, flip the matching key to
`True`.

Everything else in these files -- reference genome paths, GATK/VCF filtering
thresholds, HPO/HGNC/gnomAD-metrics resource paths -- is untouched from
ngstk.
