# Usage

## Overview

This is a two-stage workflow:

1. **Alignment + calling** -- run [nf-core/sarek](https://nf-co.re/sarek) (not
   vendored here; invoked as a normal upstream nf-core pipeline) to go from
   FASTQ to a germline VCF, using its WES/targeted-sequencing mode.
2. **Annotation + reporting** -- this repo (`main.nf`) takes that VCF,
   re-annotates it with ngstk's exome resource set via vcfanno, and produces
   ngstk's per-sample TSV variant report.

`nf-core/exoseq` was not used: it's archived (last release 2019, DSL1, broken
on Nextflow >22.10.6) and nf-core's own repo description now points users to
sarek instead.

## Stage 1: nf-core/sarek

```bash
nextflow run nf-core/sarek -r 3.5.1 \
    --input samplesheet.csv \
    --outdir results/sarek \
    --genome GATK.GRCh38 \
    --wes \
    --intervals /path/to/exome_targets.bed \
    --tools haplotypecaller,snpeff \
    -profile singularity
```

Points that matter for stage 2 compatibility:

- **`--tools` must include `snpeff`.** `bin/reportVariants.py` (vendored
  from ngstk) currently only implements the SnpEff `ANN` INFO-field parser
  (`annotation_tool=snpeff`); a VEP-annotated VCF will not report correctly.
- **`--wes` + `--intervals`** restricts calling to your capture/exome
  regions; use the same target BED you configure in stage 2 (see below) to
  keep the two stages consistent.
- Use `--genome GATK.GRCh37` for GRCh37 samples.

See the [sarek docs](https://nf-co.re/sarek/latest/docs/usage) for full
parameter details (this pipeline does not wrap or re-document sarek itself).

## Stage 2: this pipeline

```bash
nextflow run main.nf \
    --input 'results/sarek/variant_calling/haplotypecaller/*/*.filtered.vcf.gz' \
    --genome GRCh38 \
    --resources_dir /path/to/ngstk/resources \
    --reporting_file reporting_file.csv \
    --outdir results/wholeexome \
    -profile singularity
```

| Parameter | Required | Notes |
|---|---|---|
| `--input` | yes | Glob to sarek's germline VCF(s); each must have a `.tbi` alongside it. |
| `--genome` | yes | `GRCh37` or `GRCh38`. |
| `--resources_dir` | yes | Root of the ngstk resource tree (see below). Exported as `$RESOURCES`; both the vendored `.ini` config and the vcfanno `.toml` reference it. |
| `--exome_target_set` | no | `default` (both builds), or for GRCh37 also `ucsc_all_exons` / `ucsc_coding` -- mirrors ngstk's three GRCh37 exome-region variants. |
| `--reporting_file` | no | ngstk-style CSV grouping sample IDs into one report per line (trio detection etc. -- see ngstk's `dnaseq_exome/README.md`). Omit to get one combined report of all samples in the VCF. |
| `--variant_caller` | no | Default `GATK`, matching sarek's HaplotypeCaller output. |
| `--outdir` | no | Default `results`. |

### `$RESOURCES` layout

The vendored configs and vcfanno TOMLs expect, under `--resources_dir`:

```
grch38/  (or grch37/)
├── broad/            # reference fasta + GATK bundle (not used in stage 2 directly)
├── clinvar/clinvar.vcf.gz(.tbi)
├── cosmic/cosmic_coding.vcf.gz, cosmic_noncoding.vcf.gz
├── dbsnp/dbsnp.vcf.gz(.tbi)
├── cadd/cadd.snvs.tsv.gz, cadd.indels.tsv.gz
├── revel/revel.tsv.gz
├── dbscsnv/dbscsnv.tsv.gz
├── gnomad/
│   ├── gnomad.exomes.chr{1..22,X,Y}.vcf.gz
│   ├── gnomad.genomes.chr{1..22,X,Y}.vcf.gz
│   └── gnomad.lof_metrics.by_gene.txt
└── ucsc/
    ├── genomic_superdups.bed.gz, repeat_masker.bed.gz, simple_repeats.bed.gz
    ├── refseq/refseq_select.bed.gz, refseq_curated_all_exons.{bed,interval_list}
    └── canonical.refseq.tsv

hgnc/hgnc_complete_set.txt
hpo/genes_for_HP_0000001.csv, diseases_inheritence.csv, mode_inheritance/diseases_for_HP_0000001.tsv
hpo_new/hpo_genes_diseases.tsv
```

This mirrors ngstk's own resource layout (`RI-MUHC-Bioinformatics/ngstk`,
`dnaseq_exome/README.md`) -- if your lab already has this tree provisioned
for ngstk, point `--resources_dir` at the same root.

## Running on Narval via SLURM

```bash
RESOURCES=/lustre06/project/rrg-jbriv/resources \
SLURM_ACCOUNT=rrg-jbriv \
./slurm/submit_wholeexome.sh \
    --input 'results/sarek/variant_calling/haplotypecaller/*/*.filtered.vcf.gz' \
    --genome GRCh38 \
    --outdir results/wholeexome
```

This submits one sbatch job to run the Nextflow head process, which then
submits one SLURM job per task via `-profile narval_slurm`
(`conf/narval_slurm.config`) -- the same two-level submission model as
ngstk's `slurm_tools/run_pipeline_via_slurm.py`, just delegated to
Nextflow's built-in slurm executor instead of ngstk's bespoke script.

## What has *not* been run

Nothing in this pipeline has been executed end-to-end -- there's no
Nextflow/container runtime or ngstk resource tree available in the
environment this was assembled in. Before trusting it on real samples:

- Run it once on a small test VCF and confirm `vcfanno` actually resolves
  every `$RESOURCES/...` path in `assets/vcfanno_grch3{7,8}.toml`.
- Confirm the container tags in `modules/local/vcfanno.nf`
  (`vcfanno:0.3.9--h1079eea_0`, `htslib:1.24--ha79157c_0`) pull cleanly in
  your environment.
- Confirm `bin/reportVariants.py` produces output you'd expect against a
  known VCF -- it's vendored unmodified from ngstk `develop`, so it should
  behave exactly as it does there, but that hasn't been re-verified here.
