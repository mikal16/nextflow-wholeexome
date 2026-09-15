# nextflow-wholeexome

Whole-exome sequencing analysis: [nf-core/sarek](https://nf-co.re/sarek) for
alignment and germline calling, then a small local Nextflow workflow that
re-annotates the VCF and reports it out using
[RI-MUHC-Bioinformatics/ngstk](https://github.com/RI-MUHC-Bioinformatics/ngstk)'s
exome pipeline pieces.

## Why sarek and not nf-co.re/exoseq

The original ask was to build on `nf-co.re/exoseq/dev`. That pipeline is
**archived**: last pushed in 2019, written in DSL1 (won't run on Nextflow
newer than 22.10.6), and nf-core's own repo description now reads *"please
consider using/contributing to nf-core/sarek"*. So this repo builds on sarek
(actively maintained, DSL2, has a WES mode) instead, and layers the ngstk
pieces on top of its output.

## What's here vs. what's vendored from ngstk

This repo does **not** vendor sarek -- it's invoked as a normal upstream
nf-core pipeline (see `docs/usage.md`, Stage 1).

It **does** vendor, essentially unmodified, from
`RI-MUHC-Bioinformatics/ngstk` (`develop` branch):

| ngstk source | Here | Purpose |
|---|---|---|
| `reportVariants.py`, `run_report_variants.py`, `common_tools.py` | `bin/` | Generates ngstk's per-sample TSV variant report (HGNC/HPO-annotated, recessive/denovo detection, etc.) from an annotated VCF. |
| `vcfanno_config_grch37.toml`, `vcfanno_config_grch38.toml` | `assets/` | vcfanno resource set: ClinVar, COSMIC, dbSNP, gnomAD exomes+genomes, CADD, REVEL, dbscSNV, SpliceAI, PrimateAI, UCSC segdup/repeat-masker/refseq. |
| `dnaseq_exome/dnaseq_exome_grch3{7,8}_config*.ini` | `conf/` | Reference paths, filter thresholds, and which report columns to emit. Patched to add 11 report-column keys the vendored `run_report_variants.py` expects that these `.ini` files didn't define -- see `conf/README_ngstk_configs.md`. |
| `slurm_tools/dnaseq_exome_slurm_config.ini` (time/mem/cpu per step) | `conf/base.config` | Per-process SLURM resource requests. |
| `slurm_tools/run_pipeline_via_slurm.py` pattern (sbatch launches a head job that submits one job per pipeline step) | `slurm/submit_wholeexome.sh` + `conf/narval_slurm.config` | Same two-level submission model, using Nextflow's own slurm executor instead of ngstk's bespoke per-step script. |

`main.nf` is new: it's the glue that chains vcfanno annotation, bgzip/tabix,
and the ngstk report step into one workflow, parameterized by `--genome`
(GRCh37/GRCh38) and `--exome_target_set` (ngstk's UCSC-all-exons vs.
UCSC-coding-only GRCh37 variants).

## Quick start

```bash
# Stage 1: sarek (external)
nextflow run nf-core/sarek -r 3.5.1 \
    --input samplesheet.csv --outdir results/sarek \
    --genome GATK.GRCh38 --wes --intervals targets.bed \
    --tools haplotypecaller,snpeff -profile singularity

# Stage 2: this pipeline
nextflow run main.nf \
    --input 'results/sarek/variant_calling/haplotypecaller/*/*.filtered.vcf.gz' \
    --genome GRCh38 \
    --resources_dir /path/to/ngstk/resources \
    --outdir results/wholeexome \
    -profile singularity
```

Full parameter reference, the `$RESOURCES` tree layout, and SLURM/Narval
usage: [`docs/usage.md`](docs/usage.md).

## Status

Scaffolded but **not yet run end-to-end** -- there was no Nextflow/container
runtime or ngstk resource tree available to test against while assembling
this. See "What has not been run" in `docs/usage.md` before trusting it on
real samples.

## Internal use note

`bin/`, `assets/` and the `conf/*.ini` files carry ngstk code and config
from the private `RI-MUHC-Bioinformatics/ngstk` repo. This repo is private
for the same reason -- check with the core before making it public or
sharing it outside the lab.
