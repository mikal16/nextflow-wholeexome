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

## What has actually been tested

The full chain (`VCFANNO` -> `BGZIP_TABIX` -> `REPORT_VARIANTS`) has been
run end-to-end with real Nextflow (26.04.6) and a real vcfanno binary
(v0.3.9, matching the pinned container tag) via `-profile standard`,
against:

- a synthetic single-sample, single-variant GATK/SnpEff-style VCF,
- a full set of dummy files standing in for every one of the 63 paths
  referenced in the real, unmodified `assets/vcfanno_grch38.toml`,
  correctly bgzipped and tabix-indexed,
- the real, unmodified `conf/dnaseq_exome_grch38_config.ini`, plus dummy
  gnomAD-metrics/canonical-RefSeq/HGNC/HPO resource files.

Confirmed working: annotations written by vcfanno (ClinVar/COSMIC/dbSNP
IDs, CADD, REVEL, dbscSNV, segdup/repeat flags, UCSC RefSeq-select) flow
through `bgzip`/`tabix` into `bin/reportVariants.py` and come out correctly
in the final `report_variants.tsv`, alongside gnomAD z-scores and HPO
disease terms pulled from the report-side resource files. All four
`--genome`/`--exome_target_set` combinations were also confirmed to select
the correct `.ini`/`.toml` pair (`GRCh37`: `default`, `ucsc_all_exons`,
`ucsc_coding`; `GRCh38`: `default`).

This testing found and fixed four real bugs that a syntax read-through
hadn't caught:

1. **`nextflow.config` used a top-level `def` statement** to pick
   `conf/grch3{7,8}.config` -- current Nextflow's stricter parser rejects
   bare statements outside `params`/`process`/etc. blocks. Fixed by folding
   it into the ternary directly.
2. **`bin/run_report_variants.py` had no shebang.** Nextflow's `bin/`
   auto-`$PATH` mechanism executes scripts directly by their executable
   bit; without `#!/usr/bin/env python3`, bash tried to run it as a shell
   script and failed on the first `import` line. One line added; see the
   comment in `modules/local/report_variants.nf`.
3. **vcfanno does not expand `$RESOURCES` itself.** The original design
   set the `RESOURCES` env var and left `"$RESOURCES/..."` literal in the
   TOML, assuming vcfanno would do shell-style expansion the way ngstk's
   own Python config parsing does. It doesn't -- it tries to open a path
   containing the literal four characters `$RES...` and fails. Fixed by
   `sed`-substituting the real path into the TOML before invoking vcfanno
   (`modules/local/vcfanno.nf`).
4. **`publishDir` referencing `meta.id`** (a per-task input variable) as a
   plain interpolated string is evaluated once at process-definition time,
   before `meta` exists, and throws `No such variable: meta`. Fixed by
   switching to the dynamic-directive closure form,
   `publishDir path: { "...${meta.id}" }, mode: 'copy'`.

`common_tools.py` (imported by `reportVariants.py`) also turned out to
import `pyfiglet`, `sample_sheet` and `xlsxwriter` at module level, even
though the report step never calls them -- they're only used by ngstk's
other scripts that share this file. `modules/local/environment.yml` didn't
list them; it does now.

Not yet exercised by this testing, so still worth checking before trusting
this on real samples:

- Real (non-dummy) reference data of realistic size/format -- e.g. multi-
  contig gnomAD VCFs with proper `##INFO` header declarations. The dummy
  resource files here mostly lacked those, so vcfanno logged (non-fatal)
  "Info Error: ... not found in header" warnings that a real resource file
  wouldn't produce.
- Multi-sample / trio VCFs and the `--reporting_file` trio/denovo/recessive
  logic (`DENOVO`, `RECESSIVE`, `PHASING` columns) -- the test VCF had one
  sample and no family structure.
- Actually pulling the pinned containers (`vcfanno:0.3.9--h1079eea_0`,
  `htslib:1.24--ha79157c_0`) under Docker/Singularity -- this testing ran
  the same binary versions directly on `$PATH`, not through the containers.
- The GRCh37 `.ini`/`.toml` pair's actual annotation output -- confirmed
  the pipeline selects the right files for all three GRCh37
  `--exome_target_set` values, but didn't build a second 63-file dummy
  resource tree to run annotation through them (would be identical
  mechanics to the GRCh38 run already verified).
- `slurm/submit_wholeexome.sh` / `-profile narval_slurm` -- no SLURM
  scheduler available to test against.
