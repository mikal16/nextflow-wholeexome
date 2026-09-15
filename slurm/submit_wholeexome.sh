#!/bin/bash
# Submits the Nextflow head job for this pipeline to SLURM on Narval, in the
# same two-level pattern ngstk's slurm_tools/run_pipeline_via_slurm.py uses:
# one small sbatch job runs Nextflow itself, which then submits one SLURM
# job per pipeline task via the narval_slurm profile (conf/narval_slurm.config).
#
# Usage:
#   RESOURCES=/lustre06/project/rrg-jbriv/resources \
#   SLURM_ACCOUNT=rrg-jbriv \
#   ./slurm/submit_wholeexome.sh \
#       --input 'results/sarek/annotation/*/*_snpEff.ann.vcf.gz' \
#       --genome GRCh38 \
#       --outdir results/wholeexome
#
# Any extra arguments are passed straight through to `nextflow run main.nf`.

#SBATCH --job-name=nextflow-wholeexome
#SBATCH --time=24:00:00
#SBATCH --mem=4G
#SBATCH --cpus-per-task=1
#SBATCH --output=logs/nextflow-wholeexome_%j.log

set -euo pipefail

: "${RESOURCES:?Set RESOURCES to the root of the ngstk resource tree}"
: "${SLURM_ACCOUNT:?Set SLURM_ACCOUNT to your Narval RAP ID (e.g. rrg-jbriv)}"

mkdir -p logs

# Everything this pipeline shells out to -- nextflow, vcfanno, bgzip/tabix,
# and the venv for bin/run_report_variants.py -- needs to be on $PATH here,
# before `nextflow run`, because SLURM propagates this job's environment to
# the per-task jobs Nextflow submits underneath it. Adjust these to however
# you provisioned each one on your account (see docs/usage.md, "Setting up
# Narval" for one concrete way to do it):

# module load StdEnv/2023 java/17 nextflow  # if provided as modules; else self-install
# module load StdEnv/2023 htslib            # bgzip / tabix
# export PATH="/path/to/vcfanno_dir:$PATH"  # static vcfanno binary (no module on Narval)
# source /path/to/venv/bin/activate         # pandas/numpy/cyvcf2/pyfiglet/sample-sheet/xlsxwriter

nextflow run main.nf \
    -profile narval_slurm \
    --slurm_account "${SLURM_ACCOUNT}" \
    --resources_dir "${RESOURCES}" \
    -resume \
    "$@"
