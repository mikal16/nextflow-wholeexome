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
#       --input 'results/sarek/variant_calling/haplotypecaller/*/*.filtered.vcf.gz' \
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

# Adjust to however Nextflow is provisioned on your account, e.g.:
#   module load nextflow
# or activate the same venv ngstk documents in README_Python_Installation.md.

nextflow run main.nf \
    -profile narval_slurm \
    --slurm_account "${SLURM_ACCOUNT}" \
    --resources_dir "${RESOURCES}" \
    -resume \
    "$@"
