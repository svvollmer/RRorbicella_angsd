#!/usr/bin/env bash
#SBATCH --job-name=lfmm_orb
#SBATCH --partition=short
#SBATCH --time=2-00:00:00
#SBATCH --mem=128G
#SBATCH --cpus-per-task=8
#SBATCH --array=0-0
#SBATCH --output=/work/vollmer/orbicella_genomics/logs/lfmm_orb_%a_%j.log

set -euo pipefail

module load R/4.4.0
export R_LIBS_USER=~/R/x86_64-pc-linux-gnu-library/4.4

WORKDIR=/work/vollmer/orbicella_genomics
PROJDIR=/projects/vollmer/coral-angsd-pipeline
PHENO_FILE=/work/vollmer/RR_heat_tolerance_v2/Results/08b_orbicella_blups_LEA.csv

# Runs: orb_combined × blup_scaled (K=3)
# Extend array and these lists to add per-species or additional phenotype runs
SPECIES=(  orb_combined )
PHENOS=(   blup_scaled  )
K_VALS=(   3            )

SP=${SPECIES[$SLURM_ARRAY_TASK_ID]}
PHENO=${PHENOS[$SLURM_ARRAY_TASK_ID]}
K=${K_VALS[$SLURM_ARRAY_TASK_ID]}

OUTDIR="$WORKDIR/results/lfmm/${SP}"
mkdir -p "$WORKDIR/logs" "$OUTDIR"

echo "Running LFMM: species=$SP pheno=$PHENO K=$K"

Rscript "$PROJDIR/scripts/gea/run_lfmm_thermal.R" \
    --species    "$SP" \
    --K          "$K" \
    --K_max      5 \
    --beagle     "$WORKDIR/results/angsd/all.beagle.gz" \
    --bamlist    "$WORKDIR/results/bamlists/all.bamlist" \
    --blup_file  "$PHENO_FILE" \
    --pheno_col  "$PHENO" \
    --exclude    "$WORKDIR/results/lfmm/orb_exclude.txt" \
    --outdir     "$OUTDIR" \
    --min_maf    0.05 \
    --max_miss   0.10 \
    --samp_miss  0.30 \
    --threads    "$SLURM_CPUS_PER_TASK"
