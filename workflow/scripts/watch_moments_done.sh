#!/usr/bin/env bash
# =============================================================================
# watch_moments_done.sh — Poll for apal_vs_acer_fl moments completion,
# then auto-run summarize steps and launch SMC++ plot.
#
# Run from login node:
#   nohup bash /projects/vollmer/coral-angsd-pipeline/workflow/scripts/watch_moments_done.sh \
#       > /work/vollmer/acropora_genomics/logs/watcher_moments.log 2>&1 &
# =============================================================================

set -euo pipefail

JOBID=51651910
WORK=/work/vollmer/acropora_genomics
PIPELINE=/projects/vollmer/coral-angsd-pipeline
PYTHON=/home/s.vollmer/.conda/envs/snakemake2/bin/python
SINGULARITY_IMG=/work/vollmer/acropora_genomics/.snakemake/singularity/084f8dcfcfe42901022fa6e05278b43b.simg

INFERENCE_DIR=$WORK/results/demography/inference
OUTPUT_JSON=$INFERENCE_DIR/apal_vs_acer_fl.moments.json
SMCPP_APAL=$WORK/results/smcpp/fit/apal/model.final.json
SMCPP_ACER=$WORK/results/smcpp/fit/acer/model.final.json
SMCPP_PLOT_PNG=$WORK/results/smcpp/plots/smcpp_Ne.png
SMCPP_APAL_CSV=$WORK/results/smcpp/plots/apal_Ne.csv
SMCPP_ACER_CSV=$WORK/results/smcpp/plots/acer_Ne.csv

POLL_INTERVAL=300  # seconds between checks

echo "========================================================"
echo "  watch_moments_done.sh started: $(date)"
echo "  Watching SLURM job: $JOBID"
echo "  Waiting for: $OUTPUT_JSON"
echo "========================================================"

# ── Poll loop ──────────────────────────────────────────────
while true; do
    # Check if output file exists (job finished successfully)
    if [[ -f "$OUTPUT_JSON" ]]; then
        echo "[$(date)] OUTPUT DETECTED: $OUTPUT_JSON"
        break
    fi

    # Check if job is still in queue/running
    JOB_STATE=$(squeue -j "$JOBID" -h -o '%T' 2>/dev/null || echo "GONE")
    if [[ "$JOB_STATE" == "GONE" || -z "$JOB_STATE" ]]; then
        echo "[$(date)] Job $JOBID no longer in queue."
        if [[ ! -f "$OUTPUT_JSON" ]]; then
            echo "[$(date)] ERROR: Job gone but output not found — job may have failed or hit wall limit."
            echo "[$(date)] Check: sacct -j $JOBID --format=JobID,State,ExitCode,Elapsed"
            sacct -j "$JOBID" --format=JobID,State,ExitCode,Elapsed 2>/dev/null || true
            exit 1
        fi
        break
    fi

    echo "[$(date)] Job $JOBID state=$JOB_STATE — waiting..."
    sleep "$POLL_INTERVAL"
done

# ── Step 1: Print best-fit results ────────────────────────
echo ""
echo "[$(date)] ── Step 1: Best-fit results ──"
BESTFIT=$INFERENCE_DIR/apal_vs_acer_fl.bestfit.txt
if [[ -f "$BESTFIT" ]]; then
    echo "Best-fit summary:"
    cat "$BESTFIT"
else
    echo "No bestfit.txt found — reading JSON directly:"
    $PYTHON -c "
import json, sys
with open('$OUTPUT_JSON') as f:
    d = json.load(f)
best = min(d['models'], key=lambda m: -m['log_likelihood'])
print(f\"Best model: {best['model']}\")
print(f\"Log-likelihood: {best['log_likelihood']:.2f}\")
for k,v in best['params'].items():
    print(f\"  {k}: {v:.6g}\")
" 2>/dev/null || $PYTHON -c "import json; d=json.load(open('$OUTPUT_JSON')); print(json.dumps(d, indent=2))"
fi

# ── Step 2: Run summarize across all 7 comparisons ────────
echo ""
echo "[$(date)] ── Step 2: Summarize all moments comparisons ──"
SUMMARIZE=$PIPELINE/workflow/scripts/summarize_moments.py
if [[ -f "$SUMMARIZE" ]]; then
    $PYTHON "$SUMMARIZE" \
        --inference-dir "$INFERENCE_DIR" \
        --out "$WORK/results/demography/moments_summary.tsv" \
        2>&1 && echo "[$(date)] Summary written: $WORK/results/demography/moments_summary.tsv" \
              || echo "[$(date)] WARNING: summarize_moments.py failed — check manually"
else
    echo "[$(date)] WARNING: summarize_moments.py not found at $SUMMARIZE — skipping"
fi

# ── Step 3: SMC++ plot (both species) ─────────────────────
echo ""
echo "[$(date)] ── Step 3: SMC++ Ne(t) plot ──"
if [[ -f "$SMCPP_APAL" && -f "$SMCPP_ACER" ]]; then
    mkdir -p "$WORK/results/smcpp/plots"

    singularity exec "$SINGULARITY_IMG" \
        smc++ plot \
        -g 5 --csv \
        "$SMCPP_PLOT_PNG" \
        "$SMCPP_APAL" "$SMCPP_ACER" \
        2>&1 && echo "[$(date)] SMC++ plot written: $SMCPP_PLOT_PNG" \
              || echo "[$(date)] WARNING: smc++ plot failed"

    # Copy CSVs to docs/figures for plotting script
    FIGS=$WORK/docs/figures
    mkdir -p "$FIGS"
    [[ -f "$SMCPP_APAL_CSV" ]] && cp "$SMCPP_APAL_CSV" "$FIGS/apal_Ne.csv" && echo "[$(date)] Copied apal_Ne.csv"
    [[ -f "$SMCPP_ACER_CSV" ]] && cp "$SMCPP_ACER_CSV" "$FIGS/acer_Ne.csv" && echo "[$(date)] Copied acer_Ne.csv"
else
    echo "[$(date)] WARNING: SMC++ model files not found — skipping plot"
    [[ ! -f "$SMCPP_APAL" ]] && echo "  Missing: $SMCPP_APAL"
    [[ ! -f "$SMCPP_ACER" ]] && echo "  Missing: $SMCPP_ACER"
fi

# ── Step 4: Regenerate annotated Ne plot ──────────────────
echo ""
echo "[$(date)] ── Step 4: Regenerate annotated SMC++ figure ──"
PLOT_SCRIPT=$PIPELINE/workflow/scripts/plot_smcpp.py
if [[ -f "$PLOT_SCRIPT" && -f "$FIGS/apal_Ne.csv" && -f "$FIGS/acer_Ne.csv" ]]; then
    $PYTHON "$PLOT_SCRIPT" \
        --apal "$FIGS/apal_Ne.csv" \
        --acer "$FIGS/acer_Ne.csv" \
        --out  "$FIGS/smcpp_Ne_annotated.png" \
        2>&1 && echo "[$(date)] Annotated figure written: $FIGS/smcpp_Ne_annotated.png" \
              || echo "[$(date)] WARNING: plot_smcpp.py failed — check manually"
else
    echo "[$(date)] Skipping annotated plot (missing script or CSVs)"
fi

# ── Done ──────────────────────────────────────────────────
echo ""
echo "========================================================"
echo "  watch_moments_done.sh COMPLETE: $(date)"
echo ""
echo "  Next steps:"
echo "    1. Review moments_summary.tsv and update RESULTS.md"
echo "    2. git add/commit figures and results"
echo "    3. Run: bash run.sh 7  (if SMC++ plot needs Snakemake)"
echo "========================================================"
