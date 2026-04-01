#!/usr/bin/env bash
# =============================================================================
# run.sh — Launch RRorbicella_angsd pipeline
#
# Usage:
#   bash run.sh <segment> [--dry-run]
#
# Segments:
#   1   — BAM ingestion (local_bam symlink + merge_bam merge) + QC gate
#   2a  — SNP discovery + GL + relatedness (stops at clone gate)
#   2b  — Subset BEAGLE after clone_approve.py (completes Segment 2)
#   3   — LD pruning, PCA, admixture
#   4   — SAF, SFS, theta, FST
# =============================================================================

set -euo pipefail

PIPELINE=/projects/vollmer/RRorbicella_angsd
SNAKEMAKE=/home/s.vollmer/.conda/envs/snakemake2/bin/snakemake
PROJECT="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROFILE="$PIPELINE/profiles/discovery"

DRY_RUN=""
SEGMENT=""
EXTRA_ARGS=()

while [[ $# -gt 0 ]]; do
    case "$1" in
        1|2a|2b|2c|3|4) SEGMENT="$1" ;;
        --dry-run|-n) DRY_RUN="--dry-run" ;;
        *) EXTRA_ARGS+=("$1") ;;
    esac
    shift
done

if [[ -z "$SEGMENT" ]]; then
    echo "Usage: bash run.sh <1|2a|2b|3|4> [--dry-run]"
    echo "   1 = BAM ingestion (local_bam symlink + merge_bam merge) + QC gate"
    echo "  2a = SNP discovery + GL + relatedness (stops at clone gate)"
    echo "  2b = subset BEAGLE after clone_approve.py"
    echo "  2c = grouped ngsRelate (per species)"
    echo "   3 = LD + PCA + admixture"
    echo "   4 = SAF + SFS + FST"
    exit 1
fi

TIMESTAMP=$(date +%Y%m%d_%H%M%S)
LOG="$PROJECT/logs/segment${SEGMENT}_${TIMESTAMP}.log"
mkdir -p "$PROJECT/logs"

# Segment → Snakefile + targets
case "$SEGMENT" in
    1)
        SNAKEFILE="$PIPELINE/workflow/Snakefile.align"
        TARGETS=""
        SEG_NAME="bam_ingest"
        ;;
    2a)
        SNAKEFILE="$PIPELINE/workflow/Snakefile.snps"
        TARGETS="results/relatedness/clones_report.txt"
        SEG_NAME="snps_phase1"
        ;;
    2b)
        SNAKEFILE="$PIPELINE/workflow/Snakefile.snps"
        TARGETS=""   # full default targets
        SEG_NAME="snps_phase2"
        ;;
    2c)
        SNAKEFILE="$PIPELINE/workflow/Snakefile.snps"
        TARGETS="results/relatedness/grouped/clones_report.txt"
        SEG_NAME="relatedness_grouped"
        ;;
    3)
        SNAKEFILE="$PIPELINE/workflow/Snakefile.structure"
        TARGETS=""
        SEG_NAME="structure"
        ;;
    4)
        SNAKEFILE="$PIPELINE/workflow/Snakefile.diversity"
        TARGETS=""
        SEG_NAME="diversity"
        ;;
esac

echo "============================================================"
echo "  Project:  $PROJECT"
echo "  Segment:  $SEGMENT ($SEG_NAME)"
echo "  Profile:  $PROFILE"
echo "  Snakefile: $SNAKEFILE"
[[ -n "$DRY_RUN" ]] && echo "  MODE: DRY RUN"
echo "============================================================"

module load singularity/3.10.3
cd "$PROJECT"

$SNAKEMAKE \
    --snakefile "$SNAKEFILE" \
    --profile "$PROFILE" \
    --directory "$PROJECT" \
    $DRY_RUN \
    $TARGETS \
    "${EXTRA_ARGS[@]+${EXTRA_ARGS[@]}}" \
    2>&1 | tee "$LOG"

EXIT_CODE="${PIPESTATUS[0]}"

if [[ "$EXIT_CODE" -eq 0 ]]; then
    echo "============================================================"
    echo "  Segment $SEGMENT ($SEG_NAME) COMPLETE"
    echo "============================================================"
    case "$SEGMENT" in
        2a)
            echo ""
            echo "  Next — review SNP filters:"
            echo "    python $PIPELINE/workflow/scripts/filter_select.py"
            echo "    → writes results/angsd/filter_params.yaml"
            echo "    Then re-run: bash run.sh 2a  (to run pass 2 + relatedness)"
            echo ""
            echo "  Then — review clones:"
            echo "    python $PIPELINE/workflow/scripts/clone_approve.py"
            echo "    → writes results/relatedness/unrelated_samples.txt"
            echo "    Then run: bash run.sh 2b"
            ;;
        2b)
            echo ""
            echo "  Segment 2 complete. Run structure: bash run.sh 3"
            ;;
        2c)
            echo ""
            echo "  Grouped ngsRelate complete."
            echo "  Compare: results/relatedness/grouped/clones_report.txt"
            echo "       vs: results/relatedness/clones_report.txt"
            echo "  Review: python $PIPELINE/workflow/scripts/clone_approve.py --results results/relatedness/grouped"
            ;;
        3)
            echo ""
            echo "  Run lineage assignment gate:"
            echo "    python $PIPELINE/workflow/scripts/lineage_assign.py"
            echo "  Then: bash run.sh 4"
            ;;
        4)
            echo "  Run demography: bash run.sh 6"
            ;;
    esac
else
    echo "============================================================"
    echo "  Segment $SEGMENT FAILED (exit $EXIT_CODE) — Log: $LOG"
    echo "============================================================"
fi

exit "$EXIT_CODE"
