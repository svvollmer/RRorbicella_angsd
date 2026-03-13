#!/usr/bin/env bash
# =============================================================================
# run_segment.sh — Launch a pipeline segment with S3 sync + auto-stop
#
# Usage:
#   bash scripts/run_segment.sh <segment> [--profile aws|local|slurm] [--dry-run]
#
# Segments:
#   1  align     — download/trim/map/QC   (Snakefile.align)
#   2  snps      — SNP discovery/GL       (Snakefile.snps)
#   3  structure — LD/PCA/admixture       (Snakefile.structure)
#   4  diversity — SAF/SFS/FST            (Snakefile.diversity)
#   5  report    — annotation/HTML report (Snakefile.report)
#
# Environment variables:
#   S3_BUCKET   — S3 bucket for sync (default: coral-angsd-728009587639)
#   INSTANCE_ID — EC2 instance ID for auto-stop (auto-detected if on EC2)
# =============================================================================

set -euo pipefail

# ---------------------------------------------------------------------------
# Config
# ---------------------------------------------------------------------------
S3_BUCKET="${S3_BUCKET:-coral-angsd-728009587639}"
PIPELINE_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
LOG_DIR="$PIPELINE_DIR/logs"
PROFILE="aws"
DRY_RUN=""

# ---------------------------------------------------------------------------
# Segment definitions
# ---------------------------------------------------------------------------
declare -A SNAKEFILES=(
    [1]="workflow/Snakefile.align"
    [2]="workflow/Snakefile.snps"
    [3]="workflow/Snakefile.structure"
    [4]="workflow/Snakefile.diversity"
    [5]="workflow/Snakefile.report"
)

declare -A SEGMENT_NAMES=(
    [1]="align"
    [2]="snps"
    [3]="structure"
    [4]="diversity"
    [5]="report"
)

# Files to sync FROM S3 before each segment (in addition to pipeline code)
declare -A S3_SYNC_IN=(
    [1]="reference/"
    [2]="results/bams/ results/qc/ results/filtering/ results/angsd/nonrepeat_sites.txt results/angsd/nonrepeat_sites.txt.bin results/angsd/nonrepeat_sites.txt.idx results/angsd/pass1.mafs.gz results/angsd/pass1_snps.txt results/angsd/pass1_snps.txt.bin results/angsd/pass1_snps.txt.idx results/angsd/depth_thresholds.txt reference/"
    [3]="results/angsd/all.unrelated.beagle.gz results/angsd/all.mafs.gz results/relatedness/unrelated_samples.txt"
    [4]="results/bams/ results/angsd/nonrepeat_sites.txt results/angsd/nonrepeat_sites.txt.bin results/angsd/nonrepeat_sites.txt.idx results/angsd/depth_thresholds.txt results/filtering/ results/relatedness/unrelated_samples.txt results/admixture/lineage_assignments.txt reference/"
    [5]="results/pca/ results/admixture/ results/ld/ results/diversity/ results/fst/ results/qc/ results/filtering/ results/relatedness/ results/heterozygosity/ results/angsd/pass1_snps.txt results/angsd/all.mafs.gz"
)

# ---------------------------------------------------------------------------
# Parse arguments
# ---------------------------------------------------------------------------
SEGMENT=""
while [[ $# -gt 0 ]]; do
    case "$1" in
        1|2|3|4|5) SEGMENT="$1" ;;
        align)      SEGMENT=1 ;;
        snps)       SEGMENT=2 ;;
        structure)  SEGMENT=3 ;;
        diversity)  SEGMENT=4 ;;
        report)     SEGMENT=5 ;;
        --profile)  PROFILE="$2"; shift ;;
        --dry-run)  DRY_RUN="--dry-run" ;;
        *) echo "Unknown argument: $1"; exit 1 ;;
    esac
    shift
done

if [[ -z "$SEGMENT" ]]; then
    echo "Usage: bash scripts/run_segment.sh <1-5|align|snps|structure|diversity|report> [--profile aws|local|slurm] [--dry-run]"
    exit 1
fi

SNAKEFILE="${SNAKEFILES[$SEGMENT]}"
SEG_NAME="${SEGMENT_NAMES[$SEGMENT]}"
TIMESTAMP=$(date +%Y%m%d_%H%M%S)
LOG_FILE="$LOG_DIR/segment${SEGMENT}_${SEG_NAME}_${TIMESTAMP}.log"

mkdir -p "$LOG_DIR"

echo "============================================================"
echo "  Segment $SEGMENT: $SEG_NAME"
echo "  Profile: $PROFILE"
echo "  Snakefile: $SNAKEFILE"
echo "  Log: $LOG_FILE"
echo "============================================================"

# ---------------------------------------------------------------------------
# Detect EC2 instance ID (for auto-stop)
# ---------------------------------------------------------------------------
detect_instance_id() {
    if [[ -n "${INSTANCE_ID:-}" ]]; then
        echo "$INSTANCE_ID"
        return
    fi
    # Try IMDSv2
    TOKEN=$(curl -s -X PUT "http://169.254.169.254/latest/api/token" \
        -H "X-aws-ec2-metadata-token-ttl-seconds: 21600" 2>/dev/null || true)
    if [[ -n "$TOKEN" ]]; then
        curl -s -H "X-aws-ec2-metadata-token: $TOKEN" \
            "http://169.254.169.254/latest/meta-data/instance-id" 2>/dev/null || true
    fi
}

auto_stop() {
    if [[ "$PROFILE" != "aws" ]]; then
        return
    fi
    IID=$(detect_instance_id)
    if [[ -n "$IID" ]]; then
        echo "Auto-stopping instance $IID..."
        aws ec2 stop-instances --instance-ids "$IID" || \
            echo "WARNING: ec2:StopInstances failed — stop instance manually."
    fi
}

# ---------------------------------------------------------------------------
# S3 sync helpers
# ---------------------------------------------------------------------------
s3_sync_in() {
    if [[ "$PROFILE" != "aws" ]]; then
        return
    fi
    echo "--- Syncing inputs from S3 ---"
    local s3_prefix="s3://${S3_BUCKET}/results-production"
    local pipeline_prefix="s3://${S3_BUCKET}/pipeline"

    # Sync pipeline code first
    aws s3 sync "$pipeline_prefix/" "$PIPELINE_DIR/" \
        --exclude '.snakemake/*' --exclude 'results/*' \
        || echo "WARNING: pipeline sync failed"

    # Sync segment-specific inputs
    local items="${S3_SYNC_IN[$SEGMENT]}"
    for item in $items; do
        if [[ "$item" == reference/* ]]; then
            # Reference lives in reference/ prefix directly
            aws s3 sync "s3://${S3_BUCKET}/reference/" "$PIPELINE_DIR/reference/" \
                || echo "WARNING: reference sync failed"
        elif [[ "$item" == */ ]]; then
            # Directory
            local local_dir="$PIPELINE_DIR/$item"
            mkdir -p "$local_dir"
            aws s3 sync "$s3_prefix/$item" "$local_dir" \
                --exclude '*.fastq' --exclude '*.fastq.gz' \
                || echo "WARNING: sync failed for $item"
        else
            # Single file
            local local_path="$PIPELINE_DIR/$item"
            mkdir -p "$(dirname "$local_path")"
            aws s3 cp "$s3_prefix/$item" "$local_path" \
                2>/dev/null || echo "NOTE: $item not in S3 (may not exist yet)"
        fi
    done

    # Fix ANGSD index timestamps to prevent "stale index" error (Bug 22)
    for idx_file in \
        "$PIPELINE_DIR/results/angsd/nonrepeat_sites.txt.bin" \
        "$PIPELINE_DIR/results/angsd/nonrepeat_sites.txt.idx" \
        "$PIPELINE_DIR/results/angsd/pass1_snps.txt.bin" \
        "$PIPELINE_DIR/results/angsd/pass1_snps.txt.idx" \
        "$PIPELINE_DIR/reference/"*.fai; do
        [[ -f "$idx_file" ]] && touch "$idx_file"
    done
    echo "--- S3 sync complete ---"
}

s3_sync_out() {
    if [[ "$PROFILE" != "aws" ]]; then
        return
    fi
    echo "--- Syncing results to S3 ---"
    local s3_prefix="s3://${S3_BUCKET}/results-production"
    aws s3 sync "$PIPELINE_DIR/results/" "$s3_prefix/" \
        --exclude '*.fastq' --exclude '*.fastq.gz' \
        --exclude '*.raw.bam' --exclude '*.dedup.bam' \
        || echo "WARNING: results sync to S3 failed"
    echo "--- S3 sync complete ---"
}

# ---------------------------------------------------------------------------
# Snakemake cores / jobs by segment and profile
# ---------------------------------------------------------------------------
declare -A CORES_AWS=(  [1]=16 [2]=64 [3]=32 [4]=64 [5]=4)
declare -A JOBS_AWS=(   [1]=8  [2]=4  [3]=4  [4]=4  [5]=1)
# SLURM: --cores is for local (run:) rules only; --jobs controls SLURM submissions
declare -A CORES_SLURM=([1]=8  [2]=8  [3]=8  [4]=8  [5]=4)
declare -A JOBS_SLURM=( [1]=40 [2]=20 [3]=20 [4]=40 [5]=4)

if [[ "$PROFILE" == "aws" ]]; then
    CORES="${CORES_AWS[$SEGMENT]}"
    JOBS="${JOBS_AWS[$SEGMENT]}"
    SM_EXTRA="--use-singularity --singularity-args '--bind /home'"
elif [[ "$PROFILE" == "slurm" ]]; then
    CORES="${CORES_SLURM[$SEGMENT]}"
    JOBS="${JOBS_SLURM[$SEGMENT]}"
    SM_EXTRA=""   # singularity args come from profiles/slurm/config.yaml
else
    CORES=16
    JOBS=2
    SM_EXTRA=""
fi

# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------
cd "$PIPELINE_DIR"

s3_sync_in

echo "--- Starting snakemake ---"
snakemake \
    --snakefile "$SNAKEFILE" \
    --profile "profiles/$PROFILE" \
    --cores "$CORES" \
    --jobs "$JOBS" \
    --keep-going \
    --latency-wait 60 \
    --rerun-triggers mtime \
    $SM_EXTRA \
    $DRY_RUN \
    2>&1 | tee "$LOG_FILE"

EXIT_CODE="${PIPESTATUS[0]}"

s3_sync_out

if [[ "$EXIT_CODE" -eq 0 ]]; then
    echo "============================================================"
    echo "  Segment $SEGMENT ($SEG_NAME) COMPLETE"
    echo "============================================================"
    # Print gate instructions for segments that have them
    case "$SEGMENT" in
        1)
            echo ""
            echo "Gate 1 (QC):"
            echo "  python workflow/scripts/qc_approve.py"
            echo "  → writes results/qc/samples_approved.txt"
            echo "  Then: bash scripts/run_segment.sh 2"
            ;;
        2)
            echo ""
            echo "If pass 1 just finished, run the SNP filter gate first:"
            echo "  python workflow/scripts/filter_select.py"
            echo "  → writes results/angsd/filter_params.yaml"
            echo "  Then re-run: bash scripts/run_segment.sh 2"
            echo ""
            echo "If pass 2 finished, run the clone gate:"
            echo "  python workflow/scripts/clone_approve.py"
            echo "  → writes results/relatedness/unrelated_samples.txt"
            echo "  Then re-run: bash scripts/run_segment.sh 2  (to subset beagle)"
            ;;
        3)
            echo ""
            echo "Gate 3 (Lineage):"
            echo "  python workflow/scripts/lineage_assign.py"
            echo "  → shows K log-likelihoods, select K, writes lineage_assignments.txt"
            echo "  Then: bash scripts/run_segment.sh 4"
            ;;
        4)
            echo ""
            echo "Segment 4 complete. Start the report:"
            echo "  bash scripts/run_segment.sh 5"
            ;;
        5)
            echo ""
            echo "Pipeline complete. Report at: results/report.html"
            ;;
    esac
else
    echo "============================================================"
    echo "  Segment $SEGMENT ($SEG_NAME) FAILED (exit $EXIT_CODE)"
    echo "  Log: $LOG_FILE"
    echo "============================================================"
fi

auto_stop
exit "$EXIT_CODE"
