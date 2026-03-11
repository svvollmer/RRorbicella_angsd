#!/bin/bash
# =============================================================================
# run_pipeline.sh — Run the coral ANGSD pipeline on EC2 and clean up
#
# - Verifies all host tools are present before starting
# - Runs snakemake with the AWS profile
# - On success: syncs results to S3, then stops the instance
# - On failure: syncs logs to S3, then stops the instance
#
# Usage (from pipeline directory):
#   bash scripts/run_pipeline.sh
#   bash scripts/run_pipeline.sh --dry-run     # snakemake dry-run only
# =============================================================================
set -euo pipefail

PIPELINE_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
S3_BUCKET="coral-angsd-728009587639"
S3_RESULTS="s3://${S3_BUCKET}/results-production"
S3_LOGS="s3://${S3_BUCKET}/logs"
SNS_ARN="arn:aws:sns:us-east-1:728009587639:coral-pipeline-alerts"
LOG_FILE="$PIPELINE_DIR/logs/snakemake.log"
DRY_RUN=0

[[ "${1:-}" == "--dry-run" ]] && DRY_RUN=1

log()  { echo "[$(date '+%Y-%m-%d %H:%M:%S')] $*" | tee -a "$LOG_FILE"; }
die()  { echo "ERROR: $*" >&2 | tee -a "$LOG_FILE"; }

# ---------------------------------------------------------------------------
# Get EC2 instance ID (IMDSv2)
# ---------------------------------------------------------------------------
get_instance_id() {
    local token
    token=$(curl -s -X PUT "http://169.254.169.254/latest/api/token" \
        -H "X-aws-ec2-metadata-token-ttl-seconds: 60" 2>/dev/null) || true
    if [[ -n "$token" ]]; then
        curl -s -H "X-aws-ec2-metadata-token: $token" \
            "http://169.254.169.254/latest/meta-data/instance-id" 2>/dev/null || echo ""
    fi
}

# ---------------------------------------------------------------------------
# Notify via SNS
# ---------------------------------------------------------------------------
notify() {
    local subject="$1"
    local message="$2"
    aws sns publish \
        --topic-arn "$SNS_ARN" \
        --subject "Coral Pipeline: ${subject}" \
        --message "$message" \
        --region us-east-1 || true
}

# ---------------------------------------------------------------------------
# Detect pipeline stop reason and notify
# ---------------------------------------------------------------------------
check_status() {
    local exit_code="$1"

    if [[ -f "$PIPELINE_DIR/results/qc/qc_distribution.txt" ]] && \
       [[ ! -f "$PIPELINE_DIR/results/qc/samples_approved.txt" ]]; then
        notify "QC Checkpoint — Action Required" \
"Pipeline halted at QC gate. Review QC report then:
  1. Start the instance
  2. ssh -i ~/.ssh/coral-pipeline.pem ubuntu@<ip>
  3. cd coral-angsd-pipeline && python workflow/scripts/qc_approve.py
  4. bash scripts/run_pipeline.sh"

    elif [[ -f "$PIPELINE_DIR/results/relatedness/clones_report.txt" ]] && \
         [[ ! -f "$PIPELINE_DIR/results/relatedness/unrelated_samples.txt" ]]; then
        notify "Clone Checkpoint — Action Required" \
"Pipeline halted at clone/relatedness gate. Review clones report then:
  1. Start the instance
  2. ssh -i ~/.ssh/coral-pipeline.pem ubuntu@<ip>
  3. cd coral-angsd-pipeline && python workflow/scripts/clone_approve.py
  4. bash scripts/run_pipeline.sh"

    elif [[ $exit_code -eq 0 ]]; then
        notify "Pipeline Complete" \
"All steps finished successfully.
Results: s3://${S3_BUCKET}/results-production/"

    else
        notify "Pipeline Failed (exit ${exit_code})" \
"Snakemake exited with error code ${exit_code}.
Check logs: s3://${S3_BUCKET}/logs/"
    fi
}

# ---------------------------------------------------------------------------
# Stop this instance — called on both success and failure
# ---------------------------------------------------------------------------
stop_instance() {
    local instance_id
    instance_id=$(get_instance_id)

    log "Syncing results to S3..."
    aws s3 sync "$PIPELINE_DIR/results/" "${S3_RESULTS}/" \
        --exclude "*.fastq" --exclude "*.fastq.gz" --exclude "tmp/*" || true
    aws s3 cp "$LOG_FILE" "${S3_LOGS}/snakemake-$(date '+%Y%m%d-%H%M%S').log" || true

    if [[ -n "$instance_id" ]]; then
        log "Stopping instance ${instance_id}"
        aws ec2 stop-instances --instance-ids "$instance_id" \
            --region us-east-1 > /dev/null || true
    else
        log "WARNING: Could not retrieve instance ID — stop it manually to avoid charges."
    fi
}

# ---------------------------------------------------------------------------
# Pre-flight checks
# ---------------------------------------------------------------------------
preflight() {
    log "Running pre-flight checks..."
    local fail=0

    check() {
        local label="$1"; shift
        if "$@" &>/dev/null; then
            log "  [OK]  $label"
        else
            log "  [FAIL] $label"
            fail=1
        fi
    }

    check "singularity/apptainer"   command -v singularity
    check "snakemake"                command -v snakemake
    check "fasterq-dump"             command -v fasterq-dump
    check "ngsLD"                    command -v ngsLD
    check "aws cli"                  aws --version
    check "S3 bucket accessible"     aws s3 ls "s3://${S3_BUCKET}/"
    check "reference genome"         test -f "$PIPELINE_DIR/reference/apalmata_genome.fasta"
    check "samples CSV"              test -f "$PIPELINE_DIR/config/samples.csv"
    check "Snakefile"                test -f "$PIPELINE_DIR/workflow/Snakefile"
    check "disk space (>200G free)" \
        bash -c '[[ $(df --output=avail / | tail -1) -gt 209715200 ]]'

    if [[ $fail -ne 0 ]]; then
        log "Pre-flight FAILED — fix the issues above. Stopping instance."
        stop_instance "preflight-failure"
        exit 1
    fi

    log "Pre-flight checks passed."
}

# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------
cd "$PIPELINE_DIR"
mkdir -p logs

# Activate conda
source "$HOME/miniforge3/etc/profile.d/conda.sh"
conda activate base
export PATH="/usr/local/bin:$HOME/miniforge3/bin:$PATH"

log "========================================================"
log "Coral ANGSD pipeline starting"
log "Instance: $(get_instance_id)"
log "Dry run: $DRY_RUN"
log "========================================================"

preflight

# Build snakemake command
# - local_conda_env="" disables conda PATH injection (use Singularity on AWS)
# - samples_csv overrides whatever is in config.yaml
SNAKEMAKE_CMD=(
    snakemake
    --snakefile workflow/Snakefile
    --profile profiles/aws
    --config local_conda_env="" samples_csv=config/samples.csv
)
[[ $DRY_RUN -eq 1 ]] && SNAKEMAKE_CMD+=(--dry-run)

log "Running: ${SNAKEMAKE_CMD[*]}"

if [[ $DRY_RUN -eq 1 ]]; then
    "${SNAKEMAKE_CMD[@]}" >> "$LOG_FILE" 2>&1 || true
    log "Dry run complete — instance will not be stopped."
    exit 0
fi

# ---------------------------------------------------------------------------
# Run pipeline with auto-retry on transient failures
# Snakemake caches completed jobs, so each restart picks up where it left off.
# Gates (QC, clone) break the loop so the instance can be stopped for review.
# ---------------------------------------------------------------------------
MAX_RETRIES=20
RETRY_WAIT=300  # seconds between retries
attempt=0
EXIT_CODE=0

while true; do
    if [[ $attempt -gt 0 ]]; then
        log "Retry $attempt/$MAX_RETRIES — waiting ${RETRY_WAIT}s before restarting Snakemake..."
        sleep "$RETRY_WAIT"
    fi

    EXIT_CODE=0
    "${SNAKEMAKE_CMD[@]}" >> "$LOG_FILE" 2>&1 || EXIT_CODE=$?

    # Success — done
    [[ $EXIT_CODE -eq 0 ]] && break

    # QC gate — stop for human review
    if [[ -f "$PIPELINE_DIR/results/qc/qc_distribution.txt" ]] && \
       [[ ! -f "$PIPELINE_DIR/results/qc/samples_approved.txt" ]]; then
        log "QC gate reached — stopping instance for review."
        break
    fi

    # Clone gate — stop for human review
    if [[ -f "$PIPELINE_DIR/results/relatedness/clones_report.txt" ]] && \
       [[ ! -f "$PIPELINE_DIR/results/relatedness/unrelated_samples.txt" ]]; then
        log "Clone gate reached — stopping instance for review."
        break
    fi

    attempt=$((attempt + 1))
    if [[ $attempt -gt $MAX_RETRIES ]]; then
        log "All $MAX_RETRIES retries exhausted — giving up."
        break
    fi

    log "Snakemake exited with code $EXIT_CODE (attempt $attempt/$MAX_RETRIES). Retrying..."
done

check_status "$EXIT_CODE"
stop_instance
[[ $EXIT_CODE -eq 0 ]] || exit $EXIT_CODE
