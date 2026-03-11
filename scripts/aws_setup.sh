#!/bin/bash
# =============================================================================
# aws_setup.sh — Bootstrap an EC2 instance for the coral ANGSD pipeline
#
# Tested on: Ubuntu 22.04 LTS, us-east-1
# Recommended: c6i.16xlarge (64 vCPU, 128 GB RAM)
#
# Run once after launch from the pipeline directory:
#   bash scripts/aws_setup.sh
#
# Idempotent — safe to re-run if interrupted.
# =============================================================================
set -euo pipefail

PIPELINE_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
MINIFORGE_DIR="$HOME/miniforge3"
SRA_VERSION="3.3.0"
APPTAINER_VERSION="1.3.4"

log() { echo "[$(date '+%H:%M:%S')] $*"; }
die() { echo "ERROR: $*" >&2; exit 1; }

# ---------------------------------------------------------------------------
# 1. System packages
# ---------------------------------------------------------------------------
log "Installing system packages..."
sudo apt-get update -qq
sudo apt-get install -y \
    build-essential wget curl git \
    libssl-dev libz-dev libgsl-dev \
    uuid-dev libgpgme-dev squashfs-tools \
    libseccomp-dev pkg-config cryptsetup-bin \
    awscli

# ---------------------------------------------------------------------------
# 2. Apptainer (Singularity)
# ---------------------------------------------------------------------------
if ! command -v apptainer &>/dev/null; then
    log "Installing Apptainer ${APPTAINER_VERSION}..."
    wget -q "https://github.com/apptainer/apptainer/releases/download/v${APPTAINER_VERSION}/apptainer_${APPTAINER_VERSION}_amd64.deb" \
        -O /tmp/apptainer.deb
    sudo apt-get install -y /tmp/apptainer.deb
    rm /tmp/apptainer.deb
else
    log "Apptainer already installed: $(apptainer --version)"
fi

# Snakemake calls 'singularity' — alias if needed
if ! command -v singularity &>/dev/null; then
    sudo ln -sf "$(which apptainer)" /usr/local/bin/singularity
fi

# ---------------------------------------------------------------------------
# 3. Miniforge3 + Snakemake
# ---------------------------------------------------------------------------
if [[ ! -f "$MINIFORGE_DIR/bin/conda" ]]; then
    log "Installing Miniforge3..."
    wget -q "https://github.com/conda-forge/miniforge/releases/latest/download/Miniforge3-Linux-x86_64.sh" \
        -O /tmp/miniforge.sh
    bash /tmp/miniforge.sh -b -p "$MINIFORGE_DIR"
    rm /tmp/miniforge.sh
else
    log "Miniforge3 already installed."
fi

source "$MINIFORGE_DIR/etc/profile.d/conda.sh"
conda activate base

if ! command -v snakemake &>/dev/null; then
    log "Installing Snakemake and dependencies..."
    mamba install -y -c conda-forge -c bioconda \
        "snakemake>=9" pandas numpy
else
    log "Snakemake already installed: $(snakemake --version)"
fi

# ---------------------------------------------------------------------------
# 4. SRA toolkit (fasterq-dump)
# Container version has SSL/TLS failures on EC2 — must use host binary.
# ---------------------------------------------------------------------------
if ! command -v fasterq-dump &>/dev/null; then
    log "Installing SRA toolkit ${SRA_VERSION}..."
    wget -q "https://ftp-trace.ncbi.nlm.nih.gov/sra/sdk/${SRA_VERSION}/sratoolkit.${SRA_VERSION}-ubuntu64.tar.gz" \
        -O /tmp/sratoolkit.tar.gz
    tar -xzf /tmp/sratoolkit.tar.gz -C "$HOME/"
    rm /tmp/sratoolkit.tar.gz
    sudo ln -sf "$HOME/sratoolkit.${SRA_VERSION}-ubuntu64/bin/fasterq-dump" /usr/local/bin/fasterq-dump
    sudo ln -sf "$HOME/sratoolkit.${SRA_VERSION}-ubuntu64/bin/vdb-config"   /usr/local/bin/vdb-config
else
    log "fasterq-dump already installed: $(fasterq-dump --version 2>&1 | head -1)"
fi

# Disable local .sra cache — saves disk, routes to NCBI free S3 in us-east-1
mkdir -p ~/.ncbi
cat > ~/.ncbi/user-settings.mkfg << 'EOF'
/repository/user/main/public/cache-disabled = "true"
/repository/remote/cloud/aws/accept-charges = "false"
EOF

# ---------------------------------------------------------------------------
# 5. ngsLD (no biocontainer — compiled from source)
# ---------------------------------------------------------------------------
if ! command -v ngsLD &>/dev/null; then
    log "Compiling ngsLD from source..."
    cd /tmp
    git clone --depth 1 https://github.com/fgvieira/ngsLD.git ngsLD_build
    cd ngsLD_build
    make -j4
    sudo cp ngsLD /usr/local/bin/ngsLD
    cd /tmp && rm -rf ngsLD_build
    log "ngsLD installed."
else
    log "ngsLD already installed at $(which ngsLD)"
fi
cd "$PIPELINE_DIR"

# ---------------------------------------------------------------------------
# 6. Restore pre-existing CRAMs from S3 (skips download/align for those samples)
# ---------------------------------------------------------------------------
S3_BUCKET="coral-angsd-728009587639"
if aws s3 ls "s3://${S3_BUCKET}/results-production/bams/" &>/dev/null; then
    log "Restoring CRAMs from S3..."
    mkdir -p results/bams
    aws s3 sync "s3://${S3_BUCKET}/results-production/bams/" results/bams/ \
        --exclude "*" --include "*.cram" --include "*.cram.crai"
    log "CRAMs restored: $(ls results/bams/*.cram 2>/dev/null | wc -l) files"
else
    log "No production CRAMs in S3 yet — starting fresh."
fi

# ---------------------------------------------------------------------------
# 7. Pre-flight verification
# ---------------------------------------------------------------------------
log "Running pre-flight checks..."
FAIL=0

check() {
    local label="$1"; shift
    if "$@" &>/dev/null; then
        echo "  [OK]  $label"
    else
        echo "  [FAIL] $label"
        FAIL=1
    fi
}

check "apptainer/singularity"  command -v singularity
check "snakemake"               command -v snakemake
check "fasterq-dump"            command -v fasterq-dump
check "ngsLD"                   command -v ngsLD
check "aws cli"                 aws --version
check "S3 bucket accessible"    aws s3 ls "s3://${S3_BUCKET}/"
check "reference genome"        test -f "$PIPELINE_DIR/reference/apalmata_genome.fasta"
check "samples CSV"             test -f "$PIPELINE_DIR/config/samples.csv"
check "disk space (>200G free)" \
    bash -c '[[ $(df --output=avail / | tail -1) -gt 209715200 ]]'

if [[ $FAIL -ne 0 ]]; then
    die "Pre-flight checks failed. Fix the issues above before running the pipeline."
fi

log "Setup complete. Run the pipeline with:"
log "  bash scripts/run_pipeline.sh"
