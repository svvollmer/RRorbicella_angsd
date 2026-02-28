#!/bin/bash
# =============================================================================
# aws_setup.sh — Bootstrap an EC2 spot instance for the coral ANGSD pipeline
#
# Tested on: Ubuntu 22.04 LTS (ami-0c7217cdde317cfec in us-east-1)
# Recommended instance: c6i.16xlarge (64 vCPU, 128 GB RAM)
#
# Run once after launch:
#   bash scripts/aws_setup.sh
# =============================================================================
set -euo pipefail

echo "=== Installing system packages ==="
sudo apt-get update -qq
sudo apt-get install -y \
    build-essential wget curl git python3-pip \
    libssl-dev uuid-dev libgpgme-dev squashfs-tools \
    libseccomp-dev pkg-config cryptsetup-bin

echo "=== Installing Apptainer (Singularity) ==="
APPTAINER_VERSION=1.3.4
wget -q "https://github.com/apptainer/apptainer/releases/download/v${APPTAINER_VERSION}/apptainer_${APPTAINER_VERSION}_amd64.deb"
sudo dpkg -i "apptainer_${APPTAINER_VERSION}_amd64.deb"
rm "apptainer_${APPTAINER_VERSION}_amd64.deb"
# Alias so Snakemake finds it as 'singularity'
sudo ln -sf "$(which apptainer)" /usr/local/bin/singularity

echo "=== Installing Snakemake ==="
pip3 install --user "snakemake==8.*" pandas numpy

echo "=== Installing SRA toolkit ==="
wget -q https://ftp-trace.ncbi.nlm.nih.gov/sra/sdk/current/sratoolkit.current-ubuntu64.tar.gz
tar -xzf sratoolkit.current-ubuntu64.tar.gz
SRA_DIR=$(ls -d sratoolkit.*-ubuntu64 | head -1)
sudo cp "${SRA_DIR}/bin/fasterq-dump" /usr/local/bin/
sudo cp "${SRA_DIR}/bin/vdb-config"   /usr/local/bin/
rm -rf sratoolkit.*

echo "=== Configuring SRA toolkit for cloud access ==="
# fasterq-dump on EC2 in us-east-1 automatically routes to NCBI's free
# S3 bucket (sra-pub-run-odp) — no download charges.
# Disable local .sra caching to save disk.
mkdir -p ~/.ncbi
cat > ~/.ncbi/user-settings.mkfg << 'EOF'
/repository/user/main/public/cache-disabled = "true"
/repository/remote/cloud/aws/accept-charges = "false"
EOF

echo "=== Setting up /data directory ==="
# If you attached a separate EBS volume (recommended: 600 GB gp3):
#   sudo mkfs.ext4 /dev/nvme1n1
#   sudo mount /dev/nvme1n1 /data
#   echo '/dev/nvme1n1 /data ext4 defaults,nofail 0 2' | sudo tee -a /etc/fstab
# Otherwise work in /home (instance store or root EBS):
sudo mkdir -p /data
sudo chown "$USER" /data

echo ""
echo "=== Setup complete. Next steps: ==="
echo ""
echo "1. Clone the repo:"
echo "   cd /data && git clone https://github.com/svvollmer/coral-angsd-pipeline.git"
echo "   cd coral-angsd-pipeline"
echo ""
echo "2. Download reference genome (public NCBI assembly GCF_964030605.1):"
echo "   cd /data/coral-angsd-pipeline"
echo "   pip install ncbi-datasets-cli"
echo "   datasets download genome accession GCF_964030605.1 \\"
echo "       --include genome,gff3 --filename reference/ncbi_dataset.zip"
echo "   cd reference && unzip -o ncbi_dataset.zip"
echo "   cp ncbi_dataset/data/GCF_964030605.1/GCF_964030605.1_*.fna apalmata_genome.fasta"
echo "   cp ncbi_dataset/data/GCF_964030605.1/genomic.gff apalmata_genes.gff"
echo "   cd .."
echo ""
echo "3. Edit config/config.yaml:"
echo "   - Comment out local_conda_env"
echo "   - Set samples_csv: config/samples.csv"
echo ""
echo "4. Run:"
echo "   ~/.local/bin/snakemake --snakefile workflow/Snakefile --profile profiles/aws"
echo ""
echo "5. If interrupted (spot), relaunch instance and re-run same command."
echo "   rerun-incomplete: true in the profile handles recovery automatically."
