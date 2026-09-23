#!/bin/bash
#SBATCH --job-name=fastqc_pool
#SBATCH --cpus-per-task=8
#SBATCH --mem=16G
#SBATCH --time=02:00:00
#SBATCH --output=/rsstu/users/r/rrellan/CERCA-Cold/PTxB73xBRBseq/logs/fastqc_pool_%A_%a.out
#SBATCH --error=/rsstu/users/r/rrellan/CERCA-Cold/PTxB73xBRBseq/logs/fastqc_pool_%A_%a.err
#SBATCH --array=1-4

# ============================================================
# FastQC on trimmed BRB-seq reads, one pool at a time
# Pools 1-4
# ============================================================

source /usr/local/apps/conda/miniconda3/26.3.2/etc/profile.d/conda.sh
conda activate /rsstu/users/r/rrellan/sara/nirwan_backup/ntanduk/seqanal

echo "=========================================="
echo "FastQC - Pool $SLURM_ARRAY_TASK_ID"
echo "=========================================="

echo "FastQC:"
which fastqc
fastqc --version

TRIMMED_DIR="/rsstu/users/r/rrellan/CERCA-Cold/PTxB73xBRBseq/trimmed"
QC_DIR="/rsstu/users/r/rrellan/CERCA-Cold/PTxB73xBRBseq/QC/after_trimming"

mkdir -p "$QC_DIR"

POOL=$SLURM_ARRAY_TASK_ID

R1="$TRIMMED_DIR/pool_${POOL}_R1.fastq.gz"
R2="$TRIMMED_DIR/pool_${POOL}_R2.fastq.gz"

echo "R1: $R1"
echo "R2: $R2"

# Check that input files exist
if [ ! -f "$R1" ] || [ ! -f "$R2" ]; then
    echo "ERROR: Input FASTQ files not found."
    exit 1
fi

echo "Running FastQC..."

fastqc \
    -t 8 \
    -o "$QC_DIR" \
    "$R1" \
    "$R2"

echo "=========================================="
echo "FastQC completed for Pool $POOL"
echo "Output directory:"
echo "$QC_DIR"
echo "=========================================="
