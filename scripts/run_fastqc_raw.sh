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

mkdir -p /rsstu/users/r/rrellan/CERCA-Cold/PTxB73xBRBseq/QC/raw_trimming

RAW_DIR="/rsstu/users/r/rrellan/CERCA-Cold/01_Incoming/Novogene_B73xPT_Cold/X202SC26083653-Z01-F001/01.RawData"
QC_DIR="/rsstu/users/r/rrellan/CERCA-Cold/PTxB73xBRBseq/QC/raw_trimming"


POOL=$SLURM_ARRAY_TASK_ID

POOL_DIR="$RAW_DIR/PT_B73_Cold_P${POOL}"

R1="$POOL_DIR/PT_B73_Cold_P${POOL}_WKDL260013279-1A_253MY2LT4_L1_1.fq.gz"
R2="$POOL_DIR/PT_B73_Cold_P${POOL}_WKDL260013279-1A_253MY2LT4_L1_2.fq.gz"

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
