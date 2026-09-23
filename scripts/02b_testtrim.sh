#!/bin/bash
#SBATCH --job-name=trim_pool1_test
#SBATCH --cpus-per-task=16
#SBATCH --mem=32G
#SBATCH --time=08:00:00
#SBATCH --output=logs/trim_pool1_%j.out
#SBATCH --error=logs/trim_pool1_%j.err
 
# ==============================================================================
# Test submission: pool 1 only, for validating 02b_trim_pools.sh before
# committing to the larger pools (3 and 4 are much bigger and will take
# proportionally longer).
#
# Usage:
#   sbatch scripts/test_trim_pool1.sh
#
# Check status:
#   squeue -u $USER
#
# After it finishes, check:
#   - logs/trim_pool1_<jobid>.out / .err  (SLURM stdout/stderr)
#   - PTxB73xBRBseq/trimmed/pool_1_trim.log  (Trimmomatic's own summary)
#   - PTxB73xBRBseq/trimmed/pool_1_R1.fastq.gz, pool_1_R2.fastq.gz (outputs)
# ==============================================================================
 
set -euo pipefail
 
mkdir -p logs
 
bash scripts/02b_trim_pools.sh 1
 
