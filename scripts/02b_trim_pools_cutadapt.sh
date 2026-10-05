#!/bin/bash
#SBATCH --job-name=trim_pool
#SBATCH --cpus-per-task=16
#SBATCH --mem=32G
#SBATCH --time=08:00:00
#SBATCH --output=logs/trim_pool_%j.out
#SBATCH --error=logs/trim_pool_%j.err
# ==============================================================================
# 02b -- Pool-level R2 trimming for STARsolo (HPC)  [v2: cutadapt, R2-only]
#
# Alithea's STARsolo command (--clipAdapterType CellRanger4) assumes Read 2
# is ~60-90 cycles. Our sequencing was 150 PE, so R2 has adapter/polyA/polyG
# overshoot past the useful cDNA. This step trims R2 only, keeping pairs in
# sync, before STARsolo demux/alignment.
#
# Why v2 (replaces Trimmomatic v1):
#   - v1 used NexteraPE-PE.fa, which contains NO TruSeq sequence. BRB-seq R2
#     reads into the P5 side = TruSeq Read 1 adapter (AGATCGGAAGAGC...), so
#     adapter was never removed (~20% residual in pool 4 post-trim FastQC).
#   - v1 LEADING/TRAILING/SLIDINGWINDOW also ran on R1. That truncated R1
#     across the polyT, and LEADING:3 could clip base 1 and silently shift
#     the barcode/UMI frame.
#   - Trimmomatic has no polyA / polyG trimming.
#
# Design decisions (v2):
#   - cutadapt uppercase options (-A, -Q) act on R2 ONLY; R1 (14 nt barcode +
#     14 nt UMI + polyT + cDNA) passes through untouched.
#   - R2 read structure on short inserts:
#         cDNA -> polyA -> UMI/BC (rc, 28 nt) -> TruSeq R1 adapter -> polyG
#     -A TruSeq adapter  : removes adapter and everything after it
#     -A A{20}           : removes from polyA tail onward (also drops the
#                          BC/UMI/adapter that follow on short inserts)
#     -A G{20}           : removes polyG no-signal runs (2-color chemistry)
#     -Q 20              : 3' quality trim, R2 only
#     -n 2               : up to 2 trimming rounds per read
#   - -m 28:20 + --pair-filter=any : drop pair if R1 < 28 (full BC+UMI) or
#     R2 < 20 (no usable cDNA, e.g. pool 4 primer/polyA artifacts)
#   - Inputs are the RAW Novogene fastqs (v1 trimmed R1 is not reusable).
#   - Outputs go to trimmed_cutadapt/ so v1 outputs are not overwritten.
#   - STARsolo downstream: set --clipAdapterType None (R2 already trimmed).
#
# Usage:
#     sbatch scripts/02b_trim_pools.sh <POOL_N>      (POOL_N = 1..4)
#
# Inputs:
#   $poolDir/PT_B73_Cold_P{N}_WKDL260013279-1A_253MY2LT4_L1_1.fq.gz
#   $poolDir/PT_B73_Cold_P{N}_WKDL260013279-1A_253MY2LT4_L1_2.fq.gz
#
# Outputs:
#   $outDir/pool_N_R1.fastq.gz        paired R1 (untrimmed, filtered only)
#   $outDir/pool_N_R2.fastq.gz        paired R2 (trimmed)
#   $outDir/pool_N_trim.log           cutadapt text report
#   $outDir/pool_N_cutadapt.json      cutadapt JSON report (MultiQC-readable)
#   $outDir/fastqc/                   post-trim FastQC (if fastqc available)
# ==============================================================================

set -euo pipefail

pool="${1:-}"
if [[ ! "$pool" =~ ^[1-4]$ ]]; then
  echo "usage: $0 <POOL_N>   (POOL_N = 1, 2, 3, or 4)" >&2
  exit 1
fi

rawDataDir="/rsstu/users/r/rrellan/CERCA-Cold/01_Incoming/Novogene_B73xPT_Cold/X202SC26083653-Z01-F001/01.RawData"
poolDir="${rawDataDir}/PT_B73_Cold_P${pool}"

outBaseDir="/rsstu/users/r/rrellan/CERCA-Cold/PTxB73xBRBseq"

R1_in="${poolDir}/PT_B73_Cold_P${pool}_WKDL260013279-1A_253MY2LT4_L1_1.fq.gz"
R2_in="${poolDir}/PT_B73_Cold_P${pool}_WKDL260013279-1A_253MY2LT4_L1_2.fq.gz"

outDir="${outBaseDir}/trimmed_cutadapt"
mkdir -p "$outDir"
R1_out="${outDir}/pool_${pool}_R1.fastq.gz"
R2_out="${outDir}/pool_${pool}_R2.fastq.gz"
log="${outDir}/pool_${pool}_trim.log"
json="${outDir}/pool_${pool}_cutadapt.json"

# TruSeq Read 1 adapter as it appears in BRB-seq R2 (orientation confirmed by
# grep on raw R2: pool1 7.6%, pool2 16%, pool3 23%, pool4 63% of first 1M reads)
truseq_r2="AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT"

for f in "$R1_in" "$R2_in"; do
  [ -f "$f" ] || { echo "ERROR: input not found: $f" >&2; exit 1; }
done

module load conda
source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate /usr/local/usrapps/maize/zglover/brbseq_environment/env

# cutadapt >= 4.0 needed for -m X:Y and --json
if ! command -v cutadapt >/dev/null 2>&1; then
  echo "ERROR: cutadapt not in env. Install with:" >&2
  echo "  conda install -p /usr/local/usrapps/maize/zglover/brbseq_environment/env -c bioconda 'cutadapt>=4.4'" >&2
  exit 1
fi
ca_ver="$(cutadapt --version)"
if (( ${ca_ver%%.*} < 4 )); then
  echo "ERROR: cutadapt $ca_ver found; need >= 4.0" >&2
  exit 1
fi

echo "=== cutadapt (R2-only trimming) for pool $pool ==="
echo "cutadapt:        $ca_ver"
echo "R1 in:           $R1_in  ($(du -h "$R1_in" | cut -f1))"
echo "R2 in:           $R2_in  ($(du -h "$R2_in" | cut -f1))"
echo "R1 out:          $R1_out"
echo "R2 out:          $R2_out"

cutadapt -j "${SLURM_CPUS_PER_TASK:-16}" \
  -A "$truseq_r2" \
  -A "A{20}" \
  -A "G{20}" \
  -Q 20 \
  -n 2 \
  -m 28:20 \
  --pair-filter=any \
  --json="$json" \
  -o "$R1_out" -p "$R2_out" \
  "$R1_in" "$R2_in" \
  2>&1 | tee "$log"

echo ""
echo "=== pool $pool trim done ==="
ls -la "$R1_out" "$R2_out"
echo ""
echo "Summary lines from log:"
grep -E "Total read pairs processed|Read 2 with adapter|Pairs that were too short|Pairs written" "$log" || tail -20 "$log"

# Post-trim FastQC (action item: confirm R1 untouched, R2 adapter/polyA/polyG ~0)
if command -v fastqc >/dev/null 2>&1; then
  mkdir -p "${outDir}/fastqc"
  echo ""
  echo "=== post-trim FastQC ==="
  fastqc -t 2 -o "${outDir}/fastqc" "$R1_out" "$R2_out"
else
  echo "NOTE: fastqc not in env; run FastQC on $R1_out and $R2_out separately."
fi
