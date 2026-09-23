#!/bin/bash
#SBATCH --job-name=trim_pool
#SBATCH --cpus-per-task=16
#SBATCH --mem=32G
#SBATCH --time=08:00:00
#SBATCH --output=logs/trim_pool_%j.out
#SBATCH --error=logs/trim_pool_%j.err 
# ==============================================================================
# 02b -- Pool-level R2 trimming for STARsolo (HPC)
#
# Alithea's STARsolo command (--clipAdapterType CellRanger4) assumes Read 2
# is ~60-90 cycles. Our sequencing was 150 PE, so R2 has ~60-90 bp of
# adapter/polyA overshoot past the useful cDNA that STARsolo's built-in
# clipping isn't tuned for. This step trims R2 with Trimmomatic PE so pairs
# stay in sync, before STARsolo demux/alignment.
#
# Design decisions:
#   - PE mode so paired reads stay lockstep even when reads get dropped
#   - ILLUMINACLIP targets Nextera+TruSeq adapters (superset covers BRB-seq)
#   - R1 is barcode+UMI (14+14 nt): MINLEN 28 keeps every read with a full
#     barcode+UMI intact. R1 will rarely lose bases -- no Illumina adapter
#     inside the first 28 nt.
#   - SLIDINGWINDOW 4:15 for quality trim on R2 3' end
#   - Only writes the paired-output files; unpaired reads dropped
#
# Usage:
#     bash scripts/02b_trim_pools.sh <POOL_N>
# where POOL_N is 1..4.
#
# Inputs:
#   $poolDir/PT_B73_Cold_P{N}_WKDL260013279-1A_253MY2LT4_L1_1.fq.gz
#   $poolDir/PT_B73_Cold_P{N}_WKDL260013279-1A_253MY2LT4_L1_2.fq.gz
#
# Outputs:
#   $outBaseDir/trimmed/pool_N_R1.fastq.gz     paired R1 (essentially unchanged)
#   $outBaseDir/trimmed/pool_N_R2.fastq.gz     paired R2 (trimmed)
#   $outBaseDir/trimmed/pool_N_trim.log        Trimmomatic summary
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
 
outDir="${outBaseDir}/trimmed"
mkdir -p "$outDir"
R1_out="${outDir}/pool_${pool}_R1.fastq.gz"
R2_out="${outDir}/pool_${pool}_R2.fastq.gz"
R1_unp="${outDir}/pool_${pool}_R1.unpaired.fastq.gz"
R2_unp="${outDir}/pool_${pool}_R2.unpaired.fastq.gz"
log="${outDir}/pool_${pool}_trim.log"
 
module load conda
source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate /usr/local/usrapps/maize/zglover/brbseq_environment/env
export _JAVA_OPTIONS="-Xmx28g" 
# Trimmomatic ships adapter FASTAs with the install. Path is fixed by the
# conda env location; confirmed present via `find` pre-flight.
adapter_file="/usr/local/usrapps/maize/zglover/brbseq_environment/env/share/trimmomatic-0.41-0/adapters/NexteraPE-PE.fa"
if [ ! -f "$adapter_file" ]; then
  echo "ERROR: adapter file not found at $adapter_file" >&2
  echo "  Search: find /usr/local/usrapps/maize/zglover/brbseq_environment/env -name 'NexteraPE-PE.fa'" >&2
  exit 1
fi
 
echo "=== Trimmomatic PE for pool $pool ==="
echo "adapter file:    $adapter_file"
echo "R1 in:           $R1_in  ($(du -h "$R1_in" | cut -f1))"
echo "R2 in:           $R2_in  ($(du -h "$R2_in" | cut -f1))"
echo "R1 out:          $R1_out"
echo "R2 out:          $R2_out"
 
trimmomatic PE -threads 16 -phred33 \
  "$R1_in" "$R2_in" \
  "$R1_out" "$R1_unp" \
  "$R2_out" "$R2_unp" \
  ILLUMINACLIP:${adapter_file}:2:30:10 \
  LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:28 \
  2>&1 | tee "$log"
 
# clean up the singleton "unpaired" files (STARsolo only needs the paired set)
rm -f "$R1_unp" "$R2_unp"
 
echo ""
echo "=== pool $pool trim done ==="
ls -la "$R1_out" "$R2_out"
echo ""
echo "Summary lines from log:"
grep -E "Input Read Pairs|Both Surviving" "$log" || tail -5 "$log"
