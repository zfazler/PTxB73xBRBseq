#!/bin/bash

# activate conda environment (samtools needs to be on PATH)
module load conda
source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate /usr/local/usrapps/maize/hdpil/hdpil

# Check if a species name is provided
if [ -z "$1" ]; then
    echo "Usage: $0 <Species Name>"
    exit 1
fi

species=$1

# Define paths based on species
rawDataDIR="/rsstu/users/r/rrellan/sara/RNA_Sequencing_raw/BZea_CLY23D1/NVS205B_RellanAlvarez/hannah/raw_reads/"
inDIR="/rsstu/users/r/rrellan/sara/RNA_Sequencing_raw/BZea_CLY23D1/NVS205B_RellanAlvarez/hannah/clean_reads/"
outDIR="/rsstu/users/r/rrellan/sara/RNA_Sequencing_raw/BZea_CLY23D1/NVS205B_RellanAlvarez/hannah/alignments/"
metadata="/rsstu/users/r/rrellan/sara/RNA_Sequencing_raw/BZea_CLY23D1/NVS205B_RellanAlvarez/hannah/metadata.txt"
summary_file="/rsstu/users/r/rrellan/sara/RNA_Sequencing_raw/BZea_CLY23D1/NVS205B_RellanAlvarez/hannah/${species}/summary_statistics_${species}.txt"

# Create the summary file if it doesn't exist
if [ ! -f "$summary_file" ]; then
    echo -e "SampleID\tRawReads\tCleanReads\tMappedReads\tUniqMappedReads" > "$summary_file"
fi

# Iterate through the metadata lines and extract sample IDs for the given species
awk -v sp="$species" -F'\t' '$2 == sp {print $1}' "$metadata" | while read -r sample_id; do
    # Skip empty sample IDs
    if [ -z "${sample_id}" ]; then
        echo "Sample ID is EMPTY, skipping..."
        continue
    fi

    # Check if the sample is already in the summary file
    if grep -q "^${sample_id}\b" "$summary_file"; then
        echo "Sample ${sample_id} already processed, skipping..."
        continue
    fi

    # Define file paths based on sample ID
    raw_file="${rawDataDIR}${sample_id}.fastq.gz"
    bam_file="${outDIR}${sample_id}_Aligned.sortedByCoord.out.bam"
    final_log_file="${outDIR}${sample_id}_Log.final.out"
    stats_file="${outDIR}${sample_id}_Aligned.sortedByCoord.out.bam.stats"

    # Check if the raw FASTQ file exists before proceeding
    if [ -f "$raw_file" ]; then
        # Count raw reads
        rawReads=$(( $(zcat "$raw_file" | wc -l) / 4 ))
    else
        echo "Raw FASTQ file for ${sample_id} does not exist, skipping raw read count..."
        rawReads="N/A"
    fi

    # Generate stats using samtools if BAM file exists
    if [ -f "$bam_file" ]; then
        samtools stats "$bam_file" > "$stats_file"

        # Extract read counts
        cleanReads=$(grep 'Number of input reads' "$final_log_file" | cut -d'|' -f2 | tr -dc '0-9')
        uniqMappedReads=$(grep 'Uniquely mapped reads number' "$final_log_file" | cut -d'|' -f2 | tr -dc '0-9')
        mappedReads=$(grep 'reads mapped:' "$stats_file" | awk '{print $4}' | tr -dc '0-9')

        # Write to summary
       	printf "%s\t%s\t%s\t%s\t%s\n" "${sample_id}" "${rawReads}" "${cleanReads}" "${mappedReads}" "${uniqMappedReads}" >> "$summary_file"
    else
        echo "BAM file for ${sample_id} does not exist, skipping..."
    fi
done

echo "Summary statistics generated/updated in $summary_file"
