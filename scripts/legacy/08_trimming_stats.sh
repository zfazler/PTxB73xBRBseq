#!/bin/bash

# Directories
raw_dir="/rsstu/users/r/rrellan/sara/RNA_Sequencing_raw/BZea_CLY23D1/NVS205B_RellanAlvarez/hannah/raw_reads"
trimmed_dir="/rsstu/users/r/rrellan/sara/RNA_Sequencing_raw/BZea_CLY23D1/NVS205B_RellanAlvarez/hannah/03_trimmed"
output_file="/rsstu/users/r/rrellan/sara/RNA_Sequencing_raw/BZea_CLY23D1/NVS205B_RellanAlvarez/hannah/trimmed_read_statistics.txt"

# Initialize the output file
echo -e "SampleID\tRawReads\tTrimmedReads\tTrimmedPercentage" > "$output_file"

# Process each raw read file
for raw_file in "$raw_dir"/*.fastq.gz; do
    # Extract sample ID
    sample_id=$(basename "$raw_file" .fastq.gz)
    trimmed_file="${trimmed_dir}/${sample_id}_trimmed.fq.gz"
    
    # Check if both files exist
    if [[ -f "$raw_file" && -f "$trimmed_file" ]]; then
        # Count reads in raw and trimmed files
        raw_reads=$(zcat "$raw_file" | wc -l)
        trimmed_reads=$(zcat "$trimmed_file" | wc -l)
        
        # Convert line counts to read counts (4 lines per read)
        raw_reads=$((raw_reads / 4))
        trimmed_reads=$((trimmed_reads / 4))
        
        # Calculate percentage trimmed
        if [[ $raw_reads -gt 0 ]]; then
            trimmed_percentage=$(awk "BEGIN {printf \"%.2f\", (($raw_reads - $trimmed_reads) / $raw_reads) * 100}")
        else
            trimmed_percentage=0
        fi

        # Append results to the output file
        echo -e "${sample_id}\t${raw_reads}\t${trimmed_reads}\t${trimmed_percentage}" >> "$output_file"
    else
        echo "Missing file(s) for ${sample_id}. Skipping..."
    fi
done

echo "Trimming statistics saved to $output_file."
