#!/bin/bash

# Original metadata file
original_metadata="/rsstu/users/r/rrellan/sara/RNA_Sequencing_raw/BZea_CLY23D1/NVS205B_RellanAlvarez/hannah/metadata_cpp.txt"
# New file to be created, without carriage return characters
cleaned_metadata="/rsstu/users/r/rrellan/sara/RNA_Sequencing_raw/BZea_CLY23D1/NVS205B_RellanAlvarez/hannah/metadata.txt"

# Remove carriage return characters (\r) and write to new file
tr -d '\r' < "$original_metadata" > "$cleaned_metadata"

echo "Carriage return characters removed. Cleaned file created at: $cleaned_metadata"
