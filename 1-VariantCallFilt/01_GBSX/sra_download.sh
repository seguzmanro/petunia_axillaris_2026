#!/bin/bash

# Define the list of accessions file
ACC_LIST="axil_acc.txt"

# Check if the accession list file exists
if [ ! -f "$ACC_LIST" ]; then
    echo "Error: Accession list file '$ACC_LIST' not found."
    exit 1
fi

echo "Starting batch download of SRA sequences from $ACC_LIST..."

# Loop through each accession in the file
while IFS= read -r ACC
do
    # Skip empty lines or lines starting with comments
    if [[ -z "$ACC" || "$ACC" =~ ^# ]]; then
        continue
    fi

    echo "--- Processing $ACC ---"

    # Use prefetch to download the .sra file (handles retries)
    # The file will be saved in the current directory
    prefetch "$ACC" --output-directory $(pwd)

    # Use fasterq-dump to convert the .sra file to FASTQ format
    # --split-files for paired-end data, --gzip for compression
    # The output will be in the current working directory
    /usr/bin/fasterq-dump "$ACC"

    pigz -p 24 "$ACC.fastq"

    # Optional: Remove the original .sra file after conversion to save space
    # The exact path might vary based on your vdb-config settings
    rm -rf "$ACC/"

done < "$ACC_LIST"

echo "Batch download finished."
