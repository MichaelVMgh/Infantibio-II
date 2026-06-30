#!/bin/bash

# Set the directories
input_dir="/path/to/your/data"  # Change this to your input directory
metaphlan_db="/path/to/your/data"  # Change this to your MetaPhlAn database directory
output_dir="/path/to/your/data"  # Change this to your output directory

# Create the output directories if they don't exist
mkdir -p "$output_dir/profile_out_probiotics"
mkdir -p "$output_dir/sam_out_probiotics"

# Loop through each R1 file in the input directory
for r1_file in "$input_dir"/*_R1_001.fastq.gz; do
    # Determine the corresponding R2 file
    r2_file="${r1_file/_R1_001.fastq.gz/_R2_001.fastq.gz}"
    
    # Extract the sample ID (assuming it is everything before the first underscore)
    sample_id=$(basename "$r1_file" | cut -d'_' -f1)
    
    # Determine the base sample name
    base_sample=$(basename "$r1_file" _R1_001.fastq.gz)
    
    # Remove the stdin_map.bowtie2out.txt file if it exists
    rm -f stdin_map.bowtie2out.txt
    
    # Check if the R2 file exists
    if [[ -f "$r2_file" ]]; then
        # Define output files
        sam_out="$output_dir/sam_out_probiotics/${base_sample}.sam"
        profile_out="$output_dir/profile_out_probiotics/${base_sample}_profile.txt"
        
        # Run MetaPhlAn command
        zcat "$r1_file" "$r2_file" | \
        metaphlan --bowtie2db "$metaphlan_db" -t rel_ab_w_read_stats --index mpa_vOct22_CHOCOPhlAnSGB_lon_subsp \
        --nproc 24 --input_type fastq -s "$sam_out" > "$profile_out"
        
        echo "Processed $base_sample"
    else
        echo "R2 file for $base_sample not found, skipping"
    fi
done

# Define the profile directory
profile_dir="$output_dir/profile_out_probiotics"

# Clean up the profile files
find "$profile_dir" -name "*_profile.txt" -type f | xargs sed -i -e '/reads processed/d'
find "$profile_dir" -name "*_profile.txt" -type f | xargs sed -i -e '/estimated_reads/d'

# Merge the MetaPhlAn tables
merge_metaphlan_tables.py "$profile_dir"/*_profile.txt > "$output_dir/merged_abundance_table_probiotics.txt"

echo "Merging completed. The merged table is saved at $output_dir/merged_abundance_table_probiotics.txt"