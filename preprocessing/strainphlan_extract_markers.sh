#!/bin/bash

# Define directories
metaphlan_dir="/path/to/your/data/my_metaphlan"
sam_dir="$metaphlan_dir/sam_out"
mkdir -p "$metaphlan_dir/compressed_sam_out" #create if it does not exist
compressed_sam_dir="$metaphlan_dir/compressed_sam_out"
mkdir -p "/path/to/your/data/my_strainphlan" #create if it does not exist
output_dir="/path/to/your/data/my_strainphlan"
mkdir -p "$output_dir/consensus_markers" #create if it does not exist
strainphlan_dir="/path/to/your/data/strainphlan"
metaphlan_db_dir="/path/to/your/data/modified_database_with_blongum_markers/mpa_vOct22_CHOCOPhlAnSGB_lon_subsp.pkl"
#cur_clade=t__SGB6939
#cur_clade=t__SGB6952 #V.dispar
cur_clade=t__SGB17256 #B.bifidum

# Compress each .sam file to .sam.bz2 file and move to compressed_sam_dir
#for sam_file in "$sam_dir"/*.sam; do
    #echo "Compressing $sam_file..."
    #bzip2 "$sam_file"
    #mv "${sam_file}.bz2" "$compressed_sam_dir"
    #echo "Moved ${sam_file}.bz2 to $compressed_sam_dir"

    # Convert .sam to .bam
    #bam_file="${sam_file%.sam}.bam"
    #echo "Converting $sam_file to $bam_file..."
    #samtools view -S -b "$sam_file" > "$bam_file"
    #echo "Converted $sam_file to $bam_file"
    #mv "${bam_file}.bam" "$bam_dir"
    #echo "Moved ${bam_file}.bam to $bam_dir"
#done


#GET MARKERS FOR THE CLADE FROM METAPHLAN DATABASE

mkdir -p "$output_dir/clade_markers"
markers_dir="$output_dir/clade_markers"

python "/path/to/your/data/extract_markers.py" -c ${cur_clade} -o ${markers_dir} -d ${metaphlan_db_dir}