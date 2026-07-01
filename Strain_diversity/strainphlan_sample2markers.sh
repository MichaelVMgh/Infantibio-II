#!/bin/bash

# Define directories
metaphlan_dir="/path/to/your/data/my_metaphlan"
sam_dir="$metaphlan_dir/sam_out"
mkdir -p "$metaphlan_dir/compressed_sam_out" #create if it does not exist
compressed_sam_dir="$metaphlan_dir/compressed_sam_out"
mkdir -p "/path/to/your/data/strainphlan_moran" #create if it does not exist
output_dir="/path/to/your/data/strainphlan_moran"
mkdir -p "$output_dir/consensus_markers_moran_with_infantibio" #create if it does not exist
metaphlan_db_dir="/path/to/your/data/modified_database_with_blongum_markers/mpa_vOct22_CHOCOPhlAnSGB_lon_subsp.pkl"
#cur_clade=s__Bifidobacterium_longum


# Compress each .sam file to .sam.bz2 file and move to compressed_sam_dir
#for sam_file in "$sam_dir"/*.sam; do
#    echo "Compressing $sam_file..."
#    bzip2 "$sam_file"
#    mv "${sam_file}.bz2" "$compressed_sam_dir"
#    echo "Moved ${sam_file}.bz2 to $compressed_sam_dir"

    # Convert .sam to .bam
#    bam_file="${sam_file%.sam}.bam"
#    echo "Converting $sam_file to $bam_file..."
#    samtools view -S -b "$sam_file" > "$bam_file"
#    echo "Converted $sam_file to $bam_file"
#    mv "${bam_file}.bam" "$bam_dir"
#    echo "Moved ${bam_file}.bam to $bam_dir"
#done

#GET CONSENSUS MARKERS

# Generate a marker file for each sample
for bz2_file in "$compressed_sam_dir"/*.sam.bz2; do
    if [[ -f "$bz2_file" ]]; then
        echo "Processing $bz2_file with sample2markers.py..."

        # Run the sample2markers.py script
        python "/path/to/your/data/MetaPhlAn/metaphlan/utils/sample2markers.py" \
            -i "$bz2_file" \
            --nproc 20 \
            -o "$output_dir/consensus_markers_moran_with_infantibio/" \
            -d "$metaphlan_db_dir"

        echo "Generated marker file for $bz2_file"
    else
        echo "No .sam.bz2 files found in $compressed_sam_dir"
    fi
done





