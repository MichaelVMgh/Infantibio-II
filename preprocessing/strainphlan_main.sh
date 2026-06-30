#!/bin/bash

# Define directories
metaphlan_db_dir="/path/to/your/data/modified_database_with_blongum_markers/mpa_vOct22_CHOCOPhlAnSGB_lon_subsp.pkl"
consensus_markers_dir="/path/to/your/data/my_strainphlan/consensus_markers"
clade_markers_dir="/path/to/your/data/my_strainphlan/clade_markers"

mkdir -p "/path/to/your/data/my_strainphlan/output_t__SGB6939_vparvula"
output_dir="/path/to/your/data/my_strainphlan/output_t__SGB6939_vparvula"

#cur_clade="t__SGB17256" #B.bifidum
cur_clade="t__SGB6939" #V.parvula

# Generate trees from alignments
#for bz2_file in "${consensus_markers_dir}"/*.json.bz2; do
    #if [[ -f "$bz2_file" ]]; then
        #echo "Processing $bz2_file with strainphlan..."

        # Remove the folder if it exists
        #rm -rf "${output_dir}/${cur_clade}_mutation_rates/mutation_rates"


        # Run StrainPhlAn
        #python "/nfs/data3/michael_metagenomics/MetaPhlAn/metaphlan/strainphlan.py" \
        #    -s "${consensus_markers_dir}"/* \
        #    -m "${clade_markers_dir}/${cur_clade}.fna" \
        #    -o "${output_dir}" \
        #    -n 4 \
        #    -c  ${cur_clade} \
        #    -d "$metaphlan_db_dir" \
        #    --phylophlan_mode accurate \
        #    --mutation_rates

        #echo "Generated output file for $bz2_file"
    #else
        #echo "No .json.bz2 files found in $consensus_markers_dir"
    #fi
#done

#strainphlan -s consensus_markers/* -m clade_markers/t__SGB4933.fna -r reference_genomes/*.fna.bz2 -o output -c t__SGB4933 --phylophlan_mode fast --nproc 4
#strainphlan -s  ${consensus_markers}/*.pkl -m ${marker_dir}/${cur_clade}.fna -o ${output_dir} -n 3 -c  ${cur_clade}  -r  ${ref_genome} --phylophlan_mode accurate --mutation_rates

#Run StrainPhlAn
python "/path/to/your/data/strainphlan.py" \
    -s ${consensus_markers_dir}/*.json.bz2 \
    -m "${clade_markers_dir}/${cur_clade}.fna" \
    -o "${output_dir}" \
    -n 4 \
    -c  ${cur_clade} \
    -d "$metaphlan_db_dir" \
    --phylophlan_mode accurate \
    --mutation_rates