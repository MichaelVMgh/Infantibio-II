#!/bin/bash

output_dir="/path/to/your/data"

python "/path/to/your/data/add_metadata_tree.py" \
    -t "${output_dir}/RAxML_bestTree.t__subsp.infantis.StrainPhlAn4.tre" \
    -f "/path/to/your/data/mapping_file_Blongum_all-tree-samples.txt" \
    -m HistoryAtTimePoint \
    --string_to_remove .fastq.bz2