#!/bin/bash


output_dir="/path/to/your/data/my_strainphlan/output_t__subsp.infantis_with_probiotics_11.12"

python "/path/to/your/data/plot_tree_graphlan_modified.py" \
    --ifn_tree "${output_dir}/RAxML_bestTree.t__subsp.infantis.StrainPhlAn4.tre.metadata" \
    --colorized_metadata HistoryAtTimePoint \
    -m HistoryAtTimePoint \
    --leaf_marker_size 60 \
    --legend_marker_size 60 \
    --figure_extension .png \
    --dpi 600
    