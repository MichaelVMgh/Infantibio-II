# Infantibio-II Shotgun Metagenomics Analysis Pipeline and Downstream Analysis





## Requirements





* KneadData v 0.10.0
* Trim Galore v 0.6.10
* MetaPhlAn 4.1.1
* StrainPhlAn4 v 4.1.1
* HUMAnN v 3.6
* MEGAHIT v 1.2.9
* Prodigal v 2.6.3
* cd-hit v 4.8.1
* MSPminer
* InterProScan v 5.55_88.0





## In-house custom scripts for the reference-based tools used in the pipeline include:



- `metaphlan4_subsp_longum_infantibio.sh`
- `strainphlan_extract_markers.sh`
- `strainphlan_sample2markers.sh`
- `strainphlan_main.sh`
- `strainphlan_add_metadata.sh`
- `strainphlan_graphlan.sh`
- `global_Bifidobacterium_strain_diversity_trees.sh` (for global Bifidobacterium strain diversity analysis with StrainPhlAn)
  - uses the `make_annotation_file_ALL.py` python and `total_metadata_merged_ALL.txt` excel files.
- Additional reference-based tools (KneadData, Trim Galore and HUMAnN) and assembly-based tools (MEGAHIT, Prodigal, cd-hit, MSPMiner, and InterProScan) were executed using an in-house pipeline implemented in Snakemake v 7.26.0.










## Downstream analysis





### Alpha- and Beta- diversities (overall species-level microbiome composition)





For alpha-diversity and beta-diversity analysis we will use the R script named `diversities.R`. For this script we would need to install the R packages "vegan" and "ggplot". We would need two input files, A) a MetaPhlAn output table with relative abundances where the bacterial species are in the columns and the samples are in the rows, loaded as "df" (the `merged_abundance_table_species_transposed.txt` file is our original output table with MetaPhlAn4 species-level taxonomic profiles and can be used for this purpose), and B) a mapping file (text file; `mapping_file.txt`) with all the relevant metadata as columns and the samples in the rows, which we load as "mf".





### Dominant genera and Bifidobacterium species





For visualizing the relative abundances of the dominant Bifidobacterium species by formula group or the stacked barplots with the top bacterial genera we use the Jupyter scripts `Bifidobacterium_boxplots_by_formula.ipynb` or `top_genera_stacked_barplots.ipynb` respectively. We would only need one input file, which would be a MetaPhlAn output table with added metadata where both species names and metadata are in the columns and the samples names are in the rows. For the stacked barplots, a similar table would be needed with the genera abundances. The table included in the `merged_abundance_table_species_transposed.txt` can be used as input for the Bifidobacterium boxplots. If not, your own MetaPhlAn transposed table with the added metadata info as columns can also be used.





### Bray-Curtis dissimilarity heatmap at 24 months





For the Bray-Curtis dissimilarity heatmap at 24 months the Python script `bray_curtis_heatmap.py` can be used with the same `merged_abundance_table_species_transposed.txt` input file as mentioned before.





### Strain Diversity Analysis of Bifidobacterium dominant species in a Global Reference Database





For the global Bifidobacterium strain diversity analysis with StrainPhlAn we use the file named `global_Bifidobacterium_strain_diversity_trees.sh` and follow the steps included there.





### Multivariate association analyses with Maaslin2 for functional profiling (HUMAnN3 data: unstratified; gene families / gene pathways)





For the multivariate association analyses using Maaslin2 (with Humann3 data) the R script named `maaslin.R` is used. For this purpose the Maaslin2 R package would need to be installed with BiocManager. Two input files are required: a table with all the gene families/pathways abundances, and a mapping file. Sample IDs should match between these two files.





### Volcano plots for visualization of significant HUMAnN3 gene families / gene pathways from Maaslin2 output





For the volcanoplots with Humann3 functional profiles (gene families or gene pathways), the R script `volcanoplots.R` is used. The Maaslin2 output file with the significant gene families/pathways (filtered by the associated ones to your desired variable; e.g., feeding) would be used as input file.





### Analysis of top ranked annotated genes (de-novo assembly data) significantly different between BL subsp. infantis cluster 1 and cluster 2 (strain-level analysis)





For the BL subsp. infantis horizontal barplots with the top significant annotated genes, the script `BLinfantis_clusters_top_genes_barplots.py` may be used. Input files should be: A) a Maaslin2 output file with the significant genes associated with the relevant variable (e.g., clade or cluster), and B) an output file from InterproScan with the annotated Pfam domains.





### Visualization of StrainPhlAn4 strain evolutionary distances (for relevant species) from infant to paired parents: related vs unrelated; impact of feeding





For the evolutionary distances we would need the Python script `evolutionary_distances.py` and two input files. The first input file should be a "full symmetric distance matrix" table (e.g., `buniformis_full_symmetric_distance_matrix.txt`) for each bacterial species where sample names are included as row names AND columns. Values are evolutionary distances between samples according to the StrainPhlAn distance matrix output (Kimura-corrected). The second file (e.g., `buniformis_with_distances_to_parents.txt`) would be a table with the samples as rows and the "Distances_to_mother" and "Distances_to_father" as two separate columns with the evolutionary distances values from each infant sample to the respective mother or father. 





## Citation





If you use this code, please cite:





**Preprint:**





Merino MV, Weißenberger S, Ermolova A, Dhilly E, Spadazzi R, Levy L, Collado MC, Hertz T, Omer H, Heidrich V, Segata N, Larsen M, Schirmer M, Haller D. Probiotic formula intervention in infants leads to colonization and competitive strain displacement independent of IgA binding. bioRxiv (2026). https://doi.org/10.64898/2026.06.03.729775





