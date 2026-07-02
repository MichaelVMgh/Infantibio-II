# Infantibio-II Shotgun Metagenomics Analysis Pipeline and Downstream Analysis





## Requirements





### Reference-based tools:





* KneadData v 0.10.0
* Trim Galore v 0.6.10
* MetaPhlAn v 4.1.1
* StrainPhlAn v 4.1.1
* HUMAnN v 3.6 (with MetaPhlAn 3.1.0 output)





### Assembly-based tools:





* MEGAHIT v 1.2.9
* Prodigal v 2.6.3
* cd-hit v 4.8.1
* MSPminer
* InterProScan v 5.55_88.0





## In-house custom scripts for the tools used in the pipeline include:



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





For alpha-diversity and beta-diversity analysis we will use the R script named `Diversities.R`. For this script we would need to install the R packages "vegan" and "ggplot". We would need two input files, A) a MetaPhlAn output table with relative abundances where the bacterial species are in the columns and the samples are in the rows (the `merged_abundance_table_species_transposed.txt` file is our original output table with MetaPhlAn4 species-level taxonomic profiles and can be used for this purpose), and B) a mapping file (`mapping_file.txt`) with all the relevant metadata as columns and the samples in the rows.





### Dominant genera (genus-level)





For visualizing the median relative abundances of top bacterial genera as stacked barplots we use the Jupyter script `top_genera_stacked_barplots.ipynb`. We would only need one input file, which would be a MetaPhlAn output table (transposed) with the taxonomic relative abundances at the genus-level as values, the bacterial genera as columns and the samples as rows. The table included in the `merged_abundance_table_genus_transposed.txt` can be used. If not, your own MetaPhlAn transposed table with the added metadata info as columns can also be used.





### Top Bifidobacterium species





For generating the relative abundances of the dominant Bifidobacterium species by formula group we use the Jupyter script `Bifidobacterium_boxplots_by_formula.ipynb` file. As previously, we would need the `merged_abundance_table_species_transposed.txt` as input file (table with the species-level relative abundances). This table would include both species names and metadata in the columns and the samples names in the rows.





### Bray-Curtis dissimilarity heatmap at 24 months





For the Bray-Curtis dissimilarity heatmap at 24 months the Python script `bray_curtis_heatmap.py` can be used with the same `merged_abundance_table_species_transposed.txt` input file as mentioned before.





### Strain Diversity Analysis





To generate these trees we would use the StrainPhlAn tool with each step in each separated script. We start with `strainphlan_extract_markers.sh` to extract reference species-specific marker genes (from the desired clade) from the MetaPhlAn database and `strainphlan_sample2markers.sh` to reconstruct consensus marker sequences for every species/strain detected in each metagenomic sample. The `strainphlan_main.sh` is the main script and builds the phylogenetic tree. Then, the script `strainphlan_add_metadata.sh` adds sample metadata to the tree while `strainphlan_graphlan.sh` is used for visualization and for generating the final figures.





### Strain Diversity Analysis of Bifidobacterium dominant species in a Global Reference Database





For the global Bifidobacterium strain diversity analysis with StrainPhlAn we use the file named `global_Bifidobacterium_strain_diversity_trees.sh` and follow the steps included there (using the StrainPhlAn tool, similarly to the step before).





### Multivariate association analyses with Maaslin2 for functional profiling (HUMAnN3 data: unstratified; gene families / gene pathways)





For the multivariate association analyses using Maaslin2 (with Humann3 data) the R script named `maaslin.R` is used. For this purpose the Maaslin2 R package would need to be installed with BiocManager. Two input files are required: a table with all the gene families/pathways abundances (`genepathways_unstratified.txt` or `genefamilies_unstratified`), and a mapping file. Sample IDs should match between these two files.





### Volcano plots for visualization of significant HUMAnN3 gene families / gene pathways from Maaslin2 output





For the volcanoplots with Humann3 functional profiles (gene families or gene pathways), the R script `volcanoplots.R` is used. The Maaslin2 output file with the significant gene families/pathways (filtered by the associated ones to your desired variable; e.g., feeding) would be used as input file.





### Analysis of top ranked annotated genes (de-novo assembly data) significantly different between BL subsp. infantis cluster 1 and cluster 2 (strain-level analysis)





For the BL subsp. infantis horizontal barplots with the top significant annotated genes, the script `BLinfantis_clusters_top_genes_barplots.py` may be used. Input files should be: A) a Maaslin2 output file with the significant genes associated with the relevant variable (e.g., clade or cluster), and B) an output file from InterproScan with the annotated Pfam domains.





### Visualization of StrainPhlAn4 strain evolutionary distances (for selected species) from infant to paired parents: related vs unrelated; impact of feeding





For the evolutionary distances we would need the Python script `evolutionary_distances.py` and two input files. The first input file should be a "full symmetric distance matrix" table (e.g., `buniformis_full_symmetric_distance_matrix.txt`) for each bacterial species where sample names are included as row names AND columns. Values are evolutionary distances between samples according to the StrainPhlAn distance matrix output (Kimura-corrected). The second file (e.g., `buniformis_with_distances_to_parents.txt`) would be a table with the samples as rows and the "Distances_to_mother" and "Distances_to_father" as two separate columns with the evolutionary distances values from each infant sample to the respective mother or father.





## Authors




- Michael A. Vig Merino - Main author - TUM - michael.vig@tum.de; michael.vig.merino@gmail.com - [ORCID](https://orcid.org/0009-0004-3180-8480)
- Svenja Weißenberger - TUM - svenja.weissenberger@tum.de - [ORCID](https://orcid.org/0009-0009-4374-0360)





## Citation





If you use this code, please cite:





**Preprint:**





Merino MV, Weißenberger S, Ermolova A, Dhilly E, Spadazzi R, Levy L, Collado MC, Hertz T, Omer H, Heidrich V, Segata N, Larsen M, Schirmer M, Haller D. Probiotic formula intervention in infants leads to colonization and competitive strain displacement independent of IgA binding. bioRxiv (2026). https://doi.org/10.64898/2026.06.03.729775





