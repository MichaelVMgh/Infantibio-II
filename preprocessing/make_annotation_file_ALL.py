import pandas as pd
import numpy as np
from matplotlib import cm
from matplotlib.colors import rgb2hex

# Configuration
# TREE_FILES = [
    # '../RAxML_bestTree.t__SGB17247.StrainPhlAn4_cFMD_newick.txt',
    # '../RAxML_bestTree.t__SGB17256.StrainPhlAn4_cFMD_newick.txt',
    # '../RAxML_bestTree.t__subsp.infantis.StrainPhlAn4_cFMD_newick.txt',
    # '../RAxML_bestTree.t__subsp.longum.StrainPhlAn4_cFMD_newick.txt'
# ]
TREE_FILES = ["subsp.infantis_MIDPOINT_ALL_ACCURATE.txt", "subsp.longum_MIDPOINT_ALL.txt", "SGB17256_MIDPOINT_ALL.txt", "SGB17247_MIDPOINT_ALL.txt"
]

SPECIES_NAMES = {
    't__SGB17247': 'Bifidobacterium breve (SGB17247)',
    't__SGB17256': 'Bifidobacterium bifidum (SGB17256)',
    't__subsp.infantis': 'Bifidobacterium longum infantis',
    't__subsp.longum': 'Bifidobacterium longum longum'
}

# Color scheme
DATASET_COLORS = {
    'HeppnerN_2024': '#FF8C00',  # Orange
    'CM_MTB': '#87CEEB',          # Light blue
    'cMD3': '#4169E1',             # Blue
    'cFMD': '#daff75', # lime yellow
    'sup_str': '#fa0505', # red
    'metaref_isolate' : '#04ff00' # radioactive green
}

CONTINENT_COLORS = {
    'Europe': '#ff9cc0',
    'Asia': '#00ffb7',
    'Africa': '#657cb8',
    'North America': '#a600ff',
    'South America': '#ff30cb',
    #'Oceania': '#bb00ff',
    #'Antarctica': '#FCBAD3'
}

WESTERNIZED_COLORS = {
    'yes': '#ededed',  # White
    'no': '#000000'    # Black
}

FEEDING_COLORS = {
    'A': '#F9BC05',
    'B': '#9B3350',
    'C': '#71884A',
    'D': '#042B44',
    'BF': '#c0c0c0',
    'infant_at_24mo': '#fff8dc' # cornsilk
}

# FOOD_COLORS = {
    # 'fermented_grains': '#ffe175',
    # 'fermented_fruits_and_vegetables': '#75ff78',
    # 'fermented_tubers_and_roots': '#e875ff',
    # 'alcohol': '#800f1c',
    # 'fermented_beverages': '#0f1580',
    # 'dairy': '#786543',
    # 'probiotics': '#5e5e5e',
    # 'fermented_meat': '#591800'
# }


def read_tree_samples(tree_file):
    """Extract sample names from tree file"""
    with open(tree_file, 'r') as f:
        tree_content = f.read()
    
    # Parse tree to extract sample names
    samples = []
    tree_cleaned = tree_content.replace(',', ' ').replace('(', '').replace(')', '').replace(';', '')
    for item in tree_cleaned.split():
        if ':' in item:
            sample_name = item.split(':')[0]
            if sample_name:
                samples.append(sample_name)
    
    return samples

def get_age_color(age, min_age, max_age):
    """Convert age to color using Greens colormap"""
    if pd.isna(age) or max_age == min_age:
        return '#a8a8a8'  # Grey for missing data
    
    # Normalize age to 0-1 range
    normalized = (age - min_age) / (max_age - min_age)
    
    # Use Greens colormap
    colormap = cm.get_cmap('Greens')
    rgba = colormap(normalized)
    
    return rgb2hex(rgba)

def create_annotation_file(tree_file, metadata_df, output_file):
    """Generate GraPhlAn annotation file"""
    
    # Extract species identifier from filename
    species_key = None
    for key in SPECIES_NAMES.keys():
        if key in tree_file or key.replace('t__', '') in tree_file:
            species_key = key
            break
    
    if not species_key:
        print(f"Warning: Could not identify species for {tree_file}")
        return
    
    # Read samples from tree
    samples_in_tree = read_tree_samples(tree_file)
    print(f"Processing {tree_file}: {len(samples_in_tree)} samples found")
    
    # Filter metadata for samples in tree
    metadata_filtered = metadata_df[metadata_df['sample_id'].isin(samples_in_tree)]
    
    # Calculate age range for gradient
    min_age = metadata_filtered['age_in_days'].min()
    max_age = metadata_filtered['age_in_days'].max()
    
    # Start building annotation text
    text = []
    
    # Header settings
    text.append(f'title\t{SPECIES_NAMES[species_key]}')
    text.append('total_plotted_degrees\t340')
    text.append('start_rotation\t90')
    text.append('class_legend_font_size\t10')
    text.append('annotation_background_separation\t0.03')
    text.append('annotation_background_width\t0.0')
    text.append('annotation_font_size\t11')
    text.append('clade_marker_size\t0')
    text.append('branch_bracket_width\t0.0')
    text.append('branch_thickness\t0.3')
    text.append('clade_marker_edge_width\t0.3')
    
    # Ring labels
    text.append('ring_label\t1\tDataset')
    text.append('ring_label\t2\tAge (days)')
    text.append('ring_label\t3\tWesternized')
    text.append('ring_label\t4\tContinent')
    text.append('ring_label\t5\tFeeding')
    #text.append('ring_label\t6\tFood metagenome')
    
    text.append('ring_label_font_size\t1\t8')
    text.append('ring_label_font_size\t2\t8')
    text.append('ring_label_font_size\t3\t8')
    text.append('ring_label_font_size\t4\t8')
    text.append('ring_label_font_size\t5\t8')
    #text.append('ring_label_font_size\t6\t8')
    
    # Ring separators
    for ring in range(1, 6):
        text.append(f'ring_internal_separator_thickness\t{ring}\t0.5')
        text.append(f'ring_external_separator_thickness\t{ring}\t0.5')
        text.append(f'ring_separator_color\t{ring}\tk')
    
    # Legend markers for datasets
    for dataset, color in DATASET_COLORS.items():
        if color == '#fa0505' or color == '#daff75' or color == '#04ff00':
            text.append(f'{dataset}\tclade_marker_shape\t*')
            text.append(f'{dataset}\tclade_marker_size\t100')
            text.append(f'{dataset}\tclade_marker_color\t{color}')
        else:
            text.append(f'{dataset}\tclade_marker_shape\ts')
            text.append(f'{dataset}\tclade_marker_size\t100')
            text.append(f'{dataset}\tclade_marker_color\t{color}')
    
    for dataset, color in CONTINENT_COLORS.items():
        text.append(f'{dataset}\tclade_marker_shape\ts')
        text.append(f'{dataset}\tclade_marker_size\t100')
        text.append(f'{dataset}\tclade_marker_color\t{color}')
    
    for dataset, color in WESTERNIZED_COLORS.items():
        text.append(f'{dataset}\tclade_marker_shape\ts')
        text.append(f'{dataset}\tclade_marker_size\t100')
        text.append(f'{dataset}\tclade_marker_color\t{color}')
    
    for dataset, color in FEEDING_COLORS.items():
        text.append(f'{dataset}\tclade_marker_shape\ts')
        text.append(f'{dataset}\tclade_marker_size\t100')
        text.append(f'{dataset}\tclade_marker_color\t{color}')
    
    # for dataset, color in FOOD_COLORS.items():
        # text.append(f'{dataset}\tclade_marker_shape\ts')
        # text.append(f'{dataset}\tclade_marker_size\t100')
        # text.append(f'{dataset}\tclade_marker_color\t{color}')
    # text.append(f'Present\tclade_marker_shape\ts')
    # text.append(f'Present\tclade_marker_size\t100')
    # text.append(f'Present\tclade_marker_color\t#969400')
    
    # text.append(f'Absent\tclade_marker_shape\ts')
    # text.append(f'Absent\tclade_marker_size\t100')
    # text.append(f'Absent\tclade_marker_color\t#FFFFFF')
        
    # Process each sample
    for sample in samples_in_tree:
        # Get sample metadata
        sample_data = metadata_filtered[metadata_filtered['sample_id'] == sample]
        
        if sample_data.empty:
            print(f"Warning: No metadata for sample {sample}")
            continue
        
        sample_data = sample_data.iloc[0]
        
        # Determine dataset group
        if sample_data['dataset'] == 'HeppnerN_2024':
            dataset_group = 'HeppnerN_2024'
        elif sample_data['dataset'] == 'CM_MTB':
            dataset_group = 'CM_MTB'
        elif sample_data['dataset'] == 'sup_str':
            dataset_group = 'sup_str'
        elif sample_data['dataset'] == 'metaref_isolate':
            dataset_group = 'metaref_isolate'
        elif sample_data['dataset'] == 'cFMD':
            dataset_group = 'cFMD'
        else:
            dataset_group = 'cMD3'
        
        # Get colors
        dataset_color = DATASET_COLORS[dataset_group]
        age_color = get_age_color(sample_data['age_in_days'], min_age, max_age)
        westernized_color = WESTERNIZED_COLORS.get(sample_data['is_westernized'], '#FFFFFF')
        continent_color = CONTINENT_COLORS.get(sample_data['continent'], '#CCCCCC')
        feeding_color = FEEDING_COLORS.get(sample_data['HistoryAtTimePoint'], '#FFFFFF')
        #food_color = FOOD_COLORS.get(sample_data['category'], '#FFFFFF')
        
        # Clade marker
        if dataset_color == '#fa0505':
            text.append(f'{sample}\tclade_marker_color\t{dataset_color}')
            text.append(f'{sample}\tclade_marker_size\t80')
            text.append(f'{sample}\tclade_marker_shape\t*')
        elif dataset_color == '#04ff00':
            text.append(f'{sample}\tclade_marker_color\t{dataset_color}')
            text.append(f'{sample}\tclade_marker_size\t40')
            text.append(f'{sample}\tclade_marker_shape\t*')
        elif dataset_color == '#daff75':
            text.append(f'{sample}\tclade_marker_color\t{dataset_color}')
            text.append(f'{sample}\tclade_marker_size\t30')
            text.append(f'{sample}\tclade_marker_shape\tp')
        else:
            text.append(f'{sample}\tclade_marker_color\t{dataset_color}')
            text.append(f'{sample}\tclade_marker_size\t9')
            text.append(f'{sample}\tclade_marker_shape\to')
        
        # Ring 1: Dataset
        text.append(f'{sample}\tring_color\t1\t{dataset_color}')
        text.append(f'{sample}\tring_width\t1\t1')
        text.append(f'{sample}\tring_height\t1\t0.75')
        
        # Ring 2: Age
        text.append(f'{sample}\tring_color\t2\t{age_color}')
        text.append(f'{sample}\tring_width\t2\t1')
        text.append(f'{sample}\tring_height\t2\t0.75')
        
        # Ring 3: Westernized
        text.append(f'{sample}\tring_color\t3\t{westernized_color}')
        text.append(f'{sample}\tring_width\t3\t1')
        text.append(f'{sample}\tring_height\t3\t0.75')
        
        # Ring 4: Continent
        text.append(f'{sample}\tring_color\t4\t{continent_color}')
        text.append(f'{sample}\tring_width\t4\t1')
        text.append(f'{sample}\tring_height\t4\t0.75')
        
        # RING 5: Feeding_Overall
        text.append(f'{sample}\tring_color\t5\t{feeding_color}')
        text.append(f'{sample}\tring_width\t5\t1')
        text.append(f'{sample}\tring_height\t5\t0.75')
        
        # RING 6: Food category
        # if food_color == "#FFFFFF":
            # text.append(f'{sample}\tring_color\t6\t{food_color}')
        # else:
            # text.append(f'{sample}\tring_color\t6\t#969400')
        # text.append(f'{sample}\tring_width\t6\t1')
        # text.append(f'{sample}\tring_height\t6\t0.75')
    
    # Write to file
    with open(output_file, 'w') as f:
        f.write('\n'.join(text))
    
    print(f"Created annotation file: {output_file}")

def main():
    # Load metadata
    metadata = pd.read_csv('../../data/add_cFMD/total_metadata_merged_ALL.csv')
    print(f"Loaded metadata: {len(metadata)} samples")
    
    # Generate annotation file for each tree
    for tree_file in TREE_FILES:
        # Extract base name for output
        base_name = tree_file.split('MIDPOINT_ALL.txt')[0] + '_cFMD_annotation.txt'
        output_file = f'{base_name}'
        
        create_annotation_file(tree_file, metadata, output_file)
    
    print("\nAll annotation files generated successfully!")

if __name__ == '__main__':
    main()