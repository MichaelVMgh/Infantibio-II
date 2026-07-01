# -*- coding: utf-8 -*-
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from scipy.spatial.distance import pdist, squareform


# Load the data
df = pd.read_csv(
    "merged_abundance_table_species_transposed.txt",
    sep="\t"
)

# Keep only the 24-month samples
df = df[df["TimePoint"] == "24mo"]


# Split metadata and taxa
metadata = df.iloc[:, :12]
taxa = df.iloc[:, 12:]


# Aggregate by Feeding_Overall (overall feeding column)
taxa_grouped = (
    taxa
    .assign(Feeding_Overall=metadata["Feeding_Overall"])
    .groupby("Feeding_Overall")
    .median()
)


# We create Bray-Curtis distance matrix
bray_dist = pdist(
    taxa_grouped.values,
    metric="braycurtis"
)

bray_matrix = pd.DataFrame(
    squareform(bray_dist),
    index=taxa_grouped.index,
    columns=taxa_grouped.index
)


# We sort the feeding groups
desired_order = ["A", "B", "C", "D", "BF"]

#We keep only the groups present
desired_order = [
    x for x in desired_order
    if x in bray_matrix.index
]

#We define the Bray Curtis matrix
bray_matrix = bray_matrix.loc[
    desired_order,
    desired_order
]

# Feeding group colors
feeding_palette = {
    "A": "#F9BC05",
    "B": "#9B3350",
    "C": "#71884A",
    "D": "#042B44",
    "BF": "#c0c0c0"
}

row_colors = [feeding_palette[x] for x in desired_order]
col_colors = [feeding_palette[x] for x in desired_order]

#Create heatmap
g = sns.clustermap(
    bray_matrix,
    cmap="Blues",
    row_colors=row_colors,
    col_colors=col_colors,
    row_cluster=False,
    col_cluster=False,
    linewidths=0.5,
    figsize=(10, 8),
    vmin=0,
    vmax=1,
    cbar_kws={"label": "Bray–Curtis dissimilarity"}
)


# Labels
g.ax_heatmap.set_xlabel(
    "Overall Feeding",
    fontsize=22
    #fontweight="bold"
)

#Remove y-axis label
g.ax_heatmap.set_ylabel("")

# Generate the tick labels
g.ax_heatmap.tick_params(
    axis="x",
    labelsize=18,
    rotation=0
)

g.ax_heatmap.tick_params(
    axis="y",
    labelsize=18
)

g.fig.subplots_adjust(right=0.80)

# Move colorbar
g.cax.set_position([
    0.86,   # x
    0.25,   # y
    0.025,  # width
    0.50    # height
])

g.cax.set_ylabel(
    "Bray–Curtis dissimilarity",
    fontsize=16
    #fontweight="bold"
)

g.cax.tick_params(
    labelsize=14
)

g.cax.yaxis.set_label_position("right")
g.cax.yaxis.tick_right()

# Save
#g.fig.savefig(
#    "bray_curtis_overall_feeding_24mo_30.05.2026.png",
#    dpi=800,
#    bbox_inches="tight"
#)
# PDF format
#g.fig.savefig(
#    "bray_curtis_overall_feeding_24mo_30.05.2026.pdf",
#    bbox_inches="tight"
#)

plt.show()

