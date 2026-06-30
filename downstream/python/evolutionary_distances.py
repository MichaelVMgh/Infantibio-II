# -*- coding: utf-8 -*-
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from statannotations.Annotator import Annotator
from scipy.stats import kruskal


# Preprocessing for subplot 1 & 2 (mothers and fathers)
df_mother = pd.read_csv("Buniformis_full_symmetric_distance_matrix.txt", sep="\t")
df_father = df_mother.copy()

def preprocess_relatedness_df(df, role):
    df = df.reset_index().melt(id_vars="sampleid", var_name="Comparison", value_name="Distance")
    df = df[df["sampleid"] != df["Comparison"]]
    df = df[df["Comparison"] != "index"]

    df["Relationship"] = df.apply(
        lambda row: "Related" if row["sampleid"].split("_")[0] == row["Comparison"].split("_")[0] else "Unrelated",
        axis=1,
    )

    df = df[df["sampleid"].str.contains(role) | df["Comparison"].str.contains(role)]
    other_role = "father" if role == "mother" else "mother"
    df = df[~df["sampleid"].str.contains(other_role) & ~df["Comparison"].str.contains(other_role)]

    df["Pair"] = df.apply(lambda row: tuple(sorted([row["sampleid"], row["Comparison"]])), axis=1)
    df = df.drop_duplicates(subset="Pair").drop(columns="Pair")

    df = df[~(df["sampleid"].str.contains(role) & df["Comparison"].str.contains(role))]

    df["TimePoint"] = df.apply(
        lambda row: row["sampleid"].split("_")[-1]
        if row["sampleid"].split("_")[-1] != role
        else row["Comparison"].split("_")[-1],
        axis=1,
    )
    df["TimePoint"] = pd.Categorical(df["TimePoint"], categories=["1mo", "3mo", "7mo", "12mo", "24mo"], ordered=True)
    return df

df1 = preprocess_relatedness_df(df_mother, role="mother")  # Subplot 1
df_f = preprocess_relatedness_df(df_father, role="father")  # Subplot 2


# Preprocessing for subplot 4 & 5
df2 = pd.read_csv("df_buniformis_with_distances_to_parents.txt", sep="\t")
df2 = df2[~df2["TimePoint"].isin(["mother", "father", "12mo", "24mo"])]
df2["Feeding"] = df2["HistoryAtTimePoint"].apply(lambda x: "BF" if x == "BF" else "Formula")
df2["Feeding"] = pd.Categorical(df2["Feeding"], categories=["BF", "Formula"], ordered=True)


#We create 5 vertical subplots
fig, axes = plt.subplots(5, 1, figsize=(12, 45), gridspec_kw={"hspace": 0.5})

#We define the colors
npg_colors = ["#E64B35", "#4DBBD5"]
relationship_order = ["Related", "Unrelated"]

# Subplot 1 (Mothers: related vs unrelated)
sns.boxplot(
    ax=axes[0],
    x="Relationship",
    y="Distance",
    data=df1,
    order=relationship_order,
    showfliers=False,
    width=0.6,
    palette=npg_colors,
    zorder=1,
)

#Generate Scatter: different alpha per group
for group, alpha in zip(["Related", "Unrelated"], [0.50, 0.15]):
    sns.stripplot(
        ax=axes[0],
        data=df1[df1["Relationship"] == group],
        x="Relationship",
        y="Distance",
        order=relationship_order,
        color="black",
        size=3,
        jitter=True,
        alpha=alpha,
        zorder=2,
    )

#Stats
related = df1[df1["Relationship"] == "Related"]["Distance"]
unrelated = df1[df1["Relationship"] == "Unrelated"]["Distance"]
kw_stat, kw_pval = kruskal(related, unrelated)

xlim = axes[0].get_xlim()
ylim = axes[0].get_ylim()
axes[0].text(
    xlim[1] - 0.05,
    ylim[1] * 0.95,
    f"Kruskal-Wallis p = {kw_pval:.3e}",
    ha="right",
    va="top",
    fontsize=20,
    bbox=dict(facecolor="white", edgecolor="none", boxstyle="round,pad=0.3"),
)

axes[0].set_xlabel("Relationship", fontsize=24)
axes[0].set_ylabel("Evolutionary Distance", fontsize=24)
axes[0].tick_params(labelsize=20)
axes[0].set_title("Distances: Related vs Unrelated (mothers)", fontsize=26)


# Subplot 2 (fathers: related vs unrelated)
sns.boxplot(
    ax=axes[1],
    x="Relationship",
    y="Distance",
    data=df_f,
    order=relationship_order,
    showfliers=False,
    width=0.6,
    palette=npg_colors,
    zorder=1,
)

# Scatter: different alpha per group
for group, alpha in zip(["Related", "Unrelated"], [0.50, 0.15]):
    sns.stripplot(
        ax=axes[1],
        data=df_f[df_f["Relationship"] == group],
        x="Relationship",
        y="Distance",
        order=relationship_order,
        color="black",
        size=3,
        jitter=True,
        alpha=alpha,
        zorder=2,
    )

# Stats
related_f = df_f[df_f["Relationship"] == "Related"]["Distance"]
unrelated_f = df_f[df_f["Relationship"] == "Unrelated"]["Distance"]
kw_stat_f, kw_pval_f = kruskal(related_f, unrelated_f)

xlim = axes[1].get_xlim()
ylim = axes[1].get_ylim()
axes[1].text(
    xlim[1] - 0.05,
    ylim[1] * 0.95,
    f"Kruskal-Wallis p = {kw_pval_f:.3e}",
    ha="right",
    va="top",
    fontsize=20,
    bbox=dict(facecolor="white", edgecolor="none", boxstyle="round,pad=0.3"),
)

axes[1].set_xlabel("Relationship", fontsize=24)
axes[1].set_ylabel("Evolutionary Distance", fontsize=24)
axes[1].tick_params(labelsize=20)
axes[1].set_title("Distances: Related vs Unrelated (fathers)", fontsize=26, pad=20)

# Subplot 3 (time points)
df1_related = df1[df1["Relationship"] == "Related"].copy()
#df1_related["TimePoint"] = df1_related["TimePoint"].str.replace("mo", "")
timepoint_order = ["1mo", "3mo", "7mo", "12mo", "24mo"]
palette_timepoints = {
    "1mo": "lightgrey",
    "3mo": "teal",
    "7mo": "aquamarine",
    "12mo": "darkviolet",
    "24mo": "cornsilk",
}

sns.boxplot(
    ax=axes[2],
    data=df1_related,
    x="TimePoint",
    y="Distance",
    order=timepoint_order,
    palette=palette_timepoints,
    showfliers=False,
    width=0.6,
)

sns.stripplot(
    ax=axes[2],
    data=df1_related,
    x="TimePoint",
    y="Distance",
    order=timepoint_order,
    color="black",
    size=8,
    jitter=True,
    alpha=0.6,
)

pairs = [("1mo", "24mo")]
annot5 = Annotator(axes[2], pairs, data=df1_related, x="TimePoint", y="Distance", order=timepoint_order)
annot5.configure(
    test="Mann-Whitney",
    text_format="simple",
    loc="outside",
    verbose=0,
    fontsize=20,
    comparisons_correction="BH",
    line_offset_to_group=0.1,
    line_height=0.03,
)
annot5.apply_and_annotate()

axes[2].set_title("Distances to paired mothers", fontsize=26, pad=70)
axes[2].set_xlabel("Time point (months)", fontsize=24)
axes[2].set_ylabel("Evolutionary Distance", fontsize=24)
axes[2].tick_params(labelsize=20)


# Subplot 4 (Distances to Mothers by feeding)
sns.boxplot(
    ax=axes[3],
    data=df2,
    x="Feeding",
    y="Distances_to_mother",
    palette={"BF": "#c0c0c0", "Formula": "darkturquoise"},
    linewidth=2,
    showfliers=False,
    width=0.6,
)

sns.stripplot(
    ax=axes[3],
    data=df2,
    x="Feeding",
    y="Distances_to_mother",
    jitter=True,
    alpha=0.6,
    marker="o",
    size=10,
    color="black",
)

annot2 = Annotator(axes[3], [("BF", "Formula")], data=df2, x="Feeding", y="Distances_to_mother")
annot2.configure(
    test="Mann-Whitney",
    text_format="simple",
    loc="outside",
    line_offset_to_group=0.05,
    line_height=0.03,
    verbose=0,
    fontsize=26,
    comparisons_correction="BH",
)
axes[3].set_ylim(0, df2["Distances_to_mother"].max() * 1.3)
annot2.apply_and_annotate()

axes[3].set_xlabel("Feeding", fontsize=24)
axes[3].set_ylabel("Evolutionary Distance", fontsize=24)
axes[3].tick_params(labelsize=20)
axes[3].set_title("Distance to paired mothers", fontsize=26)


# Subplot 5 (Distances to fathers by feeding)
sns.boxplot(
    ax=axes[4],
    data=df2,
    x="Feeding",
    y="Distances_to_father",
    palette={"BF": "#c0c0c0", "Formula": "darkturquoise"},
    linewidth=2,
    showfliers=False,
    width=0.6,
)

sns.stripplot(
    ax=axes[4],
    data=df2,
    x="Feeding",
    y="Distances_to_father",
    jitter=True,
    alpha=0.6,
    marker="o",
    size=8,
    color="black",
)

ymax = df2["Distances_to_father"].max()
axes[4].set_ylim(0, ymax * 1.3)

annot3 = Annotator(axes[4], [("BF", "Formula")], data=df2, x="Feeding", y="Distances_to_father")
annot3.configure(
    test="Mann-Whitney",
    text_format="simple",
    loc="inside",
    line_offset_to_group=0.05,
    line_height=0.03,
    verbose=0,
    fontsize=26,
    comparisons_correction="BH",
)
annot3.apply_and_annotate()

axes[4].set_xlabel("Feeding", fontsize=24)
axes[4].set_ylabel("Evolutionary Distance", fontsize=24)
axes[4].tick_params(labelsize=20)
axes[4].set_title("Distance to paired fathers", fontsize=26)

# Save
sns.despine()
pos3 = axes[2].get_position()
axes[2].set_position([pos3.x0, pos3.y0 - 0.01, pos3.width, pos3.height])
plt.subplots_adjust(top=0.92)
#plt.savefig("Buniformis_combined_5plots_30.05.2026.pdf", dpi=800, bbox_inches="tight")
plt.show()


