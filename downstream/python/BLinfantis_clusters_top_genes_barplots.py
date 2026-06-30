# -*- coding: utf-8 -*-
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np

# Load significant genes: output from Maaslin, where "clade" refers to the "cluster" variable (1 or 2)
df_significant_genes = pd.read_csv("significant_results_clade.txt", sep="\t")
df_significant_genes.rename(columns={"feature": "protein_id"}, inplace=True)
df_significant_genes["protein_id_clean"] = df_significant_genes["protein_id"].str.replace(".", "_", regex=False)

# Load Pfam annotations
df01 = pd.read_csv("FG_interproscan_Pfam.tsv", sep="\t")
df_pfams = df01.copy()
df_pfams["protein_id_clean"] = df_pfams["protein_id"].str.replace(":", "_", regex=False) \
                                                     .str.replace("+", "_", regex=False) \
                                                     .str.replace("-", "_", regex=False)

# Merge on the normalized column
df_significant_genes_with_pfams = df_significant_genes.merge(
    df_pfams, 
    on="protein_id_clean", 
    how="left"
)

# Save to file
#df_significant_genes_with_pfams.to_csv("significant_results_clade_with_pfams.txt", sep="\t", index=False)


# We define a score: high coef + low qval
df_significant_genes_with_pfams["score"] = np.sign(df_significant_genes_with_pfams["coef"]) * -np.log10(df_significant_genes_with_pfams["qval"])
df_ranked = df_significant_genes_with_pfams.sort_values("score", ascending=False)
df_ranked.to_csv("significant_results_clade_with_pfams_and_scores.txt", sep="\t", index=False) #we add the scores to the significant genes file (pfams also added)
#df_significant_genes_with_pfams.to_csv("significant_results_clade_with_pfams.txt", sep="\t", index=False)

# We can load the updated file or directly use df_ranked
df_ranked = pd.read_csv("significant_results_clade_with_pfams_and_scores_updated.txt", sep="\t")


df_annotated = df_ranked.dropna(subset=["FG"]) # Drop rows with missing FG

top_fg_max_scores = ( # Top = max score per FG
    df_annotated.groupby("FG")["score"].max()
    .sort_values(ascending=False)
    .head(16)
)

bottom_fg_min_scores = ( # Bottom = min score per FG
    df_annotated.groupby("FG")["score"].min()
    .sort_values(ascending=True)
    .head(17)
)


extreme_fg_scores = pd.concat([top_fg_max_scores, bottom_fg_min_scores]) # Combine top & bottom
extreme_fg_scores = extreme_fg_scores.sort_values(ascending=False)


# Annotation mapping
fg_annotations = {
    "PF05713": "Bacterial mobilisation protein (MobC)",
    "PF03432": "MobA/VirD2-like, nuclease domain",
    "PF07275": "Antirestriction protein (ArdA)",
    "PF02195": "ParB N-terminal domain",
    "PF03837": "RecT family",
    "PF00929": "Exonuclease",
    "PF08240": "Alcohol dehydrogenase GroES-like domain",
    "PF00188": "Cysteine-rich secretory protein family",
    "PF03412:::PF00664:::PF00005": "Peptidase C39; ABC TM; ABC transporter",
    "PF03412:::PF00664": "Peptidase C39; ABC TM",
    "PF13575:::PF05147": "Lanthionine synthetase C-like protein",
    "PF00933:::PF01915:::PF14310":"Glycosyl hydrolase family 3C/3N; Fibronectin type III-like",
    "PF02837:::PF00703:::PF02836":"Glycosyl hydrolases family 2, sugar binding and TIM barrel domains",
    "PF02055:::PF17189":"Glycosyl hydrolase family 30 TIM-barrel and beta sandwich domains",
    "PF00528":"Binding-protein-dependent transport system inner membrane component",
    "PF00440:::PF13977":"Bacterial regulatory proteins, tetR family; BetI-type transcriptional repressor",
    "PF14867":"Lantibiotic alpha",
    "PF13276:::PF00665":"HTH-like domain; Integrase core domain",
    "PF17998:::PF16364":"Cell surface antigen C terminal domains",
    "PF13518":"Helix-turn-helix domain",
    "PF13730":"Helix-turn-helix domain",
    "PF00041":"Fibronectin type III domain",
    "PF13333":"Integrase core domain",
    "PF00665:::PF13333":"Integrase core domain",
    "PF13518:::PF13276:::PF00665":"HTH-like domain; Integrase core domain",
    "PF00358":"phosphoenolpyruvate-dependent sugar phosphotransferase system, EIIA 1",
    "PF13276:::PF00665:::PF13333":"HTH-like domain; Integrase core domain",
    "PF21088:::PF00924":"Mechanosensitive ion channel",
    "PF01695":"IstB-like ATP binding protein",
    "PF14691":"Dihydroprymidine dehydrogenase domain II",
    "PF00005:::PF12848:::PF00005":"ABC transporter",
    "PF13749:::PF01978":"ATP-dependent DNA helicase recG; TrmB",
    "PF05977":"Transmembrane secretion effector",
    "PF01142":"tRNA pseudouridine synthase D (TruD)",
    "PF13350:::PF00149":"Tyrosine phosphatase; Calcineurin-like phosphoesterase",
    "PF00528:::PF00528":"Binding-protein-dependent transport system inner membrane component",
    "PF01418:::PF01380":"HTH rpiR family; SIS domain",
    "PF00665":"Integrase core domain",
    "PF02481":"DNA recombination-mediator protein A",
    "PF14690:::PF01610":"zinc-finger of transposase; Transposase",
    "PF20038":"Helix-turn-helix domain",
    "PF13470":"PIN domain",
    "PF01053":"Cys/Met metabolism PLP-dependent enzyme",
    "PF01569:::PF00149:::PF02872":"PAP2; Calcineurin-like phosphoesterase",
    "PF12728":"Helix-turn-helix domain",
    "PF00892:::PF00892":"EamA-like transporter",
    "PF00005":"ABC transporter",
    "PF00072:::PF00486":"Response regulator (RR)",
    "PF02687":"FtsX-like permease",
    "PF00512:::PF02518":"His kinase / HSP90-like ATPase",
    "PF07690":"MFS transporter",
    "PF13443":"Cro/C1 HTH domain",
    "PF00375":"Sodium:dicarboxylate symporter",
    "PF12704:::PF02687":"MacB-like + FtsX-like",
    "PF07751":"Abi-like protein",
    "PF00440":"TetR regulator",
    "PF01610":"Transposase",
    "PF00196":"LuxR family regulator",
    "PF13343":"Extracellular solute-binding protein",
    "PF13657":"HipA N-terminal domain",
    "PF00436":"SSB protein",
    "PF13350:::PF00149":"Tyrosine phosphatase; Calcineurin-like phosphoesterase"
}


# Assign functional Category to each PFAM
fg_to_category = {
    "PF05713": "Mobilization",
    "PF03432": "Mobilization",
    "PF07275": "Mobilization",
    "PF01695": "Mobilization",

    "PF02195": "Recombination",
    "PF03837": "Recombination",
    "PF00929": "Recombination",
    "PF13749:::PF01978": "Recombination",
    "PF02481": "Recombination",
    "PF00436": "Recombination",

    "PF00528": "Transport",
    "PF00528:::PF00528": "Transport",
    "PF03412:::PF00664:::PF00005": "Transport",
    "PF03412:::PF00664": "Transport",
    "PF00005": "Transport",
    "PF00005:::PF12848:::PF00005": "Transport",
    "PF00375": "Transport",
    "PF07690": "Transport",
    "PF02687": "Transport",
    "PF12704:::PF02687": "Transport",
    "PF13343": "Transport",

    "PF00933:::PF01915:::PF14310": "Carbohydrate degradation",
    "PF02837:::PF00703:::PF02836": "Carbohydrate degradation",
    "PF02055:::PF17189": "Carbohydrate degradation",
    "PF00358": "Carbohydrate degradation",

    "PF00188": "Cell surface",
    "PF17998:::PF16364": "Cell surface",
    "PF21088:::PF00924": "Cell surface",
    "PF00041": "Cell surface",

    "PF00440:::PF13977": "Regulation",
    "PF13518": "Regulation",
    "PF13730": "Regulation",
    "PF01418:::PF01380": "Regulation",
    "PF00072:::PF00486": "Regulation",
    "PF00512:::PF02518": "Regulation",
    "PF13443": "Regulation",
    "PF00440": "Regulation",
    "PF00196": "Regulation",
    "PF20038": "Regulation",
    "PF12728": "Regulation",

    "PF13276:::PF00665": "Integrase/Transposase",
    "PF13333": "Integrase/Transposase",
    "PF00665:::PF13333": "Integrase/Transposase",
    "PF13518:::PF13276:::PF00665": "Integrase/Transposase",
    "PF13276:::PF00665:::PF13333": "Integrase/Transposase",
    "PF00665": "Integrase/Transposase",
    "PF14690:::PF01610": "Integrase/Transposase",
    "PF01610": "Integrase/Transposase",

    "PF13575:::PF05147": "Antimicrobial activity",
    "PF14867": "Antimicrobial activity",
    "PF05977": "Antimicrobial activity",
    "PF07751": "Antimicrobial activity",

    "PF08240": "Metabolism",
    "PF14691": "Metabolism",
    "PF01142": "Metabolism",
    "PF13350:::PF00149": "Metabolism",
    "PF01053": "Metabolism",
    "PF01569:::PF00149:::PF02872": "Metabolism",

    "PF13470": "Stress responses",
    "PF13657": "Stress responses"
}


# Now we assign colors for each function
category_colors = {
    "Mobilization": "#e41a1c",
    "Recombination": "#377eb8",
    "Transport": "#4daf4a",
    "Carbohydrate degradation": "#984ea3",
    "Cell surface": "#ff7f00",
    "Regulation": "#a65628",
    "Integrase/Transposase": "#f781bf",
    "Antimicrobial activity": "#999999",
    "Metabolism": "#66c2a5",
    "Stress responses": "#fc8d62",
}


# Color per bar
bar_colors = [
    category_colors.get(fg_to_category.get(fg, "Metabolism"), "#333333")
    for fg in extreme_fg_scores.index
]


# Labels
fg_labels = [fg_annotations.get(fg, fg) for fg in extreme_fg_scores.index]


# Horizontal bar plot
sns.set_style("white")
plt.figure(figsize=(58, 34))  # Wider than tall

sns.barplot(
    y=fg_labels,                 # labels on y-axis
    x=extreme_fg_scores.values,  # values on x-axis
    palette=bar_colors,
    orient="h"
)


# X-axis and y-axis labels
plt.xlabel("Score: sign(log2FC) × -log10(qval)", fontsize=80)
plt.ylabel("Functional Group (FG)", fontsize=80)
plt.title(
    "Top & Bottom BL subsp. infantis genes with annotated FGs (Cluster 1 vs Cluster 2)",
    fontsize=60,
    loc="left",
    pad=60
)

plt.xticks(fontsize=60)
plt.yticks(fontsize=60)

# Legend
handles = [plt.Rectangle((0,0),1,1,color=col) for col in category_colors.values()]

plt.legend(
    handles,
    category_colors.keys(),
    title="Functional category",
    fontsize=50,
    title_fontsize=55,
    loc="lower right"
)

plt.tight_layout()
#Save the figure
plt.savefig(
    "top_bottom_genes_with_pfams_cluster1vs2_COLORED_WITH_LEGEND_HORIZONTAL_LARGE_LABELS_01.12.2025.png",
    dpi=600,
    facecolor="white"
)
plt.show()

