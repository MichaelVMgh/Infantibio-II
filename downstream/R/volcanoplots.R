##################################### new BF vs formula-fed Volcano plot with highlighted Carbohydrate degradation pathways ############################################
#Volcano plot
# Volcano plot
setwd("/path/to/your/data/bf_vs_non-bf/gene_pathways/Feeding-without-24mo_antib_bm_sex")
# Load output from Maaslin with significant genes (gene families or pathways) related to feeding (BF vs Formula // A vs B,C,D)
volcano_df <- read.csv("all_results_feeding.txt", sep="\t")

library(ggplot2)
#install.packages("ggrepel")
library(ggrepel)
#install.packages("ggtext")
library(ggtext)

# Filter significant points
significant_points <- volcano_df[volcano_df$qval < 0.25, ]

# Pathways to label
pathways_to_label <- c(
  "P124.PWY..Bifidobacterium.shunt",
  "GLYCOCAT.PWY..glycogen.degradation.I",
  "PWY.4041...gamma..glutamyl.cycle",
  "PWY0.1061..superpathway.of.L.alanine.biosynthesis",
  "P164.PWY..purine.nucleobases.degradation.I..anaerobic.",
  "PWY.5030..L.histidine.degradation.III",
  "PWY.5505..L.glutamate.and.L.glutamine.biosynthesis",
  "PWY.7456...beta...1.4..mannan.degradation",
  "PWY.4984..urea.cycle",
  "PWY.5189..tetrapyrrole.biosynthesis.II..from.glycine.",
  "HOMOSER.METSYN.PWY..L.methionine.biosynthesis.I",
  "PWY0.1338..polymyxin.resistance",
  "METSYN.PWY..superpathway.of.L.homoserine.and.L.methionine.biosynthesis"
)

# Clean labels
pathway_labels <- c(
  "P124.PWY..Bifidobacterium.shunt" = "italic(Bifidobacterium)~shunt",
  "GLYCOCAT.PWY..glycogen.degradation.I" = '"Glycogen degradation"',
  "PWY.4041...gamma..glutamyl.cycle" = '"Gamma glutamyl cycle"',
  "PWY0.1061..superpathway.of.L.alanine.biosynthesis" = '"Superpathway of L-alanine biosynthesis"',
  "P164.PWY..purine.nucleobases.degradation.I..anaerobic." = '"Purine nucleobases degradation"',
  "PWY.5030..L.histidine.degradation.III" = '"L-histidine degradation"',
  "PWY.5505..L.glutamate.and.L.glutamine.biosynthesis" = '"L-glutamate and L-glutamine biosynthesis"',
  "PWY.7456...beta...1.4..mannan.degradation" = '"Beta-1,4 Mannan degradation"',
  "PWY.4984..urea.cycle" = '"Urea cycle"',
  "PWY.5189..tetrapyrrole.biosynthesis.II..from.glycine." = '"Tetrapyrrole biosynthesis"',
  "HOMOSER.METSYN.PWY..L.methionine.biosynthesis.I" = '"L-methionine biosynthesis"',
  "PWY0.1338..polymyxin.resistance" = '"Polymyxin resistance"',
  "METSYN.PWY..superpathway.of.L.homoserine.and.L.methionine.biosynthesis" = '"Superpathway of L-homoserine and L-methionine biosynthesis"'
)

highlighted_points <- significant_points[significant_points$feature %in% pathways_to_label, ]
highlighted_points$pretty_label <- pathway_labels[highlighted_points$feature]

# Pathway groups
blue_pathways <- c(
  "P124.PWY..Bifidobacterium.shunt",
  "PWY.7851..coenzyme.A.biosynthesis.II..eukaryotic.",
  "COA.PWY..coenzyme.A.biosynthesis.I..prokaryotic.",
  "COA.PWY.1..superpathway.of.coenzyme.A.biosynthesis.III..mammals.",
  "PWY.4041...gamma..glutamyl.cycle",
  "PWY0.1061..superpathway.of.L.alanine.biosynthesis"
)
green_pathways <- c("GLYCOCAT.PWY..glycogen.degradation.I")

volcano_df$color_group <- factor(
  ifelse(volcano_df$feature %in% blue_pathways, "B. longum",
         ifelse(volcano_df$feature %in% green_pathways, "B. breve", "Other")),
  levels = c("B. longum", "B. breve", "Other")
)

# Shapes
triangle_pathways <- c(
  "P124.PWY..Bifidobacterium.shunt",
  "GLYCOCAT.PWY..glycogen.degradation.I"
)

volcano_df$shape_group <- factor(
  ifelse(volcano_df$feature %in% triangle_pathways,
         "Carbohydrate degradation", "Other"),
  levels = c("Carbohydrate degradation", "Other")
)

# Volcano plot
p <- ggplot(volcano_df, aes(x = coef, y = -log10(qval))) +
  geom_hline(yintercept = -log10(0.25), col = "red", linewidth = 1) +
  annotate("text", x = min(volcano_df$coef) * 1.2, y = 0.8, #pvalue label
           label = "adj_pvalue = 0.25", col = "red", size = 18, hjust = 0) +
  
  geom_text_repel(
    data = highlighted_points,
    aes(x = coef, y = -log10(qval), label = pretty_label),
    parse = TRUE,
    force = 30,
    force_pull = 0.4,
    box.padding = 0.5,
    point.padding = 4.0,
    max.overlaps = Inf,
    min.segment.length = 0,
    segment.size = 0.8,
    segment.curvature = 0,
    segment.angle = 50,
    segment.ncp = 1,
    size = 15, #size of pathway labels
    inherit.aes = FALSE,
    max.iter = 5000,
    max.time = 3,
    nudge_y = 1
  ) +
  

geom_point(
  aes(
    fill = color_group,
    size = dominant_median_contribution,
    shape = shape_group
  ),
  color = "black",
  stroke = 0.6,
  alpha = 0.9
) +
  
  labs(
    x = "Log2 Fold Change (breast-fed vs. formula-fed)",
    y = "-Log10 Adjusted P-Value",
    fill = "Dominant species",
    size = "Median CPM abundance"
  ) +
  
  scale_fill_manual(
    values = c(
      "B. longum" = "slateblue",
      "B. breve" = "#5f9e6e",
      "Other" = "black"
    ),
    labels = c(
      "B. longum" = "*B. longum*",
      "B. breve" = "*B. breve*",
      "Other" = "Other"
    )
  ) +
  
  scale_size_continuous(
    range = c(4, 24),
    breaks = c(1, 200, 400, 600)
  ) +
  
 
# Shape in legend
scale_shape_manual(
  values = c(
    "Carbohydrate degradation" = 24,  # filled triangle
    "Other" = 21                      # filled circle
  )
) +
  
#Guides
guides(
  fill = guide_legend(
    override.aes = list(shape = 22, color = "black", size = 10),
    order = 1,
    title.position = "top"
  ),
  size = guide_legend(
    override.aes = list(shape = 21, fill = "grey", color = "black"),
    order = 2,
    title.position = "top"
  ),
  shape = guide_legend(
    title = "Functional category",
    override.aes = list(size = 10, fill = "grey"),
    order = 3,
    title.position = "top"
  )
) +
  
  theme_classic() +
  theme(
    axis.title.x = element_markdown(size = 60, face = "bold"),
    axis.title.y = element_markdown(size = 60, face = "bold"),
    axis.text = element_markdown(size = 52),
    plot.title = element_markdown(size = 70, hjust = 0.5, face = "bold"),
    legend.title = element_markdown(size = 52, face = "bold"),
    legend.text = element_markdown(size = 48),
    legend.position = c(0.05, 1.02),
    legend.justification = c("left", "top"),
    legend.box = "vertical",
    plot.margin = margin(t = 40, r = 20, b = 20, l = 20)
  )

# Save
#ggsave("volcano_plot_gene_pathways_babies_feeding_prev-filtered-no24mo_19.12.2025.png",
#       plot = p, width = 33, height = 23, dpi = 800, limitsize=FALSE)





######################## new Supplemented vs Placebo Volcano plot with Carbohydrate degradation highlighted #################################




# Volcano plot
setwd("/path/to/your/data/supplement_vs_placebo_29.07/gene_families")
volcano_df <- read.csv("all_results_supplementation_updated.txt", sep="\t")

library(ggplot2)
library(ggrepel)
library(ggtext)   # for markdown/HTML text formatting

# Filter significant points
significant_points <- volcano_df[volcano_df$qval < 0.05, ]

# To label (gene families in this case)
pathways_to_label <- c(
  # "RXN.12280...expasy..Isoamylase..3.2.1.68.",
  # "RXN.14380...expasy..Isoamylase..3.2.1.68.",
  # "RXN.4301...expasy..Isoamylase..3.2.1.68.",
  # "RXN.21218...expasy..Bleomycin.hydrolase..3.4.22.40.",
  # "RXN.21219...expasy..Bleomycin.hydrolase..3.4.22.40.",
  # "RXN.21220...expasy..Bleomycin.hydrolase..3.4.22.40.",
  # "RXN.18913...expasy..Alpha.mannosidase..3.2.1.24.",
  # "RXN.15035...expasy..Lysophospholipase..3.1.1.5.",
  # "RXN.16225...expasy..Lysophospholipase..3.1.1.5.",
  # "RXN.19311...expasy..Lysophospholipase..3.1.1.5.",
  # "RXN.19312...expasy..Lysophospholipase..3.1.1.5.",
  # "RXN.21069...expasy..Lysophospholipase..3.1.1.5.",
  # "RXN.14281...expasy..Alpha.glucosidase..3.2.1.20.",
  # "RXN.14282...expasy..Alpha.glucosidase..3.2.1.20.",
  # "RXN.14283...expasy..Alpha.glucosidase..3.2.1.20.",
  # "RXN.15910...expasy..Alpha.glucosidase..3.2.1.20.",
  # "RXN0.5183...expasy..Alpha.glucosidase..3.2.1.20.",
  #"X6.PHOSPHOFRUCTO.2.KINASE.RXN...expasy..6.phosphofructo.2.kinase..2.7.1.105.",
  "X6.PHOSPHOFRUCTO.2.KINASE.RXN...expasy..6.phosphofructo.2.kinase..2.7.1.105.",
  "RXN.9297...expasy..N.acetylhexosamine.1.kinase..2.7.1.162.",
  "X3.6.4.6.RXN...expasy..Vesicle.fusing.ATPase..3.6.4.6.",
  "X3.2.1.68.RXN...expasy..Isoamylase..3.2.1.68.",
  "OLEATE.HYDRATASE.RXN...expasy..Oleate.hydratase..4.2.1.53.",
  "RXN.13685...expasy..Adrenodoxin.NADP....reductase..1.18.1.6.",
  "X3.5.4.30.RXN...expasy..dCTP.deaminase..dUMP.forming...3.5.4.30.",
  "X3.4.22.40.RXN...expasy..Bleomycin.hydrolase..3.4.22.40.",
  "X3.1.3.46.RXN...expasy..Fructose.2.6.bisphosphate.2.phosphatase..3.1.3.46.",
  "RXN.16707...expasy..Prokaryotic.ubiquitin.like.protein.ligase..6.3.1.19.",
  "X3.2.1.24.RXN...expasy..Alpha.mannosidase..3.2.1.24.",
  "RXN.12466...expasy..tRNA..adenine.58..N.1...methyltransferase..2.1.1.220.",
  "RXN.13028...expasy..Starch.synthase..maltosyl.transferring...2.4.99.16.",
  "ACP.S.ACETYLTRANSFER.RXN...expasy...Acyl.carrier.protein..S.acetyltransferase..2.3.1.38.",
  "ACYLPHOSPHATASE.RXN...expasy..Acylphosphatase..3.6.1.7.",
  "X3.4.21.92.RXN...expasy..Endopeptidase.Clp..3.4.21.92.",
  "X6.3.5.7.RXN...expasy..Glutaminyl.tRNA.synthase..glutamine.hydrolyzing...6.3.5.7.",
  "LYSOPHOSPHOLIPASE.RXN...expasy..Lysophospholipase..3.1.1.5.",
  "TDCEACT.RXN...expasy...Formate.C.acetyltransferase..activating.enzyme..1.97.1.4.",
  "MALTODEXGLUCOSID.RXN...expasy..Alpha.glucosidase..3.2.1.20.",
  "GLUCOSE.6.PHOSPHATE.1.EPIMERASE.RXN...expasy..Glucose.6.phosphate.1.epimerase..5.1.3.15.",
  "L.FUCONATE.HYDRATASE.RXN...expasy..L.fuconate.dehydratase..4.2.1.68.",
  "ATP.ADENYLYLTRANSFERASE.RXN...expasy..ATP.adenylyltransferase..2.7.7.53.",
  "D.XYLULOSE.REDUCTASE.RXN...expasy..D.xylulose.reductase..1.1.1.9.",
  "RXN.12093...expasy..7.cyano.7.deazaguanine.synthase..6.3.4.20."
)

# Clean names
pathway_labels <- c(
  "X6.PHOSPHOFRUCTO.2.KINASE.RXN...expasy..6.phosphofructo.2.kinase..2.7.1.105." = "6-phosphofructo 2-kinase",
  "RXN.9297...expasy..N.acetylhexosamine.1.kinase..2.7.1.162." = "N-acetylhexosamine 1-kinase",
  "X3.6.4.6.RXN...expasy..Vesicle.fusing.ATPase..3.6.4.6." = "Vesicle fusing ATPase",
  "X3.2.1.68.RXN...expasy..Isoamylase..3.2.1.68." = "Isoamylase",
  "OLEATE.HYDRATASE.RXN...expasy..Oleate hydratase..4.2.1.53." = "Oleate hydratase",
  "RXN.13685...expasy..Adrenodoxin.NADP....reductase..1.18.1.6." = "Adrenodoxin NADP reductase",
  "X3.5.4.30.RXN...expasy..dCTP.deaminase..dUMP.forming...3.5.4.30." = "dCTP deaminase",
  "X3.4.22.40.RXN...expasy..Bleomycin.hydrolase..3.4.22.40." = "Bleomycin hydrolase",
  "X3.1.3.46.RXN...expasy..Fructose.2.6.bisphosphate.2.phosphatase..3.1.3.46." = "Fructose 2,6-biphosphate 2-phosphatase",
  "RXN.16707...expasy..Prokaryotic.ubiquitin.like.protein.ligase..6.3.1.19." = "Ubiquitin ligase",
  "X3.2.1.24.RXN...expasy..Alpha.mannosidase..3.2.1.24." = "Alpha mannosidase",
  "RXN.12466...expasy..tRNA..adenine.58..N.1...methyltransferase..2.1.1.220." = "Methyltransferase",
  "RXN.13028...expasy..Starch.synthase..maltosyl.transferring...2.4.99.16." = "Starch synthase",
  "ACP.S.ACETYLTRANSFER.RXN...expasy...Acyl.carrier.protein..S.acetyltransferase..2.3.1.38." = "Acetyltransferase",
  "ACYLPHOSPHATASE.RXN...expasy..Acylphosphatase..3.6.1.7." = "Acyl Phosphatase",
  "X3.4.21.92.RXN...expasy..Endopeptidase.Clp..3.4.21.92." = "Endopeptidase",
  "X6.3.5.7.RXN...expasy..Glutaminyl.tRNA.synthase..glutamine.hydrolyzing...6.3.5.7." = "Glutaminyl tRNA synthase",
  "LYSOPHOSPHOLIPASE.RXN...expasy..Lysophospholipase..3.1.1.5." = "Lysophospholipase",
  "TDCEACT.RXN...expasy...Formate.C.acetyltransferase..activating.enzyme..1.97.1.4." = "Formate C-acetyltransferase",
  "MALTODEXGLUCOSID.RXN...expasy..Alpha.glucosidase..3.2.1.20." = "Alpha glucosidase",
  "GLUCOSE.6.PHOSPHATE.1.EPIMERASE.RXN...expasy..Glucose.6.phosphate.1.epimerase..5.1.3.15." = "Glucose 6-phosphate 1-epimerase",
  "L.FUCONATE.HYDRATASE.RXN...expasy..L.fuconate.dehydratase..4.2.1.68." = "L-fuconate hydratase",
  "ATP.ADENYLYLTRANSFERASE.RXN...expasy..ATP.adenylyltransferase..2.7.7.53." = "ATP Adenylyl transferase",
  "D.XYLULOSE.REDUCTASE.RXN...expasy..D.xylulose.reductase..1.1.1.9." = "D-Xylulose reductase",
  "RXN.12093...expasy..7.cyano.7.deazaguanine.synthase..6.3.4.20." = "PreQ0 synthase"
)

highlighted_points <- significant_points[significant_points$feature %in% pathways_to_label, ]
highlighted_points$pretty_label <- pathway_labels[highlighted_points$feature]

# Families for blue
blue_pathways <- c(
  "X6.PHOSPHOFRUCTO.2.KINASE.RXN...expasy..6.phosphofructo.2.kinase..2.7.1.105.",
  "RXN.9297...expasy..N.acetylhexosamine.1.kinase..2.7.1.162.",
  "X3.6.4.6.RXN...expasy..Vesicle.fusing.ATPase..3.6.4.6.",
  "X3.2.1.68.RXN...expasy..Isoamylase..3.2.1.68.",
  "RXN.12280...expasy..Isoamylase..3.2.1.68.",
  "RXN.14380...expasy..Isoamylase..3.2.1.68.",
  "RXN.4301...expasy..Isoamylase..3.2.1.68.",
  "OLEATE.HYDRATASE.RXN...expasy..Oleate.hydratase..4.2.1.53.",
  "RXN.13685...expasy..Adrenodoxin.NADP....reductase..1.18.1.6.",
  "X3.5.4.30.RXN...expasy..dCTP.deaminase..dUMP.forming...3.5.4.30.",
  "X3.4.22.40.RXN...expasy..Bleomycin.hydrolase..3.4.22.40.",
  "RXN.21218...expasy..Bleomycin.hydrolase..3.4.22.40.",
  "RXN.21219...expasy..Bleomycin.hydrolase..3.4.22.40.",
  "RXN.21220...expasy..Bleomycin.hydrolase..3.4.22.40.",
  "X3.1.3.46.RXN...expasy..Fructose.2.6.bisphosphate.2.phosphatase..3.1.3.46.",
  "RXN.16707...expasy..Prokaryotic.ubiquitin.like.protein.ligase..6.3.1.19.",
  "X3.2.1.24.RXN...expasy..Alpha.mannosidase..3.2.1.24.",
  "RXN.18913...expasy..Alpha.mannosidase..3.2.1.24.",
  "RXN.12466...expasy..tRNA..adenine.58..N.1...methyltransferase..2.1.1.220.",
  "RXN.13028...expasy..Starch.synthase..maltosyl.transferring...2.4.99.16.",
  "ACP.S.ACETYLTRANSFER.RXN...expasy...Acyl.carrier.protein..S.acetyltransferase..2.3.1.38.",
  "ACYLPHOSPHATASE.RXN...expasy..Acylphosphatase..3.6.1.7.",
  "X3.4.21.92.RXN...expasy..Endopeptidase.Clp..3.4.21.92.",
  "X6.3.5.7.RXN...expasy..Glutaminyl.tRNA.synthase..glutamine.hydrolyzing...6.3.5.7.",
  "LYSOPHOSPHOLIPASE.RXN...expasy..Lysophospholipase..3.1.1.5.",
  "RXN.15035...expasy..Lysophospholipase..3.1.1.5.",
  "RXN.16225...expasy..Lysophospholipase..3.1.1.5.",
  "RXN.19311...expasy..Lysophospholipase..3.1.1.5.",
  "RXN.19312...expasy..Lysophospholipase..3.1.1.5.",
  "RXN.21069...expasy..Lysophospholipase..3.1.1.5.",
  "TDCEACT.RXN...expasy...Formate.C.acetyltransferase..activating.enzyme..1.97.1.4.",
  "RXN.12093...expasy..7.cyano.7.deazaguanine.synthase..6.3.4.20.",
  "X6.PHOSPHOFRUCTO.2.KINASE.RXN...expasy..6.phosphofructo.2.kinase..2.7.1.105.",
  "L.FUCONATE.HYDRATASE.RXN...expasy..L.fuconate.dehydratase..4.2.1.68."
)
#Families for green
green_pathways <- c(
  "MALTODEXGLUCOSID.RXN...expasy..Alpha.glucosidase..3.2.1.20.",
  "RXN.14281...expasy..Alpha.glucosidase..3.2.1.20.",
  "RXN.14282...expasy..Alpha.glucosidase..3.2.1.20.",
  "RXN.14283...expasy..Alpha.glucosidase..3.2.1.20.",
  "RXN.15910...expasy..Alpha.glucosidase..3.2.1.20.",
  "RXN0.5183...expasy..Alpha.glucosidase..3.2.1.20.",
  "GLUCOSE.6.PHOSPHATE.1.EPIMERASE.RXN...expasy..Glucose.6.phosphate.1.epimerase..5.1.3.15.",
  "ATP.ADENYLYLTRANSFERASE.RXN...expasy..ATP.adenylyltransferase..2.7.7.53.",
  "D.XYLULOSE.REDUCTASE.RXN...expasy..D.xylulose.reductase..1.1.1.9."
)

volcano_df$color_group <- factor(
  ifelse(volcano_df$feature %in% blue_pathways, "B. longum",
         ifelse(volcano_df$feature %in% green_pathways, "B. breve", "Other")),
  levels = c("B. longum", "B. breve", "Other")
)

# Shapes
carb_deg_genes <- c(
  "MALTODEXGLUCOSID.RXN...expasy..Alpha.glucosidase..3.2.1.20.",
  "RXN.14281...expasy..Alpha.glucosidase..3.2.1.20.",
  "RXN.14282...expasy..Alpha.glucosidase..3.2.1.20.",
  "RXN.14283...expasy..Alpha.glucosidase..3.2.1.20.",
  "RXN.15910...expasy..Alpha.glucosidase..3.2.1.20.",
  "RXN0.5183...expasy..Alpha.glucosidase..3.2.1.20.",
  "X3.2.1.24.RXN...expasy..Alpha.mannosidase..3.2.1.24.",
  "RXN.18913...expasy..Alpha.mannosidase..3.2.1.24.",
  "X3.2.1.68.RXN...expasy..Isoamylase..3.2.1.68.",
  "RXN.4301...expasy..Isoamylase..3.2.1.68.",
  "RXN.14380...expasy..Isoamylase..3.2.1.68.",
  "RXN.12280...expasy..Isoamylase..3.2.1.68.",
  "GLUCOSE.6.PHOSPHATE.1.EPIMERASE.RXN...expasy..Glucose.6.phosphate.1.epimerase..5.1.3.15.",
  "L.FUCONATE.HYDRATASE.RXN...expasy..L.fuconate.dehydratase..4.2.1.68.",
  "X3.1.3.46.RXN...expasy..Fructose.2.6.bisphosphate.2.phosphatase..3.1.3.46."
)

volcano_df$shape_group <- factor(
  ifelse(volcano_df$feature %in% carb_deg_genes,
         "Carbohydrate degradation",
         "Other"),
  levels = c("Carbohydrate degradation", "Other")
)

# Volcano plot
p <- ggplot(volcano_df, aes(x = coef, y = -log10(qval))) +
  geom_hline(yintercept = -log10(0.05), col = "red", linewidth = 1) +
  annotate("text", x = min(volcano_df$coef) * 1.2, y = 1.2, 
           label = "adj_pvalue = 0.05", col = "red", size = 20, hjust = 0) +
  
  geom_text_repel(
    data = highlighted_points,
    aes(x = coef, y = -log10(qval), label = pretty_label),
    parse = FALSE,
    force = 50,
    force_pull = 0.5,
    box.padding = 0.5,
    point.padding = 4.0,
    max.overlaps = Inf,
    min.segment.length = 0,
    segment.size = 0.8,
    segment.curvature = 0,
    segment.angle = 50,
    segment.ncp = 1,
    size = 15,
    inherit.aes = FALSE,
    max.iter = 5000,
    max.time = 3,
    nudge_y = 1
  ) +
  
  geom_point(
    aes(fill = color_group, size = dominant_median_contribution, shape = shape_group),
    color = "black", stroke = 0.6, alpha = 0.9
  ) +
  
  labs(
    x = "Log2 Fold Change (Placebo (A) vs. Supplementation (B, C, D))",
    y = "-Log10 Adjusted P-Value",
    fill = "Dominant species",
    size = "Median CPM abundance"
  ) +
  
  scale_fill_manual(
    values = c(
      "B. longum" = "slateblue",
      "B. breve" = "#5f9e6e",
      "Other" = "black"
    ),
    labels = c(
      "B. longum" = "*B. longum*",
      "B. breve" = "*B. breve*",
      "Other" = "Other"
    )
  ) +
  
  scale_shape_manual(
    values = c(
      "Carbohydrate degradation" = 24,  # triangle
      "Other" = 21                      # circle
    )
  ) +
  
  scale_size_continuous(
    range = c(4, 24),
    breaks = c(1, 10, 20, 30, 40, 50)
  ) +
  
  guides(
    fill = guide_legend(
      override.aes = list(shape = 22, color = "black", size = 10),
      order = 1,
      title.position = "top"
    ),
    size = guide_legend(
      override.aes = list(shape = 21, fill = "grey", color = "black"),
      order = 2,
      title.position = "top"
    ),
    shape = guide_legend(
      override.aes = list(size = 10, fill = "grey"),
      order = 3,
      title = "Functional category",
      title.position = "top"
    )
  ) +
  
  theme_classic() +
  theme(
    axis.title.x = element_markdown(size = 70, face = "bold"),
    axis.title.y = element_markdown(size = 70, face = "bold"),
    axis.text = element_markdown(size = 62),
    plot.title = element_markdown(size = 80, hjust = 0.5, face = "bold"),
    legend.title = element_markdown(size = 62, face = "bold"),
    legend.text = element_markdown(size = 58),
    legend.position = c(0.02, 0.77),
    legend.justification = c("left", "top"),
    legend.box = "vertical",
    plot.margin = margin(t = 40, r = 20, b = 20, l = 20)
  )

# Save
#ggsave("volcano_plot_gene_families_babies_supplementation_prev-filtered-no24mo_updated_19.12.2025.png",
#       plot = p, width = 36, height = 26, dpi = 800, limitsize=FALSE)