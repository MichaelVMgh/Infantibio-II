#load libraries
library("vegan")
#install.packages("fossil")
#library("fossil")
library("stringr")
library("tidyverse")
library("ggpubr")
library("ggplot2")
library("dplyr")
library(grid)

#load df
getwd()
setwd("C:/Users/micha/Desktop/modified_figures")
df = read.table("merged_abundance_table_species_babies.txt", sep="\t", header=TRUE, row.names = 1, check.names = FALSE)
mf = read.table("mapping_file_original_babies.txt", sep="\t", header=TRUE)
#assign rownames using the sample ids
rownames(mf)=mf$sampleid
#assign variables to factors
mf$TimePoint <- factor(mf$TimePoint, levels=c("1mo", "3mo", "7mo", "12mo", "24mo"))
mf$HistoryAtTimePoint <- factor(mf$HistoryAtTimePoint, levels=c("A", "B", "C", "D", "BF", "infant_24mo"))


#BETA DIVERSITY
#calculate the Bray-Curtis dissimilarity
braycurtis <- vegdist(t(df), "bray")
#Prepare the data for the PCoA plot
#calculate the principal coordinates (axis) and the eigenvalue to infer the variation explained by each principal coordinate (PCs)
#set k=2 to retrieve the first two PCs for a 2D visualization
tmp_pcoa <- cmdscale(braycurtis, k=2, eig=TRUE)

#In this scatter plot each sample will represented by a point, the distance between the points indicates how similar the community composition between the samples is
#Save the coordinates for each sample as a dataframe
my_points = data.frame(tmp_pcoa$points)
colnames(my_points) = c("PC1", "PC2")

# Add the Timepoint information for each sample
my_points$TimePoint = mf[row.names(my_points), "TimePoint"]
# We will also save the Feeding of each sample to the dataframe
#my_points$Feeding_Overall = mf[rownames(my_points),]$Feeding_Overall
# We will also save the HistoryAtTimePoint of each sample to the dataframe
my_points$HistoryAtTimePoint = mf[rownames(my_points),]$HistoryAtTimePoint
# We will also save the TypeOfBirth of each sample to the dataframe
my_points$TypeOfBirth = mf[rownames(my_points),]$TypeOfBirth
# We will also save the Sex of each sample to the dataframe
my_points$Sex = mf[rownames(my_points),]$Sex
my_points$richness = mf[rownames(my_points),]$richness
my_points$shannon <- mf[rownames(my_points),]$shannon
my_points$bifido <- mf[rownames(my_points),]$g__Bifidobacterium

# Calculate how much variation is explained by the first and second PC
# In case some eigenvalues are negative, the abs() function will return the absolute value. These values should be small and this should _not_ be the case for the first few PCs.
pcoa_percentage = round(tmp_pcoa$eig*100 / sum(abs(tmp_pcoa$eig)), 1)
pcoa_percentage = paste(colnames(my_points), "(", paste( as.character(pcoa_percentage), "%", ")", sep="") )
pcoa_percentage = str_replace(pcoa_percentage, "axis", "PCoA")
head(my_points)

# Set factor levels directly to match histp_colors names
my_points$HistoryAtTimePoint <- factor(my_points$HistoryAtTimePoint, 
                                       levels = c("A", "B", "C", "D", "BF", "infant_24mo"))

timepoint_colors <- c(
  "1mo" = "lightgrey",
  "3mo" = "#008080",
  "7mo" = "aquamarine",
  "12mo" = "darkviolet",
  "24mo" = "cornsilk"
)


# Define the custom colors for the plot
histp_colors <- c(
  "A" = "#F9BC05",
  "B" = "#9B3350",
  "C" = "#71884A",
  "D" = "#042B44",
  "BF" = "#c0c0c0",
  "infant_24mo" = "#FFF8DC"
)

# Create a custom legend labels to replace "infant_24mo" with "Infant at month 24"
legend_labels <- c(
  "A" = "A",
  "B" = "B",
  "C" = "C",
  "D" = "D",
  "BF" = "BF",
  "infant_24mo" = "Infant at month 24"  # Rename this category
)




# PCoA PLOT WITH BIFIDOBACTERIUM RELATIVE ABUNDANCE BY COLOR

library(ggplot2)
library(scales)  # For gradientn

# Use a green palette
green_palette <- c("#f7fcf5", "#e5f5e0", "#c7e9c0", "#a1d99b", "#74c476", "#41ab5d", "#238b45", "#006d2c", "#00441b")

# Create the PCoA plot
plt <- ggplot(data = my_points, aes(x = PC1, y = PC2, 
                                    fill = bifido,  # Map the fill to HistoryAtTimePoint
                                    #size = shannon, 
                                    label = row.names(my_points))) +
  geom_point(alpha = 0.85, shape = 21, stroke = 1.2, color = "black", size=6) +  # Set edge color to black
  theme_bw() +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.text.x = element_text(size = 24),
    axis.text.y = element_text(size = 24),
    axis.title.x = element_text(size = 30),
    axis.title.y = element_text(size = 30),
    plot.title = element_text(size = 33, hjust = 0.5),
    legend.title = element_text(size = 27),
    legend.text = element_text(size = 24),
    legend.key.size = unit(1, "cm"),  # Adjust size of legend markers
    legend.key.width = unit(1.5, "cm"),
    plot.caption = element_text(size = 14, hjust = 1.4)  # Adjust caption size
  ) +
  #scale_fill_manual(values = histp_colors, labels = legend_labels) +  # Use histp_colors for the fill and apply custom labels
  #scale_size_continuous(range = c(1, 20),
  #                      breaks = c(1, 2, 3, 4)) +
  scale_fill_gradientn(
    colours = green_palette,
    name = "Bifidobacterium rel. abundance"
  )+
  
  labs(
    title = "PCoA of all infants (Bray-Curtis)",
    x = pcoa_percentage[1],
    y = pcoa_percentage[2],
    fill = "Bifidobacterium rel. abundance"
  ) +
  guides(
    #fill = guide_legend(override.aes = list(size = 5), order = 1)
    fill = guide_colorbar(barwidth = 1.5, barheight = 10, order = 1)# Increase size of legend color markers (circles)
  )

# Save
#png("pcoa_bifido_22.05.png", width = 15, height = 10, units = "in", res = 800)
#ggsave("pcoa_bifido_22.05.png", plot = plt, width = 15, height = 10, dpi = 800)
print(plt)
dev.off()




# PCoA PLOT WITH FEEDING GROUPS AND SHANNON DIVERSITY INDEX

library(ggplot2)
library(scales)  # For gradientn


# Create the PCoA plot
plt <- ggplot(data = my_points, aes(x = PC1, y = PC2, 
                                    fill = HistoryAtTimePoint,  # Map the fill to HistoryAtTimePoint
                                    size = shannon, 
                                    label = row.names(my_points))) +
  geom_point(alpha = 0.85, shape = 21, stroke = 1.2, color = "black") +  # Set edge color to black
  theme_bw() +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.text.x = element_text(size = 24),
    axis.text.y = element_text(size = 24),
    axis.title.x = element_text(size = 30),
    axis.title.y = element_text(size = 30),
    plot.title = element_text(size = 33, hjust = 0.5),
    legend.title = element_text(size = 27),
    legend.text = element_text(size = 24),
    legend.key.size = unit(1, "cm"),  # Adjust size of legend markers
    legend.key.width = unit(1.5, "cm"),
    plot.caption = element_text(size = 14, hjust = 1.4)  # Adjust caption size
  ) +
  scale_fill_manual(values = histp_colors, labels = legend_labels) +
  scale_size_continuous(range = c(1, 20),
                        breaks = c(1, 2, 3, 4)) +
  #scale_fill_gradientn(
  #  colours = green_palette,
  #  name = "Bifidobacterium rel. abundance"
  #)+
  
  labs(
    title = "PCoA of all infants (Bray-Curtis)",
    x = pcoa_percentage[1],
    y = pcoa_percentage[2],
    fill = "Feeding group", 
    size = "Shannon index",
    #caption = caption_text 
  ) +
  guides(
    fill = guide_legend(override.aes = list(size = 5), order = 1)
    #  fill = guide_colorbar(barwidth = 1.5, barheight = 10, order = 1)
  )

# Save
#png("pcoa_bifido_22.05.png", width = 15, height = 10, units = "in", res = 800)
#ggsave("pcoa_histp_shannon_30.05.2026.png", plot = plt, width = 15, height = 10, dpi = 800)
print(plt)
dev.off()






# PCoA PLOT WITH SAMPLING TIME POINTS AND RICHNESS

timepoint_colors <- c(
  "1mo" = "#ea593c",
  "3mo" = "#af5b9f",
  "7mo" = "#901440",
  "12mo" = "#571845",
  "24mo" = "#7bca9f"
)

# Replace time point labels
legend_labels <- c(
  "1mo" = "1",
  "3mo" = "3",
  "7mo" = "7",
  "12mo" = "12",
  "24mo" = "24"
)


# Create the plot
plt <- ggplot(data = my_points, aes(x = PC1, y = PC2, 
                                    fill = TimePoint,
                                    size = richness, 
                                    label = row.names(my_points))) +
  geom_point(alpha = 0.85, shape = 21, stroke = 1.2, color = "black") +
  theme_bw() +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.text.x = element_text(size = 24),
    axis.text.y = element_text(size = 24),
    axis.title.x = element_text(size = 30),
    axis.title.y = element_text(size = 30),
    plot.title = element_text(size = 33, hjust = 0.5),
    legend.title = element_text(size = 27),
    legend.text = element_text(size = 24),
    legend.key.size = unit(1, "cm"),
    legend.key.width = unit(1.5, "cm"),
    #plot.caption = element_text(size = 21, hjust = 1.4)  # Style the caption
  ) +
  scale_fill_manual(values = timepoint_colors, labels = legend_labels) +
  scale_size_continuous(range = c(1, 20), breaks = c(12, 100, 200, 300)) +
  labs(
    title = "PCoA of all infants (Bray-Curtis)",
    x = pcoa_percentage[1],
    y = pcoa_percentage[2],
    fill = "Time point (months)",
    size = "Richness (number of species)",
    #caption = caption_text  
  ) +
  guides(
    fill = guide_legend(override.aes = list(size = 5), order = 1)
  )

# Save
#ggsave("pcoa_tp_richness_final_22.05.png", plot = plt, width = 15, height = 10, dpi = 800)