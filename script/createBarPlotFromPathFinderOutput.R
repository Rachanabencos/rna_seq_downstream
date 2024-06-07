# Load necessary library
library(ggplot2)

# Read input file name from command line arguments
args <- commandArgs(trailingOnly = TRUE)
input_file <- args[1]
inputColor <- args[2]
# Check if the correct number of arguments are provided
if (length(args) != 2) {
  cat("Usage: Rscript createBarPlotFromPathFinderOutput.R DK_vs_DW_DOWNREG_DEGs_enriched_kegg_new.tsv brown\n")
  quit(status = 1)
}
# Read the data from the input file
data <- read.table(input_file, header = TRUE, sep = '\t')
data$log2FoldEnrichment <- log2(data$Fold_Enrichment)
colnames(data)[6] <- 'p.value'
colnames(data)[8] <- 'Genes'
# Calculate the number of genes
data$Gene_Count <- sapply(strsplit(data$Genes, ", "), length)

# Create the barplot using ggplot2
p <- ggplot(data, aes(x = reorder(Term_Description, -log2FoldEnrichment), y = log2FoldEnrichment)) +
  geom_bar(stat = "identity", position = "dodge", color = inputColor, fill = inputColor, width = 0.5) +  # Decrease the width of the bars
  geom_text(aes(label = paste("(", "n =", Gene_Count, ",", "p.value =", p.value, ")")), 
            position = position_dodge(width = 0.9), vjust = -0.5, hjust = -0.01, size = 4, fontface = "bold") + # Set fontface to "bold" for the text
  coord_flip() +  # Flip the coordinates to create horizontal bars
  labs(x = "", y = "Fold Enrichment") +
  scale_y_continuous(expand = c(0, 0), limits = c(0, max(data$log2FoldEnrichment) + 1)) +  # Manipulate y-axis scale
  theme_classic() + # Apply classic theme
  theme(legend.position = "none",          # Remove the legend 
        axis.text.x = element_text(size = 12,  face = "bold"),  # Set x-axis text to bold
        axis.text.y = element_text(size = 12, face = "bold"),
        axis.title = element_text(size = 14, face = "bold"))

# Get the base name of the input file
base_name <- tools::file_path_sans_ext(basename(input_file))

# Save the plot in SVG format with the base name of the input file
ggsave(filename = paste0(base_name, "barPlot.svg"), plot = p, width = 16, height = 24)

# Save the plot in PDF format with the base name of the input file
ggsave(filename = paste0(base_name, "barPlot.pdf"), plot = p, width = 16, height = 24)

# Save the plot in PNG format with the base name of the input file
ggsave(filename = paste0(base_name, "barPlot.png"), plot = p, width = 16, height = 24)
