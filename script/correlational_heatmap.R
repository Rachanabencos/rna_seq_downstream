sampleDists <- as.matrix(cor( t(assay(rld)) ))
write.table(sampleDists,file = "images/rlog_Normalized_Values.txt",sep = "\t")
ht = Heatmap(
  sampleDists,
  col = col4(100),
  heatmap_legend_param  =  list(title = NULL)
)
pdf(file = "images/SampleDistanceRlog.pdf", height = 5, width = 5.5)
draw(ht,heatmap_legend_side  =  "left")
dev.off()



# Calculate Pearson correlation matrix
pearson_corr <- cor(t(assay(rld)))

# Calculate Spearman correlation matrix
spearman_corr <- cor(t(assay(rld)), method = "spearman")

# Combine the correlation matrices
combined_corr <- (pearson_corr + spearman_corr) / 2  # Average of Pearson and Spearman correlations

# Write the combined correlation matrix to a file
write.table(combined_corr, file = "images/rlog_Normalized_Values_combined.txt", sep = "\t")

# Create a heatmap with correlational values in each box
ht <- Heatmap(
  combined_corr,
  col = col4(100),
  heatmap_legend_param = list(title = NULL),
  top_annotation = HeatmapAnnotation(
    annotation_legend_param = list(title = "Correlation"),
    annotation_name_gp = gpar(fontsize = 8),
    annotations = as.data.frame(combined_corr)
  )
)

# Save the heatmap as a PDF
pdf(file = "images/SampleDistanceRlog_combined.pdf", height = 5, width = 5.5)
draw(ht, heatmap_legend_side = "left")
dev.off()
