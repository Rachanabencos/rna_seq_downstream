#Load libraries 
library(DESeq2)
library(ComplexHeatmap)
library(dplyr)
library(viridis)

###Draw Heatmap function 

draw_comparison_heatmap <- function(base, contrast, sample_info, normalized_counts, res_df) {
  # Filter sample info
  samplesInfo <- sample_info[sample_info$condition == base | sample_info$condition == contrast,]
  
  # Filter normalized counts
  normalized_counts <- normalized_counts |> as.data.frame() |> 
    dplyr::select(c(samplesInfo$Sample))
  
  # Filter data based on log2FoldChange and padj
  filtered_data <- subset(res_df, abs(log2FoldChange) >= 1 & padj <= 0.05)
  
  # Separate into upregulated and downregulated genes
  upregulated_genes <- subset(filtered_data, log2FoldChange > 0)
  downregulated_genes <- subset(filtered_data, log2FoldChange < 0)
  
  # Select the top 25 upregulated and downregulated genes based on log2FoldChange
  top_upregulated <- head(upregulated_genes[order(-upregulated_genes$log2FoldChange), ], 25)
  top_downregulated <- head(downregulated_genes[order(downregulated_genes$log2FoldChange), ], 25)
  merged_df <- rbind(top_upregulated, top_downregulated)
  
  # Select normalized count for top DEG
  top50<- normalized_counts[rownames(normalized_counts) %in% merged_df$Row.names,]
  top50 <- merge(top50, merged_df, by.x = "row.names", by.y = "Row.names", all.x = TRUE)
  top50$label <- paste(top50$Row.names, top50$external_gene_name, sep = "_")
  rownames(top50) <- top50$label
  top50 <- subset(top50, select = -c(Row.names, baseMean, log2FoldChange, lfcSE, stat,   
                                     pvalue, padj,ensembl_transcript_id_version, 
                                     ensembl_transcript_id,
                                     external_gene_name, uniprotswissprot, label))
  
  top_annotation_colors <- c("#1B9E77", "#E72988")  # Adjust colors as needed
  # Column colour
  splitColumn <- samplesInfo$condition
  
  p <- ComplexHeatmap::Heatmap(as.matrix(log2(top50+1)), name = "Normalized Count",
                               cluster_rows = TRUE,
                               cluster_columns = TRUE,
                               top_annotation = HeatmapAnnotation(foo = anno_block(gp = gpar(fill = top_annotation_colors))),
                               column_split = splitColumn,
                               show_row_dend = TRUE,
                               row_title_rot = 0,
                               rect_gp = gpar(col = 'black'))
  
  # Save the plot in PDF format
  base_name <- paste0(contrast, "_vs_", base)
  pdf(paste0("images/comparison_heatmap/", base_name, "_heatMap.pdf"), width = 8, height = 12)
  draw(p)
  dev.off()
  
  # Save the plot in PNG format
  png(paste0("images/comparison_heatmap/", base_name, "_heatMap.png"), width = 8, height = 12, res = 1200, units = 'in')
  draw(p)
  dev.off()
}
