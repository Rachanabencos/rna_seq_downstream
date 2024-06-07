#Load the required libraries

library(DESeq2)
library(ComplexHeatmap)
library(dplyr)
library(viridis)


#plot heatmap for all samples using the top 50 variable genes 
plotheatmap <- function(norm.counts, sample_info, t2g) {
  topVarGenes <- head(order(SparseArray::rowVars(SummarizedExperiment::assay(norm.counts)), 
                            decreasing = TRUE), 50)
  mat <- SummarizedExperiment::assay(norm.counts)[topVarGenes, ]
  mat <- as.data.frame(mat)
  mat$Gene_ID <- rownames(mat)
  mat <- dplyr::left_join(mat, t2g, by = c("Gene_ID" = "ensembl_gene_id"))
  mat$label <- paste(mat$Gene_ID, mat$external_gene_name, sep = "_")
  mat <- mat[!duplicated(mat$label), ]
  rownames(mat) <- mat$label
  cols_to_remove <- c("Gene_ID", "ensembl_transcript_id_version", 
                      "ensembl_transcript_id", "external_gene_name", 
                      "uniprotswissprot", "label")
  
  mat <- mat[, !names(mat) %in% cols_to_remove]
  
  mat <- dplyr::mutate_all(mat, function(x) as.numeric(as.character(x)))
  mat <- mat - rowMeans(mat)
  
  num_samples <- ncol(mat)
  top_annotation_colors <- viridis(num_samples)
  
  
  splitColumn <- sample_info$condition
  
  p <- ComplexHeatmap::Heatmap(as.matrix(mat), name = "Normalized Count",
                               cluster_rows = TRUE,
                               cluster_columns = FALSE,
                               top_annotation = HeatmapAnnotation(foo = anno_block(gp = grid::gpar(fill = top_annotation_colors))),
                               column_split = splitColumn,
                               show_row_dend = TRUE,
                               row_title_rot = 0,
                               column_title_rot = 0,
                               column_title_gp = gpar(fontsize = 10),
                               rect_gp = gpar(col = 'black')
  )
  
  return(p)
}

