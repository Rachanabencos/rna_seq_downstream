library(ComplexHeatmap)
library(dplyr)
library(tidyr)
library(ggplot2)
library(RColorBrewer)
#' Import the normalized count matrix
#' Rownames should be unique names 
#' It should have column with Ensemble gene ID concatenate with Gene symbol with _
#' Example
#' GenesID	DK-I_S50	DK-II_S54	DS-I_S51	DS-II_S55	DW-I_S387	DW-II_S388	external_gene_name	gene_ID
#' FBgn0000008	10.77714942	10.79683411	10.95488108	10.93677795	10.93157248	10.97613673	a	FBgn0000008_a
#' FBgn0000014	10.22072153	10.60321555	10.31405303	10.47865259	10.58712727	10.31634571	abd-A	FBgn0000014_abd-A
#' FBgn0000015	9.996561042	10.15745983	10.03664075	10.05550865	10.2317263	10.0541099	Abd-B	FBgn0000015_Abd-B
#' FBgn0000017	11.07765495	11.0340167	11.25842922	11.30287259	11.49128691	11.31571239	Abl	FBgn0000017_Abl
#' FBgn0000018	10.99866692	11.0340167	11.00409899	11.01726781	10.94464057	11.04539805	abo	FBgn0000018_abo


args <- commandArgs(trailingOnly = TRUE)
normalizedCount <- args[1]
# Check if the correct number of arguments are provided
if (length(args) != 4) {
  cat("Error: Please follow the below format! \n")
  cat("Usage: Rscript pathwayHeatmap.R DK_DW_normalizedCounts.csv 
      DK_vs_DW_UPREG_DEGs_enriched_kegg_new.tsv 
      DK_DW_Column_annotation.txt 0.005\n")
  quit(status = 1)
}

countMatrix <- read.csv(normalizedCount, 
                        header = T, sep = ',')
rownames(countMatrix) <- countMatrix$gene_ID



#' Import KEGG pathway output table
#' ID	Term_Description	Fold_Enrichment	occurrence	support	lowest_p	highest_p	Genes
# dme03010	Ribosome	6.190410282	10	0.020204082	6.33E-09	6.33E-09	RpS21, RpS25, RpL24-like, RpL24, RpL29, RpL36
# dme04624	Toll and Imd signaling pathway	6.696898396	10	0.095246734	1.57E-05	1.57E-05	CecA1, Jra, Rel
# dme04080	Neuroactive ligand-receptor interaction	13.64183007	10	0.020204082	0.003913805	0.003913805	alphaTry, epsilonTry, CG30031, deltaTry
# dme04214	Apoptosis - fly	2.407381776	10	0.067649282	0.007913189	0.015776845	Jra

KEGGPathwayinput<- args[2]


pvalue <- as.numeric(commandArgs(trailingOnly = TRUE)[3])



KEGG_data <- read.table(KEGGPathwayinput, 
                        header = T, sep = '\t')
colnames(KEGG_data)[8] <- 'Genes'


#filter the pathways based on the p value 
filtered_data <- KEGG_data %>%filter(highest_p <= pvalue)

top_terms <- filtered_data %>%
  top_n(15, wt = highest_p)

# separate multiple genes for each term into multi-line 
result <-  top_terms%>%
  tidyr::separate_rows(Genes, sep = ", ") %>%
  dplyr::select(Term_Description, Symbol = Genes)

# Merge both Pathway and normalized count matrix
merged_df <- merge(result, countMatrix, by.x = "Symbol", by.y = "external_gene_name",
                   all.x = FALSE, all.y = FALSE)
#Assign rownames
merged_df$Gene_ID_1 <- make.names(merged_df$Gene_ID, unique = TRUE)
rownames(merged_df) <- merged_df$Gene_ID_1
colnames(merged_df)[2] <- 'Pathway'
merged_df_sorted<-merged_df[order(merged_df$Pathway, 
                                  decreasing = TRUE), ] 
# Keep only the top 20 pathways
# merged_df_sorted <- head(merged_df[order(merged_df$Pathway, decreasing = TRUE), ], 0)

#Create a dataframe without alphabets only numerals : to create a matrix
dataFrame <- subset(merged_df_sorted, 
                    select = -c(Symbol, Pathway, Gene_ID, 
                                 Gene_ID_1)) |> as.matrix()


columnAnnotatioInput <- args[4]
annotation_col <- read.table(columnAnnotatioInput, row.names = 1, header = T)
splitColumn <- annotation_col$disease
splitRow <- merged_df_sorted$Pathway





#Create gaps


if(length(unique(annotation_col$Sample)) < 3){
  top_annotation_colors <-c("#1B9E77", "#D95F02")
}else{
  top_annotation_colors <- brewer.pal(length(unique(annotation_col$Sample)), 
                                      "Set2")
}
if(merged_df_sorted$Pathway |> unique() |> length() < 3){
  left_annotation_colors <- c("#7570B3", "#E7298A")
}else{
  left_annotation_colors <- c(
  "#1f77b4", "#ff7f0e", "#2ca02c", "#d62728", "#9467bd",
  "#8c564b", "#e377c2", "#7f7f7f", "#bcbd22", "#17becf",
  "#808000", "#ffbb78", "#98df8a", "#ff9896", "#c5b0d5",
  "#c49c94", "#f7b6d2", "#c7c7c7", "#dfff00", "#9edae5",
  "#008000","#00FF00","#00FFFF","#FF00FF",
  "#f17ff5", "#800080" 
   )
  
}


p<-Heatmap(dataFrame, name = "Normalized Count",
        cluster_rows = F,
        cluster_columns = F,
        top_annotation = HeatmapAnnotation(foo = anno_block(gp = gpar(fill = top_annotation_colors))),
        column_split = splitColumn,
        left_annotation = rowAnnotation(foo = anno_block(gp = gpar(fill = left_annotation_colors))),
        row_split = splitRow,
        show_row_dend =T,
        row_title_rot = 0,
        rect_gp = gpar(col = 'black'),
        heatmap_height = unit(30, "in"),
        
        
  )

# Get the base name of the input file
base_name <- tools::file_path_sans_ext(basename(KEGGPathwayinput))


# Define the directory path
output_dir <- "data/DEG_Analysis/Enrichment_Analysis/pathfindR/pathway_heatmaps/"

# Save the plot in PDF format with the base name of the input file
pdf(paste0(output_dir, base_name, "heatMap.pdf"), width = 14, height = 34)
draw(p)
dev.off()

# Save the plot in PNG format with the base name of the input file
#png(paste0(output_dir, base_name, "heatMap.png"), width = 12, height = 18, units = "in", 
    #res = 1500)
#draw(p)
#dev.off()

