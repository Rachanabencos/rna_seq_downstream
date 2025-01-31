# RNA Seq downstream data analysis
## Postprocessing samlmon results

### Downstream RNA seq analysis after we get the quant files by salmon quant

### Perform Transcript abundance estimation by mapping transcripts to the gene by tximport

```{r echo=T}
#read the mart file from ensemble biomaRT containing columns:
#ensembl_transcript_id_version	ensembl_transcript_id	ensembl_gene_id	external_gene_name	#uniprotswissprot

t2g<- utils::read.table("data/mart_export_human.txt", header = T, sep = '\t')
#select these columns to map the transcripts id to the gene
t2g2 <- t2g |> dplyr::select(ensembl_transcript_id_version,ensembl_gene_id,external_gene_name )

#set the base directory having the quant files 
base_dir <- 'data/quant_files/'

#list all the sample-specific salmon sub-directories 
sample_ids <- list.dirs(base_dir, full.names = FALSE, 
                        recursive = FALSE)

#So what we want to do now is create paths to each quant.sf file that is in each sample_id.
#This can be done by combining the base_dir, each sample_id directory, and 'quant.sf'
#For example, the path to the first file will be data/quant_files/C1/quant.sf
salm_dirs <- sapply(sample_ids, function(id) file.path(base_dir, id, 'quant.sf'))
salm_dirs

#remove duplicates from the t2g2 dataframe
t2g2 <- t2g2[!duplicated(t2g2),]

#Importing and Aggregating Transcript-Level Abundance Estimates with tximport
txi <- tximport::tximport(salm_dirs, type = 'salmon', 
                tx2gene = t2g2, dropInfReps = TRUE,
                countsFromAbundance = 'lengthScaledTPM')

```

### Creation of DESEq object for Differential gene expression analysis and perform VST normalisation

```{r echo=T}

#Read the sample metadata file
#make sure the sample info table has sample names in the first column and condition/replicates in #the second column
#format of sample table:
#Sample	condition
#C1	Control
#C2	Control

sample_info<- utils::read.table("data/sample_info.txt", header = T, sep = "\t")
#Rename the sample info column containing the condition
colnames(sample_info)[2] <- "condition"

## Create the DESeq2 object
ddsTxi <- DESeq2::DESeqDataSetFromTximport(txi = txi,
                                   colData = sample_info,
                                   design = ~ condition)
ddsTxi <- ddsTxi[rowSums(DESeq2::counts(ddsTxi))>=10,]
ddsTxi <- DESeq2::DESeq(ddsTxi)
#to perform the normalization of DESeq count data
norm.counts <- DESeq2::varianceStabilizingTransformation(ddsTxi, blind = F)
normalized_counts <- SummarizedExperiment::assay(norm.counts)
utils::write.table(normalized_counts, "data/output_files/normalized_counts.txt", sep="\t")
```
### Plot the PCA from the VST normalised counts
```{r echo=T}
#Create PCA plot 
pca <- BiocGenerics::plotPCA(norm.counts, intgroup = c("condition")) +
  ggplot2::geom_text(size = 2.5, ggplot2::aes(label = sample_info$Sample), vjust = 2, hjust = 0, check_overlap = TRUE) + ggplot2::coord_cartesian(clip = 'off')  
grDevices::tiff("images/PCA.tiff", units="in", width=10, height=8, res=600)
print(pca)
grDevices::dev.off()
```

### Plot heatmap for all samples from the top 50 variable genes
```{r echo =T}
#To create heatmap for all samples
#Source the helper function

source("script/heatmap_all_samples_function.R")

#Assuming norm.counts, sample_info, and t2g are defined in your environment

heatmap <- plotheatmap(norm.counts, sample_info, t2g)

#define the base name
base_name <- "All_samples"  # You may define base_name according to your requirements

# Save the plot in PDF format with the base name of the input file
# Define the base name for the input file
base_name <- "All_samples"  # You may define base_name according to your requirements

# Define the file path for the PDF
pdf_file <- paste0("images/", base_name, "_heatMap.pdf")

# Open the PDF device
grDevices::pdf(pdf_file, width = 8, height = 12)

# Print a message to confirm the PDF device is open
cat("PDF device opened successfully: ", pdf_file, "\n")

# Create the heatmap plot
ComplexHeatmap::draw(heatmap)

# Close the PDF device
grDevices::dev.off()

# Print a message to confirm the PDF device is closed
cat("PDF device closed successfully\n")

# Define the file path for the PNG
png_file <- paste0("images/", base_name, "_heatMap.png")

# Open the PNG device
grDevices::png(png_file, width = 8, height = 12, res = 1200, units = 'in')

# Print a message to confirm the PNG device is open
cat("PNG device opened successfully: ", png_file, "\n")

# Create the heatmap plot
ComplexHeatmap::draw(heatmap)

# Close the PNG device
grDevices::dev.off()

# Print a message to confirm the PNG device is closed
cat("PNG device closed successfully\n")


```


### DE Analysis :to create DEG table and volcano plots 

```{r echo=T}

#Create new folders to store the DEG tables and volcano plots 
dir.create("data/DEG_Analysis", recursive = TRUE)
dir.create("images/volcano_plot", recursive = TRUE)
dir.create("images/comparison_heatmap", recursive = TRUE)

#source the helper function
source("script/volcano_excel_function.R")
source("script/comparison_heatmap_function.R")

#Read the base contrast file formatted as:
#Base	Contrast 
#Control	HIV_mono

base_contrast<- read.table("data/base_contrast.txt", header= T, sep='\t')

for (i in 1:nrow(base_contrast)) {
  base <- base_contrast$Base[i]
  contrast <- base_contrast$Contrast[i]
  
  # Perform DE analysis
  res <- DESeq2::results(ddsTxi, contrast = c("condition", contrast, base),
                 independentFiltering = TRUE, alpha = 0.05,
                 pAdjustMethod = "BH")
  
  # Convert results to dataframe
  res_df <- as.data.frame(res)
  res_df <- merge(res_df, t2g, by.x = "row.names", by.y = "ensembl_gene_id", all.x = TRUE)
  res_df <- res_df[!duplicated(res_df$Row.names),]
  
  #call the volcano function to plot volcano files and save the DEG table in excel workbook
  
  drawVolcano(contrast = contrast, base = base, 
              inputdf = res_df, log2fc = 1)
  
  #call the comparison heatmap function to plot the heatmap for comparisons based on DEGs
  
  draw_comparison_heatmap(base = base, 
                         contrast = contrast, 
                         sample_info, 
                         normalized_counts, 
                         res_df)
  }


```
### GO Analysis 

```{r echo=TRUE}


#Create a new directory for the GO Analysis
dir.create("data/DEG_Analysis/Enrichment_Analysis/GO", recursive = TRUE)


#source the function
source("script/GO_analysis_function.R")


# Get a list of files in the DEG_analysis folder
files <- list.files('data/DEG_Analysis/', pattern = '_deseq2_result.xlsx', full.names = TRUE)


# Loop over each file and call the analysis function - upregulated
for (file in files) {
  run_go_analysis(file = file, sheetnumber = 4, degType = "Up", orgdb = org.Hs.eg.db)
}


# Loop over each file and call the analysis function - downregulated
for (file in files) {
  run_go_analysis(file = file, sheetnumber = 5, degType = "Down", orgdb = org.Hs.eg.db)
}



```
### KEGG Analysis

```{r echo=T}
#Create a new directory for the KEGG Analysis
dir.create("data/DEG_Analysis/Enrichment_Analysis/KEGG", recursive = TRUE)
dir.create("data/DEG_Analysis/Enrichment_Analysis/pathfindR", recursive = TRUE)

#source the function
source("script/KEGG_function.R")

# Get a list of files in the DEG_analysis folder
files <- list.files('data/DEG_Analysis/', pattern = '_deseq2_result.xlsx', full.names = TRUE)

#make sure the RDS file of gene set lists and PIN sif file for PathfindR analysis are present
#to get the gene set list
gsets_list <- pathfindR::get_gene_sets_list(
   source = "KEGG",
   org_code = "hsa")
## Save both as RDS files for later use
saveRDS(hsa_kegg_genes, "hsa_kegg_genes.RDS")
saveRDS(hsa_kegg_descriptions, "hsa_kegg_descriptions.RDS")


#to get the PIN sif file
path_to_pin_file <- pathfindR::get_pin_file(source = "BioGRID", org = "Homo_sapiens")
write.table(path_to_pin_file, file = "hsapiens_PIN.sif",
            col.names = FALSE, row.names = FALSE,
            sep = "\t", quote = FALSE)


#perform the KEGG Analysis
# Loop over each file and call the analysis function - up-regulated
for (file in files){
  run_kegg_pathfindR(file = file, sheetnumber= 4, degType="Up", org= "hsa")
}


# Loop over each file and call the analysis function - down-regulated
for (file in files) {
  run_kegg_pathfindR(file = file, sheetnumber = 5, degType = "Down", org = "hsa")}


```


### To create barplot from the pathfindR files
```{r echo=T}
#Create a new directory for the KEGG Analysis
dir.create("data/DEG_Analysis/Enrichment_Analysis/pathfindR/bar_plots", recursive = TRUE)

# path for the script 
script_path <- "script/createBarPlotFromPathFinderOutput.R"


# Path to your input file

input_dir <- "data/DEG_Analysis/Enrichment_Analysis/pathfindR/"
input_files <- list.files(input_dir, full.names = T, pattern = "_pathfindR.xls")


# Parameters
generate_color_palette <- function(n) {
  colors <- grDevices::colorRampPalette(c("red", "blue", "green", "yellow", "purple", "orange", "pink", "cyan", "magenta", "brown"))(n)
  return(colors)
}


color_palette <- generate_color_palette(length(input_files))


#apply the function
for (i in seq_along(input_files)) {
  file <- input_files[i]
  color <- color_palette[i]
  system(paste("Rscript", script_path, file, color))
}
```

### To create pathway heatmap from the pathfindR files
```{r echo=T}


#Create a new directory for the KEGG Analysis
dir.create("data/DEG_Analysis/Enrichment_Analysis/pathfindR/pathway_heatmaps", recursive = TRUE)



# Full path to your R script
script_path <- "script/pathwaycomplexHeatmap.R"

# Path to your inputs file
normalised_file <- "data/Normalised_files/HIV_mono_vs_Control.csv"

pathfindR_file<- "data/DEG_Analysis/Enrichment_Analysis/pathfindR/HIV_mono_vs_Control_deseq2_result_Up_pathfindR.xls"

annotation_file<- "data/annotation_files/HIV_mono_vs_Control_col_annotation.txt"

pvalue<-0.005

system(paste("Rscript", script_path, normalised_file, pathfindR_file, pvalue, annotation_file))
print(p)
```
