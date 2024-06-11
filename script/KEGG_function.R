library(clusterProfiler)
library(pathfindR)
library(dplyr)


run_kegg_pathfindR <- function(file, sheetnumber, degType, org) {
  
  sheetnumber <- as.numeric(sheetnumber)
  data <- readxl::read_excel(file,sheet = sheetnumber)
  #Name of the output file
  name <- basename(file)
  name <- tools::file_path_sans_ext(name)
  name <- paste0(name, "_", degType)
  
  #Define the output directory
  output_dir <- "data/DEG_Analysis/Enrichment_Analysis/KEGG"
  
  
  
  geneID<-data$uniprotswissprot
  
  kk <- enrichKEGG(
    gene = bitr_kegg(
      geneID = geneID,
      fromType = "uniprot",
      toType = "kegg",
      organism = org
    )[,2],
    organism = org,
    pvalue = 0.05
  )
  
  kkk <- enrichKEGG(
    gene = bitr_kegg(
      geneID = geneID,
      fromType = "uniprot",
      toType = "kegg",
      organism = org
    )[,2],
    organism = org,
    pvalue = 1
  )
  
  
  
  outTable <- paste(name, ".KEGG.xls", sep="")
  write.table(kkk, file=outTable, sep="\t", eol="\n", quote=FALSE)
  
  outPNG <- paste(name, ".png", sep="")
  outPDF <- paste(name, ".pdf", sep="")
  
  # save as PNG
  outPNG <- file.path(output_dir, paste0(name, "_plot.png"))
  png(outPNG, width = 900, height = 700, units = "px", pointsize = 12)
  dotplot(kk, showCategory = 20)
  dev.off()
  
  # Save as PDF
  # Save as PDF in the output directory
  outPDF <- file.path(output_dir, paste0(name, "_plot.pdf"))
  pdf(outPDF, width = 9, height = 7)
  dotplot(kk, showCategory = 20)
  dev.off()
  
  
  
  
  hsa_kegg_genes <- readRDS('data/hsa_kegg_genes.rds')
  hsa_kegg_descriptions <- readRDS('data/hsa_kegg_descriptions.rds')
  
  df2 <- as.data.frame(data) |>
    #select(gene_names, log2FoldChange, padj) 
    dplyr::select(external_gene_name, log2FoldChange, padj) 
  df2 <- dplyr::distinct(df2, .keep_all = T)
  
  pathfindR <- pathfindR::run_pathfindR(input = df2, 
                                        n_processes = 10, 
                                        convert2alias = T,
                                        gene_sets = "Custom",
                                        custom_genes = hsa_kegg_genes,
                                        custom_descriptions = hsa_kegg_descriptions,
                                        pin_name_path = "data/hsa_pin.sif",
                                        plot_enrichment_chart = F)
  
  if(degType == "Up" || degType == "UP" || degType == "up"){
    pathfindR <- pathfindR |> dplyr::select(-c(Down_regulated))
  }else{
    pathfindR <- pathfindR |> dplyr::select(-c(Up_regulated))
  }
  
  new_dir <- "data/DEG_Analysis/Enrichment_Analysis/pathfindR"
  outTable <- file.path(new_dir, paste0(name, "_pathfindR.xls"))
  
  # Write the output table to the specified file path
  write.table(pathfindR, file=outTable, sep="\t", eol="\n", quote=FALSE)
}
  