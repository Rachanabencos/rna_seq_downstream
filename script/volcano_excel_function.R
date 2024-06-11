### Helper function for Deseq2 analysis 

#load the libraries
library(EnhancedVolcano)
library(openxlsx)
library(svglite)



### function to save DESeq2 output as XLSX
saveToExcel <- function(contrast , base, inputdf, log2fc=1)  {
  
  wb <- openxlsx::createWorkbook()
  
  # Function to add a worksheet and write data table
  addWriteDataTable <- function(wb, sheetname, data) {
    addWorksheet(wb, sheetname)
    writeDataTable(wb, sheetname, data, tableStyle = "TableStyleLight9")
  }
  
  # Original data
  sheetname <- substr(paste(contrast, 'vs', base, sep="_"), 0, 30)
  addWriteDataTable(wb, sheetname, inputdf)
  
  
  # Filtered padj < 0.05
  sheetname <- substr(paste('Filtered padj', 0.05), 0, 30)
  resp <- inputdf %>% filter(padj < 0.05)
  addWriteDataTable(wb, sheetname, resp)
  
  # Filtered log2FoldChange ± 0
  sheetname <- substr(paste('Filtered log2FoldChange±', log2fc), 0, 30)
  resplog2fc <- resp %>% filter(log2FoldChange > log2fc | log2FoldChange < -log2fc)
  addWriteDataTable(wb, sheetname, resplog2fc)
  
  # UPREG in contrast
  sheetname <- substr(paste('UPREG in', contrast), 0, 30)
  addWriteDataTable(wb, sheetname, resplog2fc %>% filter(log2FoldChange > log2fc))
  
  # DOWNREG in contrast
  sheetname <- substr(paste('DOWNREG in', contrast), 0, 30)
  addWriteDataTable(wb, sheetname, resplog2fc %>% filter(log2FoldChange < -log2fc))
  #define the output directory
  output_dir <- "data/DEG_Analysis"
  # Save workbook
  file_path <- file.path(output_dir, paste(contrast, "vs", base, "deseq2_result.xlsx", sep = "_"))
  saveWorkbook(wb, file_path, overwrite = TRUE)
}




### function to draw Volcano Plots
drawVolcano <- function(contrast, base, inputdf, log2fc) {
  
  #Select these columns: Ensemble Gene ID, external_gene_name, deg table like p, adjp, fc, basemean, uniprotID
  inputdf$label <- paste(inputdf$Row.names, inputdf$external_gene_name, sep = "_")
  
  colnames(inputdf)[1] <- "gene_ID"
  
  inputdf <- inputdf[!duplicated(inputdf),]
 
  up <- inputdf %>% filter(log2FoldChange > log2fc & padj < 0.05)
  down <- inputdf %>% filter(log2FoldChange < (-log2fc) & padj < 0.05)
  ### Extracting labels to show in volcano
  labels = inputdf$label
  g11 <- up %>% arrange(up$pvalue)
  g12 <- up %>% arrange(-up$log2FoldChange)
  g21 <- down %>% arrange(down$pvalue)
  g22 <- down %>% arrange(down$log2FoldChange)
  lab <- c(g11$label[1:3], g12$label[1:3],
           g21$label[1:3], g22$label[1:3])
  inputdf <- inputdf[c(1:7,10,11)]
  ### Volcano plot
  
  p <- EnhancedVolcano(inputdf, lab = labels,
                       selectLab = lab,
                       x = 'log2FoldChange', title = paste(contrast, "vs", base , sep = " "),
                       FCcutoff = log2fc, pCutoff = 0.05, y = 'padj',
                       labSize = 5.0, labCol = 'black', labFace = 'bold',
                       boxedLabels = TRUE, legendIconSize = 4.0,
                       drawConnectors = TRUE, widthConnectors = 1, colConnectors = 'black')
  
  base_name <- paste0(contrast, "_Vs_" ,base ,"_")
  # Define the directory path
  directory_path <- "images/volcano_plot"
  
  # Save the plot in SVG format
  ggsave(filename = file.path(directory_path, paste0("volcano_plot_", base_name, ".svg")), plot = p, width = 14, height = 16)
  
  # Save the plot in PDF format
  ggsave(filename = file.path(directory_path, paste0("volcano_plot_", base_name, ".pdf")), plot = p, width = 14, height = 16)
  
  # Save the plot in PNG format
  ggsave(filename = file.path(directory_path, paste0("volcano_plot_", base_name, ".png")), plot = p, width = 14, height = 16, dpi = 300)
  
  
  saveToExcel(contrast = contrast, base = base, inputdf = inputdf, log2fc = log2fc)
}