library(clusterProfiler)
library(pathfindR)
library(dplyr)
args <- commandArgs(trailingOnly = TRUE)
#args[1] 
#Assumes that the columns are formatted as such:
## Ensembl_GID	UniProt_ID
#OS01G0108600	Q8W0N9
#OS01G0113200	Q94E34
#args[2] : KEGG organism three letter code

file= args[1]
#pass the sheet number argument which have only signicant DEG Table
sheetnumber = args[2]

#Pass the DEG type eg Up or Down
degType = args[3]

sheetnumber <- as.numeric(sheetnumber)
data <- readxl::read_excel(file,sheet = sheetnumber)
#Pass the three letter KEGG code for the organism e.g hsa
org <- args[4]
#Name of the output file
name <- tools::file_path_sans_ext(file)
name <- paste0(name, "_",degType)


# ensembl <- biomaRt::useMart('ensembl', dataset = 'hsapiens_gene_ensembl')

# gene_symbols <- biomaRt::getBM(attributes = c("ensembl_gene_id", "uniprotswissprot"),
#                                filters = "ensembl_gene_id",
#                                values = data$Gene_ID,
#                                mart = ensembl)
geneID<-data$UniProtID

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

# Save as PNG
png(outPNG, width = 900, height = 700, units = "px", pointsize = 12)
dotplot(kk, showCategory = 20)
dev.off()

# Save as PDF
pdf(outPDF, width = 9, height = 7)
dotplot(kk, showCategory = 20)
dev.off()


# gsets_list <- pathfindR::get_gene_sets_list(
#   source = "KEGG",
#   org_code = "hsa"
# )

hsa_kegg_genes <- readRDS('hsa_kegg_genes.rds')
hsa_kegg_descriptions <- readRDS('hsa_kegg_descriptions.rds')

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
                                      pin_name_path = "hsa_PIN.sif",
                                      plot_enrichment_chart = F)

if(degType == "Up" || degType == "UP" || degType == "up"){
  pathfindR <- pathfindR |> dplyr::select(-c(Down_regulated))
}else{
  pathfindR <- pathfindR |> dplyr::select(-c(Up_regulated))
}


outTable <- paste(name, "_pathfindR.xls", sep="")
write.table(pathfindR, file=outTable, sep="\t", eol="\n", quote=FALSE)
