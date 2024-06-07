library(clusterProfiler)

library(ggplot2)
#library(org.Mm.eg.db) #for mouse
library(org.Hs.eg.db) 
#input list of genesymbol and organism database (eg org.Mm.eg.db)
library(readxl)
args <- commandArgs(trailingOnly = TRUE)
file= args[1]
#pass the sheet number argument which have only signicant DEG Table
sheetnumber = args[2]

#Pass the DEG type eg Up or Down
degType = args[3]

sheetnumber <- as.numeric(sheetnumber)
data <- readxl::read_excel(file,sheet = sheetnumber)
#pass the database e.g org.Hs.eg.db
orgdb <- args[4]
#Name of the output file
name <- tools::file_path_sans_ext(file)
name <- paste0(name, "_",degType)




gene<-data$external_gene_name

ego3 <- enrichGO(gene         = gene,
                 OrgDb         = orgdb,
                 keyType       = 'SYMBOL',
                 ont           = "CC",
                 pAdjustMethod = "BH",
                 pvalueCutoff  = 0.01,
                 qvalueCutoff  = 0.05)
# Shorten the term names on the y-axis and make Title case
ego3@result$Description <- substr(ego3@result$Description, 1, 30) # Adjust width as needed
ego3@result$Description <- tools::toTitleCase(ego3@result$Description)
ego3@result$Description <- stringr::str_wrap(ego3@result$Description, width = 60)
ego1 <- enrichGO(gene         = gene,
                 OrgDb         = orgdb,
                 keyType       = 'SYMBOL',
                 ont           = "CC",
                 pAdjustMethod = "BH",
                 pvalueCutoff  = 1,
                 qvalueCutoff  = 1)


outSVG_CC <- paste("CC_",name,".svg", sep="")

g_cc<-dotplot(ego3,showCategory=20 )+geom_point(shape = 21, aes(fill = p.adjust), alpha = 0.9, stroke = 0.5) + 
  scale_size_continuous(range = c(3, 10)) +  # Adjust size range
  scale_color_gradient(low = "blue", high = "red") +  # Blue to red color gradient for CorrectedP.Value
  scale_fill_gradient(low = "blue", high = "red") +  # Blue to red color gradient for fill
  theme_minimal() +  # Minimal theme for clarity
  labs(x = "GeneRatios", y = "") +  # Label the axes
  # ggtitle("Dot Plot of GeneRatio vs Term") +  # Add title
  theme(panel.grid.major.y = element_line(linewidth = .5, linetype = "solid", colour = 'gray'),
        panel.grid.major.x = element_line(linewidth = .5, linetype = "solid", colour = 'gray'),
        panel.border = element_rect(color = "black", fill = NA, linewidth = 1),  # Add outer border
        axis.text.y = element_text(size = 10, face = "bold"))


outPNG_CC <- paste("CC_", name, ".png", sep = "")
outPDF_CC <- paste("CC_", name, ".pdf", sep = "")

if(g_cc$data |> nrow() >=1){
  ggsave(outPNG_CC, plot = g_cc, width = 8, height = 6, dpi = 300)
  ggsave(outPDF_CC, plot = g_cc, width = 8, height = 6)
  ggsave(outSVG_CC, plot = g_cc)
}


outTable <- paste("CC_",name,".GO.xls", sep="")
write.table(ego1, file=outTable, sep="\t", eol="\n",quote=FALSE, 
            row.names = FALSE)

ego_BP <- enrichGO(gene         = gene,
                   OrgDb         = orgdb,
                   keyType       = 'SYMBOL',
                   ont           = "BP",
                   pAdjustMethod = "BH",
                   pvalueCutoff  = 0.01,
                   qvalueCutoff  = 0.05)
# Shorten the term names on the y-axis and make Title case
ego_BP@result$Description <- substr(ego_BP@result$Description, 1, 30) # Adjust width as needed
ego_BP@result$Description <- tools::toTitleCase(ego_BP@result$Description)
ego_BP@result$Description <- stringr::str_wrap(ego_BP@result$Description, width = 60)

ego_BP1 <- enrichGO(gene         = gene,
                    OrgDb         = orgdb,
                    keyType       = 'SYMBOL',
                    ont           = "BP",
                    pAdjustMethod = "BH",
                    pvalueCutoff  = 1,
                    qvalueCutoff  = 1)


outSVG_BP <- paste("BP_",name,".svg", sep="")

g_BP<-dotplot(ego_BP,showCategory=20)+geom_point(shape = 21, aes(fill = p.adjust), alpha = 0.9, stroke = 0.5) + 
  scale_size_continuous(range = c(3, 10)) +  # Adjust size range
  scale_color_gradient(low = "blue", high = "red") +  # Blue to red color gradient for CorrectedP.Value
  scale_fill_gradient(low = "blue", high = "red") +  # Blue to red color gradient for fill
  theme_minimal() +  # Minimal theme for clarity
  labs(x = "GeneRatios", y = "") +  # Label the axes
  # ggtitle("Dot Plot of GeneRatio vs Term") +  # Add title
  theme(panel.grid.major.y = element_line(linewidth = .5, linetype = "solid", colour = 'gray'),
        panel.grid.major.x = element_line(linewidth = .5, linetype = "solid", colour = 'gray'),
        panel.border = element_rect(color = "black", fill = NA, linewidth = 1),  # Add outer border
        axis.text.y = element_text(size = 10, face = "bold"))


outPNG_BP <- paste("BP_", name, ".png", sep = "")
outPDF_BP <- paste("BP_", name, ".pdf", sep = "")

if(g_BP$data |> nrow() >=1){
  ggsave(outPNG_BP, plot = g_BP, width = 8, height = 6, dpi = 300)
  ggsave(outPDF_BP, plot = g_BP, width = 8, height = 6)
  ggsave(outSVG_BP, plot = g_BP)
}



outTable <- paste("BP_",name,".GO.xls", sep="")
write.table(ego_BP1, file=outTable, sep="\t", eol="\n",quote=FALSE, 
            row.names = FALSE)

ego_MF <- enrichGO(gene         = gene,
                   OrgDb         = orgdb,
                   keyType       = 'SYMBOL',
                   ont           = "MF",
                   pAdjustMethod = "BH",
                   pvalueCutoff  = 0.01,
                   qvalueCutoff  = 0.05)

# Shorten the term names on the y-axis and make Title case
ego_MF@result$Description <- substr(ego_MF@result$Description, 1, 30) # Adjust width as needed
ego_MF@result$Description <- tools::toTitleCase(ego_MF@result$Description)
ego_MF@result$Description <- stringr::str_wrap(ego_MF@result$Description, width = 60)
ego_MF1 <- enrichGO(gene         = gene,
                    OrgDb         = orgdb,
                    keyType       = 'SYMBOL',
                    ont           = "MF",
                    pAdjustMethod = "BH",
                    pvalueCutoff  = 1,
                    qvalueCutoff  = 1)


outSVG_MF <- paste("MF_",name,".svg", sep="")

g_MF <- dotplot(ego_MF,showCategory=20) +geom_point(shape = 21, aes(fill = p.adjust), alpha = 0.9, stroke = 0.5) + 
  scale_size_continuous(range = c(3, 10)) +  # Adjust size range
  scale_color_gradient(low = "blue", high = "red") +  # Blue to red color gradient for CorrectedP.Value
  scale_fill_gradient(low = "blue", high = "red") +  # Blue to red color gradient for fill
  theme_minimal() +  # Minimal theme for clarity
  labs(x = "GeneRatios", y = "") +  # Label the axes
  # ggtitle("Dot Plot of GeneRatio vs Term") +  # Add title
  theme(panel.grid.major.y = element_line(linewidth = .5, linetype = "solid", colour = 'gray'),
        panel.grid.major.x = element_line(linewidth = .5, linetype = "solid", colour = 'gray'),
        panel.border = element_rect(color = "black", fill = NA, linewidth = 1),  # Add outer border
        axis.text.y = element_text(size = 10, face = "bold"))

outPNG_MF <- paste("MF_", name, ".png", sep = "")
outPDF_MF <- paste("MF_", name, ".pdf", sep = "")

if(g_MF$data |> nrow() >=1){
  ggsave(outPNG_MF, plot = g_MF, width = 8, height = 6, dpi = 300)
  ggsave(outPDF_MF, plot = g_MF, width = 8, height = 6)
  ggsave(outSVG_MF, plot = g_MF)
}


outTable <- paste("MF_",name,".GO.xls", sep="")
write.table(ego_MF1, file=outTable, sep="\t", eol="\n",quote=FALSE, 
            row.names = FALSE)


