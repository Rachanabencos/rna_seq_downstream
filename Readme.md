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