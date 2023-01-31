
#Load biomart
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("biomaRt")

library("biomaRt")
listMarts()

ensembl <- useMart("ensembl")
ensembl = useDataset("hsapiens_gene_ensembl",mart=ensembl)
attributes = listAttributes(ensembl)


DEGs <- read.table("C:/Users/Spyridoula/OneDrive/Desktop/counts reverse/Deliverables/Step 5/transcripts.DEGs.txt", header = TRUE)

biomart_table <- getBM(attributes= c("ensembl_gene_id_version", 
                                     "ensembl_transcript_id_version",
                                     "chromosome_name",
                                     "transcript_start",
                                      "transcript_end",
                                     "strand",
                                     "gene_biotype",
                                     "transcript_biotype"), 
                       mart = ensembl)

DEGs <- merge(DEGs, biomart_table, by.x = "row", by.y = "ensembl_transcript_id_version")



#Reorder
DEGs <- DEGs[,c(9:14,1:8)]

#Add chr in first column
DEGs$chromosome_name <- paste0("chr", substr(DEGs$chromosome_name, 1, nchar(DEGs$chromosome_name)))

write.table(DEGs, "C:/Users/Spyridoula/OneDrive/Desktop/counts reverse/transcripts.annotation.txt", sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)

#### Do the same to extract gene information only

biomart_table <- getBM(attributes= c("ensembl_gene_id_version", 
                                     "chromosome_name",
                                     "start_position",
                                     "end_position",
                                     "strand",
                                     "gene_biotype",
                                     "transcript_biotype"), 
                       mart = ensembl)

hg38 <- biomart_table[,c(2:5,1,6:7)]
#Add chr in first column
hg38$chromosome_name <- paste0("chr", substr(hg38$chromosome_name, 1, nchar(hg38$chromosome_name)))
write.table(hg38, "C:/Users/Spyridoula/OneDrive/Desktop/counts reverse/hg38.annotation.txt", sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)
