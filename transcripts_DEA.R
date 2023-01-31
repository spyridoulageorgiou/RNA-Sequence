P1_L3_R2.transcript <- read.delim("C:/Users/Spyridoula/OneDrive/Desktop/counts reverse/P1_L3_R2-transcript.txt", header=FALSE)
P2_L3_R2.transcript <- read.delim("C:/Users/Spyridoula/OneDrive/Desktop/counts reverse/P2_L3_R2-transcript.txt", header=FALSE)
P3_L3_R2.transcript <- read.delim("C:/Users/Spyridoula/OneDrive/Desktop/counts reverse/P3_L3_R2-transcript.txt", header=FALSE)
`3_2_L3.transcript` <- read.delim("C:/Users/Spyridoula/OneDrive/Desktop/counts reverse/3_2_L3-transcript.txt", header=FALSE)
`3_4_L3.transcript` <- read.delim("C:/Users/Spyridoula/OneDrive/Desktop/counts reverse/3_4_L3-transcript.txt", header=FALSE)
`3_7_L3.transcript` <- read.delim("C:/Users/Spyridoula/OneDrive/Desktop/counts reverse/3_7_L3-transcript.txt", header=FALSE)

transcripts.counts.table <- cbind(`3_2_L3.transcript`, `3_4_L3.transcript`[,3], `3_7_L3.transcript`[,3],P1_L3_R2.transcript[,3],P2_L3_R2.transcript[,3],P3_L3_R2.transcript[,3])
colnames(transcripts.counts.table) <- c("gene_id","gene_name","3_2_L3", "3_4_L3","3_7_L3", "P1_L3", "P2_L3", "P3_L3")
write.table(transcripts.counts.table, "C:/Users/Spyridoula/OneDrive/Desktop/counts reverse/transcripts.counts.txt", sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE)



#How many genes did you detect?
nrow(transcripts.counts.table[rowSums(transcripts.counts.table[,3:8]) > 0,])


##How many novel transcripts and genes did you detect? MSTRG approach
novel.genes_mstrg <- transcripts.counts.table[substr(transcripts.counts.table[,1],1,5) == "MSTRG",]
nrow(novel.genes_mstrg[rowSums(novel.genes_mstrg[,3:8]) > 0,])

###################### Step 5  ######################
library(DESeq2)
countData<- transcripts.counts.table[,c(3:8)]
rownames(countData) <- transcripts.counts.table[,1]
condition <- factor(c("paraclone","paraclone", "paraclone", "parental","parental","parental"))

dds <- DESeqDataSetFromMatrix(countData, DataFrame(condition), ~ condition)

dds <- DESeq(dds)
#Extract dataframe
res <- results(dds, contrast=c("condition","paraclone","parental"))
df <- results(dds, contrast=c("condition","paraclone","parental") ,tidy=TRUE)

DEGs.transcripts <- subset(df, padj<.01 & abs(log2FoldChange)> 1)
write.table(DEGs.transcripts, "C:/Users/Spyridoula/OneDrive/Desktop/counts reverse/transcripts.DEGs.txt", sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE)


##### Volcano plot

library(ggplot2)
#reset par
par(mfrow=c(1,1))
# Make a basic volcano plot
with(res, plot(log2FoldChange, -log10(pvalue), pch=20, main="Volcano plot", xlim=c(-10,10)))

# Add colored points: blue if padj<0.01, red if log2FC>1 and padj<0.05)
with(subset(res, padj<.01 ), points(log2FoldChange, -log10(pvalue), pch=20, col="blue"))
with(subset(res, padj<.01 & abs(log2FoldChange)> 1), points(log2FoldChange, -log10(pvalue), pch=20, col="red"))


#####Get top 100  up and down
DEGs.ordered <- DEGs.transcripts[order(DEGs.transcripts$log2FoldChange,decreasing = TRUE),]
top100 <- rbind( head(DEGs.ordered, 100), tail(DEGs.ordered, 100))

#Heatmap Matrix
library(ComplexHeatmap)
library(circlize)
#Extract normalized counts table
normalized.counts <- counts(dds, normalized=T)
normalized.counts <- cbind(normalized.counts,as.data.frame(rownames(normalized.counts)))
colnames(normalized.counts)[7] <- "gene_id"

counts.table.heatmap <- merge(normalized.counts, top100, by.x = "gene_id", by.y = "row")

###Prepare table for heatmaps 
counts <- counts.table.heatmap[,c(2:7)]
z_list<-apply(counts,1,scale)
C <-as.data.frame(z_list)
final_list <- t(C)
colnames(final_list) <- c("3_2_L3", "3_4_L3","3_7_L3", "P1_L3", "P2_L3", "P3_L3")
rownames(final_list) <- counts.table.heatmap$gene_id

Heatmap(final_list, 
        name = "z-score", 
        column_title = "Samples", 
        row_title = "Top 100 differentially expressed Transcripts", 
        cluster_rows = TRUE,
        show_row_names = FALSE
)



##### PCA

vsdata <- vst(dds, blind=FALSE)

plotPCA(vsdata)
