### Generate a final table containing all the analysis
DEGs <- read.delim("transcripts.DEGs.txt")
DEGs <- DEGs[,c(1,3,6)]
colnames(DEGs)[1:3] <- c("Gene_id","log2FC","padj")

DEGs$type <- ifelse(grepl("MSTRG", DEGs$Gene_id), "novel", "annotated")


####### Get single exons from gtf file


gtf <- read.table("/data/courses/rnaseq_course/lncRNAs/Project2/users/SGeorgiou/Data/gtf_files_reverse/merged.gtf", header = F, sep = "\t")
colnames(gtf) <- c("chromosome", "source", "feature", "start", "end", "score", "strand", "phase", "attributes")

# split the 9th column into a list of key-value pairs
split_9th_column <- strsplit(as.character(gtf[,9]), ";")

# extract the keys and values from the list
transcript_ids <- sapply(split_9th_column, function(x) gsub(" ", "", x[2]))

# create new columns in the dataframe for the keys and values
gtf$transcript_ids <- transcript_ids

gtf$transcript_id <- gsub("transcript_id", "", gtf$transcript_id)


exons <- subset(gtf, feature == "exon")


# group exons by transcript_id
exon_counts <- table(exons$transcript_id)

# select only transcripts with a count of 1 exon
single_exon_transcripts <- names(exon_counts)[exon_counts == 1]


####################### Find single exons from DEGs ##############
DEGs$single_exon <- ifelse(DEGs$Gene_id %in% single_exon_transcripts, "yes", "no")

##### Get information for FANTOM #######

FANTOM <- read.delim("/data/courses/rnaseq_course/lncRNAs/Project2/users/SGeorgiou/Data/step6/FANTOM/FANTOM_overlap.txt", header = FALSE)
polyA <- read.delim("/data/courses/rnaseq_course/lncRNAs/Project2/users/SGeorgiou/Data/step6/FANTOM/polyA_tail.overlap.txt", header = FALSE)

# split the 9th column into a list of key-value pairs
split_5th_column <- strsplit(as.character(FANTOM[,5]), ";")

# extract the keys and values from the list
transcript_ids <- sapply(split_5th_column, function(x) gsub(" ", "", x[2]))

# create new columns in the dataframe for the keys and values
FANTOM$transcript_ids <- transcript_ids

FANTOM$transcript_id <- gsub("transcript_id", "", FANTOM$transcript_ids)

####################### Find FANTOM overlap ##############
DEGs$FANTOM <- ifelse(DEGs$Gene_id %in% FANTOM$transcript_id, "yes", "no")

# split the 9th column into a list of key-value pairs
split_5th_column <- strsplit(as.character(polyA[,5]), ";")

# extract the keys and values from the list
transcript_ids <- sapply(split_5th_column, function(x) gsub(" ", "", x[2]))

# create new columns in the dataframe for the keys and values
polyA$transcript_ids <- transcript_ids

polyA$transcript_id <- gsub("transcript_id", "", polyA$transcript_ids)

####################### Find polyA overlap ##############
DEGs$polyA <- ifelse(DEGs$Gene_id %in% polyA$transcript_id, "yes", "no")


############################ CPAT analysis ################################
cpat <- read.delim("/data/courses/rnaseq_course/lncRNAs/Project2/users/SGeorgiou/Data/step6/cpat/protein_coding_potential")


### Get a coordinate name in gtf
trans <- subset(gtf, feature == "transcript")
trans$chromosome.new<-toupper(trans$chromosome)
trans$id <- paste(trans$chromosome.new,":", trans$start,"-", trans$end, sep = "")

##### Merge the 2 tables

cpat.trascripts <- merge(cpat, trans[,c("transcript_id", "id")], by.x= "coord", by.y = "id")

#Keep only transcripts with protein coding potential
pr.coding.potential <- cpat.trascripts[cpat.trascripts$coding_prob > 0.364,]

####################### Find Whether there is protein coding potential ##############
DEGs$cpat.coding.potential<- ifelse(DEGs$Gene_id %in% pr.coding.potential$transcript_id, "yes", "no")

