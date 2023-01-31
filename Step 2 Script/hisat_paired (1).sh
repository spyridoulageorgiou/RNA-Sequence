#Step 2   #  /data/courses/rnaseq_course/lncRNAs/Project2/users/SGeorgiou/Data/bam_files_reverse/  #

module load UHTS/Aligner/hisat/2.2.1
source /data/users/lfalquet/SBC07107_22/scripts/module.sh

HISAT2_INDEXES=/data/courses/rnaseq_course/lncRNAs/Project2/users/SGeorgiou/hisat2-2.2.1/grch38/genome
FASTQ=/data/courses/rnaseq_course/lncRNAs/Project2/users/SGeorgiou/Data/fastq/fastq_that_we_need/1_1_L3_R2_001_qyjToP2TB6N7.fastq.gz
CORES=8

hisat2 -p $CORES --add-chrname -x $HISAT2_INDEXES -1 3_2_L3_R1_001_DID218YBevN6.fastq.gz -2 3_2_L3_R2_001_UPhWv8AgN1X1.fastq.gz -S 3_2_L3.sam

samtools sort -@ $CORES 3_2_L3.sam -o 3_2_L3.bam


#NEW
hisat2 -p $CORES --add-chrname -x $HISAT_INDEXES -1 3_2_L3_R2_001_UPhWv8AgN1X1.fastq.gz -2 3_2_L3_R1_001_DID218YBevN6.fastq.gz -S 3_2_L3.sam

samtools sort -@ $CORES 3_2_L3.sam - o 3_2_L2.bam


hisat2 -p $CORES --add-chrname -x $HISAT2_INDEXES -1 3_4_L3_R1_001_QDBZnz0vm8Gd.fastq.gz -2 3_4_L3_R2_001_ng3ASMYgDCPQ.fastq.gz -S 3_4_L3.sam

samtools sort -@ $CORES 3_4_L3.sam -o 3_4_L3.bam



hisat2 -p $CORES --add-chrname -x $HISAT2_INDEXES -1 3_7_L3_R1_001_Tjox96UQtyIc.fastq.gz -2  3_7_L3_R2_001_f60CeSASEcgH.fastq.gz -S 3_7_L3.sam

samtools sort -@ $CORES 3_7_L3.sam -o 3_7_L3.bam


###Parental

hisat2 -p $CORES --add-chrname -x $HISAT2_INDEXES -1 P1_L3_R1_001_9L0tZ86sF4p8.fastq.gz -2  P1_L3_R2_001_yd9NfV9WdvvL.fastq.gz -S P1_L3_R2.sam

samtools sort -@ $CORES P1_L3_R2.sam -o P1_L3_R2.bam


hisat2 -p $CORES --add-chrname -x $HISAT2_INDEXES -1 P2_L3_R1_001_R82RphLQ2938.fastq.gz -2  P2_L3_R2_001_06FRMIIGwpH6.fastq.gz -S P2_L3_R2.sam

samtools sort -@ $CORES P2_L3_R2.sam -o P2_L3_R2.bam


hisat2 -p $CORES --add-chrname -x $HISAT2_INDEXES -1 P3_L3_R1_001_fjv6hlbFgCST.fastq.gz -2  P3_L3_R2_001_xo7RBLLYYqeu.fastq.gz -S P3_L3_R2.sam

samtools sort -@ $CORES P3_L3_R2.sam -o P3_L3_R2.bam



#####FOR BAM BAI FILES
samtools index file.bam file.bam.bai
