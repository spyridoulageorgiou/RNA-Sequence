#!/bin/bash

#SBATCH --mail-user=spyridoulageorgiou@gmail.com
#SBATCH --job-name="counting"
#SBATCH --nodes=1
#SBATCH --cpus-per-task=4
#SBATCH --time=30:00:00
#SBATCH --mem=25G


module load UHTS/Analysis/HTSeq/0.9.1

for FILE in /data/courses/rnaseq_course/lncRNAs/Project2/users/SGeorgiou/Data/bam_files_reverse/*.bam
do 
SAMPLE=`basename $FILE | sed s/\.bam//` 
htseq-count -f bam $FILE /data/courses/rnaseq_course/lncRNAs/Project2/users/SGeorgiou/Data/gtf_files_reverse/merged.gtf --idattr gene_id --additional-attr=gene_name -m union --stranded=no --type exon > $SAMPLE".txt" 
done

