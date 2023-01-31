# RNA-Sequence
There are three main objectives for this project. 

The first one is to annotate the protein-coding and the long noncoding genes in the cells. The first goal was met by using htseq-count which gets a BAM file and counts for each gene how many aligned reads overlap its exons.

The second objective is to identify the differentially expressed genes in the cell population. This was achieved by using DESeq which detects count outliers and removes those genes from the analysis. It also removes the genes whose mean of normalized counts is below a certain threshold that is determined by an optimization procedure.

The last one is to take the differentially expressed genes and check how they behave between the populations.
