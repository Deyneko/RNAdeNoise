# RNAdeNoise - an R function for cleaning RNA-seq data from low-reads genes

In this work, we developed a method for cleaning RNA-seq data, that improves the detection of genes with differential expression (DEGs), and specifically the  genes with low to moderate transcription (genes with high transcription statistically easy to detect). It is assumed that the measured mRNA counts consist of real and random parts. Using a data modeling approach, parameters of randomly distributed mRNA counts are identified and a number of mRNA reads, most probably originating from a technical noise, is determined individually for each dataset.

We demonstrate that the removal of this random component leads to the detection of more genes (DEGs), more significant p-values, and less statistical artifacts compared to the use of raw data or the use of filters based on fixed thresholds (for example, removal of genes with counts ≤ 3). Using our RNA-seq data on polysome profiling on Arabidopsis thaliana, a significant increase in the number of detected differentially translated regulatory genes was shown. Additionally, the method was applied to several published RNA-seq datasets covering different sequencing technologies and organisms, and in all cases, a significant increase in detected differentially expressed genes was shown. The program substitutes the widely used fixed threshold approach to remove low level mRNAs.

Please see our publication for more details: https://pubmed.ncbi.nlm.nih.gov/36384457/
