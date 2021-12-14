# ************** cleaning as a function *****************************
# input format - dataframe with columns: GeneId, Value, Value, ....
# output - dataframe in the same format
# 
# Here you can find the function itself and further examples how to use.
# *******************************************************************


RNAcountsDeNoise=function (RNAcounts, CleanStrength=0.90, RemoveZeroGenes=0) {
# INPUT: 
#  RNAcounts - dataframe of raw count values (for example output from featurecount etc.)
#  CleanStrength - exponential part to remove 
#
#  Some experimental parameters:
#  RemoveZeroGenes=0 - if set to 1 than genes with all zeroes will be removed. Genes number wiill change!
   AbsoluteNum=0; #CleanStrength can represent either a quantile of exponential part or threshould value for the exponent. 
#  When the exponent drops under this value the corresponding X (number of counts) is set to be the threshold for number of counts to remove.
#  Usage example: AbsoluteNum=1; CleanStrength=3; - the program will find x so that A*exp(-Bx)<=3. 

  
  RNAcounts1=RNAcounts;
  
  n=ncol(RNAcounts);#print (n)
  #loop on every column 
  for (L in 1:n) {
    x=RNAcounts[,L]
    x1=length(subset(x,x[]==1));x2=length(subset(x,x[]==2));x3=length(subset(x,x[]==3));x4=length(subset(x,x[]==4));
    
    z=data.frame(c(1,2,3,4),c(x1,x2,x3, x4))
    colnames(z)=c("x","y")
    T=lm(log(y)~x,z);#T
    x1=exp(T$coefficients[1]);x2=T$coefficients[2]; #cat ("Coeff=",x1,x2,"   "); # for debug only
    Thresh=-1; sum=0; TotalExp=-x1*exp(x2)/x2; 
    for (j in 1:100){
      if (AbsoluteNum == 1) {
        if (x1*exp(x2*j) <= CleanStrength ) {Thresh=j; break; }}
      else {
        sum=sum+x1*exp(x2*j); #print (sum)
        if (x1*exp(x2*j) <= (1-CleanStrength)*TotalExp)  {Thresh=j;break; }}
      
    }
    if (Thresh == -1) {print ("Increase value in loop or check your data. Noise level above 100.") }
    print (c("Subtraction value = ", Thresh))
    
    
    # De-noise data
    
    RNAcounts1[,L]=RNAcounts[,L]-Thresh;RNAcounts1[RNAcounts1[,L]<0,L]=0;
    
    # next column
  }
  # remove genes with all zeroes 
  
  if (RemoveZeroGenes) { 
    RNAcounts1$s=rowSums(RNAcounts1[,1:n])
    RNAcounts1=subset(RNAcounts1,RNAcounts1$s>0)
    RNAcounts1=RNAcounts1[,c(-(n+1))] }
  
  #Dat1$s=rowSums(Dat1[,1:6])
  #Dat1=subset(Dat1,Dat1$s>=10)
  #Dat1=Dat1[,-c(7)]
  
  return (RNAcounts1);
}




# Examples how to use RNAdeNoise

library(readxl)
library("DESeq2", lib.loc="~/R/win-library/4.0")
library("edgeR", lib.loc="C:/Program Files/R/R-4.0.3/library")
library(dplyr)

# *****************************************
# our dataset: polysomal and monosomal mRNA in Arabidopsis
# *****************************************

# input raw counts --------  our data
counts <- read.delim("countsPolyMono.txt") # our data 
counts <- counts[,c(1,7,8,9,10,11,12)];row.names(counts)=counts[,1]
Condition<-c("M", "M", "M","P", "P", "P")
Dat1=RNAcountsDeNoise(counts[,-1], CleanStrength=0.9); 
#Dat1=counts[,-1]; # to use raw data
#write.table(Dat1,file="D:/Igor/programs/TranslationEff/TranscriptomeAssembly/CountsCleaned.txt", quote = FALSE, sep = "\t", row.names = FALSE)

# DEGs detection 
# EdgeR 
y<-DGEList(counts=Dat1, group=Condition, genes=rownames(Dat1))
y<-calcNormFactors(y, method="TMM");y$samples
y<-estimateCommonDisp(y)
y<-estimateTagwiseDisp(y)
et<-exactTest(y, pair=c("M", "P"))
res<-as.data.frame(topTags(et, n=37336))
res=subset(res,(res$PValue<0.0001)&(abs(res$logFC)>=1.5)); nrow(res);
#write.table(res$genes,file="D:/Igor/programs/TranslationEff/TranscriptomeAssembly/DEG_Raw.txt", quote = FALSE, sep = "\t", row.names = FALSE,col.names = FALSE)


#DESeq2
colData=as.data.frame(Condition); 
names (colData) <-c("condition"); 
colData$condition <- factor(colData$condition);
ds <- DESeqDataSetFromMatrix(countData = Dat1, colData = colData, design = ~ condition); 
ds <- DESeq(ds); 
res <- results(ds)
res=subset(res,(res$pvalue<0.0001)&(abs(res$log2FoldChange)>=1.5)); nrow(res); 
#write.table(row.names(res),file="D:/Igor/programs/TranslationEff/TranscriptomeAssembly/DEG_DESEq_DENOISE001_pval01.txt", quote = FALSE, sep = "\t", row.names = FALSE,col.names = FALSE)



# *****************************************
# Other Datasets 1: response to Cycloastragenol, BGISEQ-500 platform, Arabidopsis. 
# Mhiri, W.et al., 2020. Transcriptomic analysis reveals responses to Cycloastragenol in Arabidopsis thaliana. PLoS One.
# *****************************************

counts <- read.delim("countsBGI15bp.txt") # BGI 
counts <- counts[,c(1,7,8)]; row.names(counts)=counts[,1]
Condition<-c("M","P")
Dat1=RNAcountsDeNoise(counts[,-1], CleanStrength=0.9);  # Cleaning !!! run only one line !!
Dat1=counts[,-1]; #Raw data 
Dat1[Dat1[,]<=3]=0 # filter counts > 3


y<-DGEList(counts=Dat1, group=Condition, genes=row.names(Dat1))
y<-calcNormFactors(y, method="TMM");y$samples
et<-exactTest(y, dispersion=0.1^2)
res<-as.data.frame(topTags(et, n=38186))
nrow(subset(res,(res$PValue<=0.002)&(res$logFC>=1.0)));nrow(subset(res,(res$PValue<=0.002)&(res$logFC<=-1.0)))

#Generate gene Datasets
res=subset(res,(res$PValue<=0.002)&(abs(res$logFC)>=1.0));nrow(res) # author= p=0.0019, FC>=1.0
#write.table(row.names(res),file="D:/Igor/programs/TranslationEff/TranscriptomeAssembly/BGI_DEG_HTS.txt", quote = FALSE, sep = "\t", row.names = FALSE,col.names = FALSE)
write.table(res, "BGI_expr.txt", append = FALSE, sep = " ", dec = ".",row.names = FALSE, col.names = TRUE)


# *****************************************
# Other Datasets 2: a circadian clock in Arabidopsis thaliana, Illumina platform, Arabidopsis. 
# Bonnot, T. and Nagel, D.H. (2021) Time of the day prioritizes the pool of translating mRNAs in response to heat stress. Plant Cell, 33, 2164-2182.
# *****************************************

counts <- read.delim("GSE158444_Counts_Time.txt")
#counts <- counts[,c("AGI","TOT.T0.C1","TOT.T0.H1","TR.T0.C1","TR.T0.H1")]; row.names(counts)=counts[,1]
counts <- counts[,c("AGI","TOT.T0.C1","TOT.T0.C2","TOT.T0.C3","TOT.T0.H1","TOT.T0.H2","TOT.T0.H3")]; row.names(counts)=counts[,1]

Condition<-c("C", "C","C","H","H","H")
Dat1=RNAcountsDeNoise(counts[,-1], CleanStrength=0.9);  # Cleaning !!! run only one line !!
Dat1=counts[,-1]; #Raw data 
Dat1[Dat1[,]<=3]=0 # filter counts > 3
Dat1$s=rowSums(Dat1[,1:6]);Dat1=subset(Dat1,Dat1$s>20);Dat1=Dat1[,-c(7)]; # Filter sum>20 used by authors
#Dat1$s=rowSums(Dat1[,4:6]);Dat1=subset(Dat1,Dat1$s>20);Dat1=Dat1[,-c(7)]

#DESEQ2
colData=as.data.frame(Condition); names (colData) <-c("condition"); 
colData$condition <- factor(colData$condition);
ds <- DESeqDataSetFromMatrix(countData = Dat1, colData = colData, design = ~ condition); 
ds <- DESeq(ds);
res <- results(ds) 
nrow(subset(res,(res$pvalue<0.0001)&(abs(res$log2FoldChange)>=1.5))); nrow(res);#res

write.table(res, "circad.txt", append = FALSE, sep = " ", dec = ".",row.names = TRUE, col.names = TRUE)




# *****************************************
# Other Datasets 3: mRNA between Alport mice and wild mice, Illumina platform, Mouse 
# Dufek, B., et al., (2020) RNA-seq analysis of gene expression profiles in isolated stria vascularis from wild-type and Alport mice reveals key pathways underling Alport strial pathogenesis. PLoS One, 15, e0237907.
# *****************************************

library(tidyverse);  
# These gene names needed to annotate the resulted table with EnsebleID, which in turn is used by DAVID for functional annotation.
GeneNames <- read.delim("mart_export(10).txt"); colnames(GeneNames)=c("EnsID","EntrID")
GeneNames$EntrID=as.character(GeneNames$EntrID)

counts <- read.delim("mouse_counts.txt")
counts <- counts[,-2]; row.names(counts)=counts[,1]
Condition<-c("M", "P")

Dat1=RNAcountsDeNoise(counts[,-1], CleanStrength=0.9);  # Cleaning !!! run only one line !!
Dat1=counts[,-1]; #Raw data 
Dat1[Dat1[,]<=3]=0 # filter counts > 3

y<-DGEList(counts=Dat1, group=Condition, genes=row.names(Dat1))
y<-calcNormFactors(y, method="TMM");y$samples
et<-exactTest(y, dispersion=0.1^2)
res<-as.data.frame(topTags(et, n=27179))
nrow(subset(res,(res$PValue<=0.0001)&(abs(res$logFC)>=1.0)));

res=subset(res,(res$PValue<=0.0001)&(abs(res$logFC)>=1.0)); nrow(res);
Dat3=left_join(res,GeneNames,by=c("genes"="EntrID"))
#counts2=counts2[!duplicated(counts2$EnsID),]; counts=counts2[,c(4,2,3)]; row.names(counts)=counts[,1]
write.table(Dat3, "mouse_DEG_DeNoise.txt", append = FALSE, sep = " ", dec = ".",row.names = FALSE, col.names = TRUE)

