<b>Here are gene lists, which are the output from filtering and DEG detection (using EdgeR and DESeq2).</b>

<p>Datasets = directory names
<br>file names = Dataset_[program]_Filter.txt

In details:
<p><b>ArabidopsisCircadian</b> - datasets from Bonnot, T. and Nagel, D.H. (2021) Time of the day prioritizes the pool of translating mRNAs in response to heat stress. Plant Cell, 33, 2164-2182.
<p>Here we used DESeq2 program as is in the original publication.
<br>USA_DeNoise.txt - filtred using RNAdeNoise
<br>USA_FIX10.txt - filtred using fixed threshold filterr counts>10
<br>USA_FIX3.txt - filtred using counts>3
<br>USA_FIX5.txt - filtred using counts>5
<br>USA_FPKM.txt - filtred using FPKM>0.3 (frequency per nucleotide per kilobase per million reads)
<br>USA_HTS.txt - filtred using HTSFilter 
<br>USA_RAW.txt - no filtering 
<br>USA_SUM20.txt - filtred using Authors' filter - sum in all samples must be above 20

<p><b>BGI_Mhiri</b> - datasets from Mhiri, W., Ceylan, M., Turgut-Kara, N., Nalbantolu, B....  (2020) Transcriptomic analysis reveals responses to Cycloastragenol in Arabidopsis thaliana. PLoS One, 15, e0242986
<br>File naming: same as above, but we used EdgeR program 

<b>Mouse</b> - datasets from Dufek, B., Meehan, D.T., Delimont, D., Wilhelm, K., Samuelson, G., Coenen, R., Madison, J., Doyle, E., Smyth, B., Phillips, G. et al. (2020) RNA-seq analysis of gene expression profiles in isolated stria vascularis from wild-type and Alport mice reveals key pathways underling Alport strial pathogenesis. PLoS One, 15, e0237907
<br>File naming: same as above, we used EdgeR program 

<b>PolysomalMonosomal</b> - our dataset on Arabidopsis [reference folowing]
<br>simialr to above, except: DEGs using DESeq2 program: 
<br>DEG_DESq_DeNoise01.txt
<br>DEG_DESq_FIX10.txt
<br>DEG_DESq_FIX3.txt
<br>DEG_DESq_FIX5.txt
<br>DEG_DESq_FPKM.txt
<br>DEG_DESq_HTS.txt
<br>DEG_DESq_Raw.txt

<p>DEGs using EdgeR program: 
<br>DEG_DeNoise01_All.txt
<br>DEG_FIX10_All.txt
<br>DEG_FIX3_All.txt
<br>DEG_FIX5_All.txt
<br>DEG_FPKM_All.txt
<br>DEG_HTS_All.txt
<br>DEG_RAW_All.txt

