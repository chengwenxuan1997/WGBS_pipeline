# WGBS_pipeline
A pipeline to process the .bismark.cov.gz file  
Follow the instruction of [GeneDMRs](https://github.com/xiaowangCN/GeneDMRs) and CompEpigen's [WGBS best practices workflow](https://compepigen.github.io/methrix_docs/articles/methrix.html).  

## Purpose 
1: Extract methylation or coverage matrices from .bismark.cov.gz files for downstream analysis.  
2: Call Differential Methylation Regions (DMRs) using [metilene](https://www.bioinf.uni-leipzig.de/Software/metilene/), [dmrseq](https://www.bioconductor.org/packages/release/bioc/html/dmrseq.html), and GeneDMRs.  
3: Annotate the DMRs using ChIPseeker  

## Notes  
1: The **WGBS_bsseq-04-27.R** code will install metilene from the compressed source package stored under **./Packages automatically**. Because the installation path of metilene was not added into the $PATH, installed metilene will not be called in the command line directly (a full path to metilene is needed).  
2: It has not been supported to use metilene in the **windows system** yet, because the metilene.exe has not been compiled successfully.  
3: It is not practical to run this code on a personal computer because of the great cost on memory.
