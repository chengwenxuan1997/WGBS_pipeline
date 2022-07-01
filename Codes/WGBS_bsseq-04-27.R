
Sys.setenv("R_REMOTES_NO_ERRORS_FROM_WARNINGS" = "true")
workdir <- "/home/a208/data/cwx/Xlu/Amanda_WGSmb"; setwd(workdir)
fig.path <- file.path(workdir,"Figures")
res.path <- file.path(workdir,"Results")
data.path <- file.path(workdir,"InputData")
comAnn.path <- file.path(workdir,"Annotation") ## path containning annotation files created by Jack
WGSmb.path <- file.path(data.path,"Amanda_WGSmb_Project_Feb2022")
pkg.path <- file.path(workdir, "Packages")

pkg.tmp.path <- file.path(pkg.path, "metilene_input")
res.tmp.path <- file.path(res.path, "bsseq") ## result provided by bsseq.

invisible(lapply(ls()[grep("path", ls())], function(x){
  if (!dir.exists(get(x))) dir.create(get(x))
}))

# if (!requireNamespace("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")
# 
# #Installing stable version from BioConductor
# #BiocManager::install("methrix")
# BiocManager::install("annotatr")
# BiocManager::install("TxDb.Mmusculus.UCSC.mm10.knownGene")
# 
# # #Installing developmental version from GitHub
# devtools::install_github("CompEpigen/methrix")
# BiocManager::install("ChIPseeker")
# BiocManager::install("TxDb.Hsapiens.UCSC.hg19.knownGene")
# BiocManager::install("TxDb.Mmusculus.UCSC.mm10.knownGene")
# BiocManager::install("dmrseq")



# Read and Calculate ------------------------------------------------------


# load R packages
library(methrix)
library(BSgenome.Hsapiens.UCSC.hg19)
library(BSgenome.Mmusculus.UCSC.mm10)
library(BiocParallel)
library(bsseq)
library(dmrseq)
options(MulticoreParam=MulticoreParam(workers = 10))


# extract bedgraph files
bmc_files <- ".*.bismark.cov.gz"
bmc_files <- dir(WGSmb.path, pattern = bmc_files)
print(basename(bmc_files))

# load sample annotation
sample_anno <- read.table(file.path(data.path,"WGSmb_sampleInfo_Feb2022.txt"),sep = "\t",row.names = NULL,check.names = F,stringsAsFactors = F,header = T)
rownames(sample_anno) <- paste0(sample_anno$sampleID, "_1_val_1_bismark_bt2_pe.deduplicated")

# First extract genome wide CpGs from the desired reference genome
mm10_cpgs <- suppressWarnings(methrix::extract_CPGs(ref_genome = "BSgenome.Mmusculus.UCSC.mm10"))

#Read the files 
meth <- methrix::read_bedgraphs(
  files = file.path(WGSmb.path,bmc_files),
  ref_cpgs = mm10_cpgs,
  pipeline = "Bismark_cov",
  # chr_idx = 1,
  # start_idx = 2,
  # end_idx = 3,
  # M_idx = 5,
  # U_idx = 6,
  # beta_idx = 4,
  #stranded = FALSE,
  zero_based = FALSE,
  #collapse_strands = FALSE, 
  coldata = sample_anno
)

# save as backup
meth.bk <- meth
saveRDS(meth.bk, file.path(res.tmp.path,"unfiltered_mm10_meth.rds"))


# filter methylation site -------------------------------------------------

cat(as.character(Sys.time()), ": filter the methylation matrix\n")
meth <- readRDS(file = file.path(res.tmp.path,"unfiltered_mm10_meth.rds"))

# Remove uncovered loci
meth <- remove_uncovered(m = meth) # Removed 974,604 [4.46%] uncovered loci of 21,867,837 sites

# Remove SNPs (Only hg19 and hg38 genomes are currently supported)
# meth <- methrix::remove_snps(m = meth)

# Remove chrM
# table(meth@elementMetadata$chr)
# chr.to.keep <- names(table(meth@elementMetadata$chr))[table(meth@elementMetadata$chr)>100]
# meth <- methrix::subset_methrix(m = meth, contigs = chr.to.keep) 
# plot_stats(plot_dat = get_stats(m = meth), what = "C", stat = "mean")
saveRDS(meth, file.path(res.tmp.path, "uncovered_removed_meth.rds"))

meth.mtx <- methrix::get_matrix(m = meth, type = "M", add_loci = T, in_granges = T)
saveRDS(meth.mtx, file.path(res.tmp.path, "uncovered_removed_meth_mtx.rds"))

# Differential methylation regions (DMRs) calling -------------------------------------------


## turn meth object to bsseq and find DMR
## https://compepigen.github.io/methrix_docs/articles/methrix.html

## bsseq find DMRs -------------------------------------------
cat(as.character(Sys.time()), ": find DMRs using bsseq\n")
bs_seq <- methrix2bsseq(m = meth)
bs_seq <- BSmooth(BSseq = bs_seq)
saveRDS(bs_seq, file.path(res.tmp.path, "bs_seq.rds"))
bs_seq.tstat <- BSmooth.tstat(BSseq = bs_seq, 
                              group1 = 1:4, 
                              group2 = 5:8, 
                              estimate.var = "group2",
                              mc.cores = 10)
bsseq.dmr <- dmrFinder(bs_seq.tstat)
saveRDS(bsseq.dmr, file.path(res.tmp.path, "bsseq.dmr.rds"))

## dmrseq Find DMRs -------------------------------------------
cat(as.character(Sys.time()), ": find DMRs using dmrseq\n")
loci.idx <- which(DelayedMatrixStats::rowSums2(getCoverage(bs_seq, type="Cov")==0) == 0)
bs_seq <- bs_seq[loci.idx, ]
bs_seq$group <- factor(x = bs_seq$group, levels = c("WTcontrol", "Daxx")) # group2 vs group1
options(MulticoreParam=MulticoreParam(workers = 5))
dmrseq.dmr <- dmrseq(bs = bs_seq,
                     cutoff = 0.05,
                     testCovariate = "group",
                     verbose = T)
saveRDS(dmrseq.dmr, file.path(res.tmp.path, "dmrseq.dmr.rds"))


## metilene Find DMRs -------------------------------------------
### write merged bedgraph for metilene computation
cat(as.character(Sys.time()), ": find DMRs using metilene\n")
methrix::write_bedgraphs(m = meth, 
                         output_dir = pkg.tmp.path, 
                         rm_NA = FALSE, 
                         metilene = TRUE,
                         multiBed = "metline_ip", 
                         phenoCol = "group", 
                         compress = FALSE)

### write command for execute metilene
if (!dir.exists(file.path(pkg.path, "metilene_v0.2-8"))){
  cat("metilene not installed, try to install it locally\n")
  if (file.exists(file.path(pkg.path, "metilene_v02-8.tar.gz"))){
    system(paste0("tar -zxvf ",
                  file.path(pkg.path, "metilene_v02-8.tar.gz "),
                  "--directory=", pkg.path),
           ignore.stdout = T, ignore.stderr = F)
    system(paste0("make --directory ", file.path(pkg.path, "metilene_v0.2-8")),
           ignore.stdout = T, ignore.stderr = T)
    # system(paste0("chmod +x ", metilene.info$metilen.filename))
  }else
    cat("warning: no metilene_** or metilene_**.tar.gz found in", pkg.path, "\n")
}

metilene.info <- list(
  metilen.filename = file.path(pkg.path, "metilene_v0.2-8", "metilene"),
  thread = 4,
  group1 = "Daxx",
  group2 = "WTcontrol",
  bedgraph.filename = file.path(pkg.tmp.path, "metline_ip.bedGraph")
)
metilene.info$output.filename = file.path(res.tmp.path, paste0("metilene.", metilene.info$group1, "_vs_", metilene.info$group2, ".tsv"))
sys.cmd <- paste(metilene.info$metilen.filename, 
                 "-t", metilene.info$thread,
                 "-a", metilene.info$group1,
                 "-b", metilene.info$group2,
                 metilene.info$bedgraph.filename,
                 ">", metilene.info$output.filename)
system(sys.cmd,
       ignore.stdout = F, ignore.stderr = T)

# Annotation --------------------------------------------------------------

cat(as.character(Sys.time()), ": annotate DMRs and matrix\n")
genome = "mm10"
library(annotatr)
library(GenomicRanges)
library(ChIPseeker)
library(TxDb.Mmusculus.UCSC.mm10.knownGene)
library(TxDb.Hsapiens.UCSC.hg19.knownGene)
library(org.Mm.eg.db)
library(org.Hs.eg.db)

## prepare database for annotatr, and save as rds under comAnn.path

# lapply(c("hg19", "mm10"), function(genome){
#   builtin_annotations()[grep(genome, builtin_annotations())]
#   annots = paste0(genome, c('_cpgs', '_basicgenes', '_genes_intergenic',
#                             '_genes_intronexonboundaries'))
#   annots %in% builtin_annotations()[grep(genome, builtin_annotations())]
#   annotations = build_annotations(genome = genome,
#                                   annotations = annots)
#   saveRDS(annotations, file.path(comAnn.path, paste0(genome, "_", "annotations.rds")))
# })


if (genome == "mm10") {
  txdb = TxDb.Mmusculus.UCSC.mm10.knownGene
  annodb <- "org.Mm.eg.db"
  annotations <- readRDS(file.path(comAnn.path, "mm10_annotations.rds"))
}else if(genome == "hg19"){
  txdb = TxDb.Mmusculus.UCSC.mm10.knownGene
  annodb <- "org.Hs.eg.db"
  annotations <- readRDS(file.path(comAnn.path, "hg19_annotations.rds"))
}

## Calculate Overlap between given range and annotation range
CalculateOverlap <- function(df){
  df$Overlap.Length <- pbapply::pbsapply(1:nrow(df), function(i){
    interval1 <- c(df$start[i], df$end[i])
    interval2 <- c(df$annot.start[i], df$annot.end[i])
    
    if (interval1[1] > interval2[2] || interval2[1] > interval1[2]){
      return(0)
    }else if (interval2[1] == interval1[1]){
      return(min(interval1[2]-interval1[1], interval2[2]- interval2[1])+1)
    }else if((interval2[2] - interval1[2])/(interval2[1] - interval1[1]) < 0){
      return(min(diff(interval1), diff(interval2))+1)
    }else{
      return(min(interval1[2]-interval2[1], interval2[2]- interval1[1])+1)
    }
  })
  return(df)
}

## a function to annotate the Grange Object using ChIPseeker and annotatr
ChIPseeker.annotate <- function(dmr.Grange, prefix = ""){
  annotations.ChIPseeker <- annotatePeak(peak = dmr.Grange,
                                         TxDb = txdb, 
                                         annoDb = annodb,
                                         tssRegion = c(-3000, 3000), 
                                         verbose = FALSE)
  annotations.ChIPseeker <- as.data.frame(annotations.ChIPseeker, row.names = NULL)
  annotations.ChIPseeker$annotationtype <- gsub(
    pattern = "(Exon|Intron)(.+)", 
    replacement = "\\1", 
    x = annotations.ChIPseeker$annotation
  )
  annotations.ChIPseeker <- dplyr::arrange(annotations.ChIPseeker, seqnames, start, end)
  annotations.ChIPseeker <- list("Merged" = annotations.ChIPseeker)
  annotations.ChIPseeker <- c(annotations.ChIPseeker, 
                              split(x = annotations.ChIPseeker$Merged, f = annotations.ChIPseeker$Merged$annotationtype))
  openxlsx::write.xlsx(x = annotations.ChIPseeker, 
                       file = file.path(res.tmp.path, paste0("ChIPseeker.", prefix, ".dmr.xlsx")),
                       row.names = F, col.names = T, overwrite = T)
  cat("write annotated DMR table to ", paste0("ChIPseeker.", prefix, ".dmr.xlsx\n"))
  return(annotations.ChIPseeker)
}
annotatr.annotate <- function(dmr.Grange, prefix = ""){
  annotations.annotatr <- annotate_regions(regions = dmr.Grange,
                                           annotations = annotations,
                                           ignore.strand = TRUE,
                                           quiet = FALSE)
  annotations.annotatr <- as.data.frame(annotations.annotatr, row.names = NULL)
  annotations.annotatr <- CalculateOverlap(annotations.annotatr)
  
  annotations.annotatr <- dplyr::arrange(annotations.annotatr, seqnames, start, end)
  annotations.annotatr <- list("Merged" = annotations.annotatr)
  annotations.annotatr <- c(annotations.annotatr,
                            split(x = annotations.annotatr$Merged, f = annotations.annotatr$Merged$annot.type))
  openxlsx::write.xlsx(x = annotations.annotatr, 
                       file = file.path(res.tmp.path, paste0("annotatr.", prefix, ".dmr.xlsx")),
                       row.names = F, col.names = T, overwrite = T)
  cat("write annotated DMR table to ", paste0("annotatr.", prefix, ".dmr.xlsx\n"))
  return(annotations.annotatr)
}

## Annotation for DMRs called by bsseq
bsseq.dmr <- readRDS(file.path(res.tmp.path, "bsseq.dmr.rds"))
bsseq.dmr.Grange <- makeGRangesFromDataFrame(df = bsseq.dmr, keep.extra.columns = T)
ChIPseeker.bsseq.dmr.Grange <- ChIPseeker.annotate(dmr.Grange = bsseq.dmr.Grange, prefix = "bsseq")
annotatr.bsseq.dmr.Grange <- annotatr.annotate(dmr.Grange = bsseq.dmr.Grange, prefix = "bsseq")

## Annotation for DMRs called by dmrseq
dmrseq.dmr.Grange <- readRDS(file.path(res.tmp.path, "dmrseq.dmr.rds"))
ChIPseeker.dmrseq.dmr.Grange <- ChIPseeker.annotate(dmr.Grange = dmrseq.dmr.Grange, prefix = "dmrseq")
annotatr.dmrseq.dmr.Grange <- annotatr.annotate(dmr.Grange = dmrseq.dmr.Grange, prefix = "dmrseq")

## Annotation for DMRs called by metilene
metilene.dmr <- read.table(file = metilene.info$output.filename, sep = "\t", header = F)
colnames(metilene.dmr) <- c("chr", "start", "end", "qvalue", "mean_diff", "CpGs", "p_MWUtest", "p_KStest", "mean_g1", "mean_g2")
metilene.dmr.Grange <- makeGRangesFromDataFrame(df = metilene.dmr, keep.extra.columns = T)
ChIPseeker.metilene.dmr.Grange <- ChIPseeker.annotate(dmr.Grange = metilene.dmr.Grange, prefix = "metilene")
annotatr.metilene.dmr.Grange <- annotatr.annotate(dmr.Grange = metilene.dmr.Grange, prefix = "metilene")


## Annotation for Methylation Matrix

meth.mtx <- annotatePeak(peak = meth.mtx,
                         TxDb=txdb, 
                         annoDb = annodb,
                         tssRegion=c(-3000, 3000), 
                         verbose=FALSE)
meth.mtx <- as.data.frame(meth.mtx)

write.table(x = meth.mtx, 
            file = file.path(res.tmp.path, "annotated.uncovered_removed_meth.mtx.txt"),
            quote = F, sep = "\t", row.names = F, col.names = T)

# Figures -----------------------------------------------------------------

cat(as.character(Sys.time()), ": start to plot\n")
## coverage distribution figure
uncovered_removed_meth <- readRDS(file.path(res.tmp.path, "uncovered_removed_meth.rds"))
meth_stats <- get_stats(m = uncovered_removed_meth)
write.table(as.data.frame(lapply(meth_stats, as.character)),
            file.path(res.tmp.path, "meth_stat.txt"),
            sep = "\t", row.names = F, col.names = T, quote = F)
meth_stats.chrM.removed <- subset(meth_stats, Chromosome != "chrM")
pdf(file.path(res.tmp.path, "CoveragePlot.pdf"), width = 12, height = 8)
lapply(unique(meth_stats$Sample_Name), function(sample_name){
  plot_stats(plot_dat = meth_stats.chrM.removed, 
             what = "C", 
             stat = "mean", 
             samples = sample_name)
})
invisible(dev.off())

## heatmap on high variable methylation site
library(ComplexHeatmap)

meth.mtx <- methrix::get_matrix(m = uncovered_removed_meth, 
                                type = "M")
var <- apply(meth.mtx, 1, sd)
HVmeth.mtx <- meth.mtx[order(var, decreasing = T)[1: 1500], ]
colnames(HVmeth.mtx) <- substr(colnames(HVmeth.mtx), 1, 9)
anncol <- data.frame(group = factor(x = sample_anno$group[match(colnames(HVmeth.mtx), sample_anno$sampleID)],
                                    levels = unique(sample_anno$group)),
                     row.names = colnames(HVmeth.mtx))
anncolor <- list()
anncolor[["group"]] <- c("Daxx" = "#B2DF8A", "WTcontrol" = "#33A02C")
pdf(file.path(res.tmp.path, "heatmap.pdf"))
pheatmap(HVmeth.mtx, 
         color = c("#6699CC","white","#FF3C38"),
         annotation_col = anncol, 
         annotation_colors = anncolor,
         name = "methylation")
invisible(dev.off())

