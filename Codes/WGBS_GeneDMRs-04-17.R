# if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
# 
# BiocManager::install(c("devtools", "clusterProfiler", "corrplot", "dplyr", "ffbase", "genomation", 
#                        "pheatmap", "plotrix", "qqman", "RCircos", "VennDiagram", "org.Mm.eg.db"))
# devtools::install_github("xiaowangCN/GeneDMRs")



Sys.setenv("R_REMOTES_NO_ERRORS_FROM_WARNINGS" = "true")
workdir <- "/home/a208/data/cwx/Xlu/Amanda_WGSmb"; setwd(workdir)
# workdir <- "i:/genomicdata/Xlu/Amanda_WGSmb";setwd(workdir)
fig.path <- file.path(workdir, "Figures")
res.path <- file.path(workdir, "Results")
data.path <- file.path(workdir, "InputData")
comAnn.path <- file.path(workdir, "Annotation") ## path containning annotation files created by Jack
WGSmb.path <- file.path(data.path, "Amanda_WGSmb_Project_Feb2022")
res.tmp.path <- file.path(res.path, "GeneDMRs") ## result provided by GeneDMRs

if (!file.exists(res.path)) { dir.create(res.path) }
if (!file.exists(fig.path)) { dir.create(fig.path) }
if (!file.exists(data.path)) { dir.create(data.path) }
if (!file.exists(comAnn.path)) { dir.create(comAnn.path) }
if (!file.exists(res.tmp.path)) { dir.create(res.tmp.path) }

# the file path of compressed bedGraph/bismark.cov file is put under WGSmb.path
# WGSmb.path/**.bedGraph.gz
# WGSmb.path/**.bismark.cov.gz
# put all sample under one folder


library(GeneDMRs)
library(ggplot2)
library(data.table)
library(dplyr)
library(ff)


# System set --------------------------------------------------------------

cat(as.character(Sys.time()), ": read files\n")
## download methdata from github
## put the folder named "methdata" under the path: system.file(package = "GeneDMRs")
system.file(package = "GeneDMRs")

## read the provided refseq, cyto, cpgifeature file #
inputrefseqfile <- Bedfile_read(paths = paste(system.file(package = "GeneDMRs"), "/methdata", sep=""), 
                                bedfile = "refseq", suffix = ".txt", feature = FALSE)
inputcytofile <- Cytofile_read()
inputcpgifeaturefile <- Bedfile_read(paths = paste(system.file(package = "GeneDMRs"), "/methdata", sep = ""),
                                     bedfile = "cpgi", feature = TRUE, featurewrite = FALSE)

## check the data is mouse (default) or human
table(inputcpgifeaturefile$chr)
table(inputcytofile$chr)


# Read file ---------------------------------------------------------------

## set 
sample.path <- list.files(path = WGSmb.path,
                          pattern = "bismark.cov.gz", full.names = T)
sample_anno <- read.table(file.path(data.path,"WGSmb_sampleInfo_Feb2022.txt"),
                          sep = "\t",row.names = NULL,check.names = F,stringsAsFactors = F,header = T)
sample_anno$path <- sample.path[match(sample_anno$sampleID, substr(basename(sample.path), 1, 9))]
controls = sample_anno$path[sample_anno$group == "WTcontrol"]
cases = sample_anno$path[sample_anno$group == "Daxx"]
inputmethfile <- Methfile_read(control_paths = controls,
                               case_paths = cases,
                               suffix = "cov.gz", 
                               WGBS = T)
## quality control
inputmethfile_QC <- Methfile_QC(inputmethfile)
table(inputmethfile_QC$chr)

# Find Differential Methylation Regions --------------------------------------------

cat(as.character(Sys.time()), ": find DMRs\n")
# allDMGs <- Quick_GeneDMRs(paths = paste(system.file(package = "GeneDMRs"), "/methdata", sep=""))
regiongeneall <- Methmean_region(inputmethfile_QC = inputmethfile_QC, 
                                 inputrefseqfile = inputrefseqfile, 
                                 chrnum = "all")
regiongeneall_Qvalue <- Logic_regression(regiongeneall)
regiongeneall_significant <- Significant_filter(regiongeneall_Qvalue)
saveRDS(regiongeneall_Qvalue, file.path(res.tmp.path, "regiongeneall_Qvalue.rds"))
openxlsx::write.xlsx(list("region_sig" = regiongeneall_significant, 
                          "region_all" = regiongeneall_Qvalue),
                     file.path(res.tmp.path, "GeneDMR.xlsx"), overwrite = T)


# Figure ------------------------------------------------------------------

cat(as.character(Sys.time()), ": output Figure\n")
library(clusterProfiler)

## summary 
pdf(file.path(res.tmp.path, "DMR_Summary.pdf"))
Group_boxplot(regiongeneall)
Chromosome_pieplot(regiongeneall_significant, title = "Distribution of DMG")
Heatmap_plot(regiongeneall_significant, title = "Methylation level (%) for genes", 
             display_numbers = FALSE)
invisible(dev.off())

## Enrichplot
Dbannotation <- "org.Mm.eg.db"
keggorganism <- "mmu"
listnum <- 20

### Enrich plot: GOgroup
eg_all <- bitr(regiongeneall$id, fromType = "REFSEQ", toType = c("ENTREZID", "SYMBOL"), OrgDb = Dbannotation)
eg_all <- dplyr::distinct(eg_all, ENTREZID, .keep_all = TRUE)
genefile <- as.vector(unlist(regiongeneall_significant[1]))
eg <- bitr(genefile, fromType = "REFSEQ", toType = c("ENTREZID", "SYMBOL"), OrgDb = Dbannotation)
eg <- dplyr::distinct(eg, ENTREZID, .keep_all = TRUE)
all(eg$REFSEQ %in% regiongeneall_significant$id)

#### Enrich plot: GOgroup for hyper-methylated gene
eg_hyper <- eg[eg$REFSEQ %in% regiongeneall_significant$id[regiongeneall_significant$Methdiff1 > 0], ]
ego_hyper <- enrichGO(gene = eg_hyper$ENTREZID, universe = eg_all$ENTREZID,
                      OrgDb = Dbannotation, ont =  "BP", 
                      readable = T)
hyper_BP <- barplot(ego_hyper, showCategory=20, title="Biological process for hyper-methylated gene")

#### Enrich plot: GOgroup for hypo-methylated
eg_hypo <- eg[eg$REFSEQ %in% regiongeneall_significant$id[regiongeneall_significant$Methdiff1 < 0], ]
ego_hypo <- enrichGO(gene = eg_hypo$ENTREZID, universe = eg_all$ENTREZID,
                     OrgDb = Dbannotation, ont =  "BP", 
                     readable = T)
hypo_BP <- barplot(ego_hypo, showCategory=20, title="Biological process for hypo-methylated gene")

#### Enrich plot: GOgroup output pdf
pdf(file.path(res.tmp.path, "BPforSignificantGene.pdf"), width = 9, height = 9)
hyper_BP + scale_y_discrete(labels=function(x) yulab.utils::str_wrap(x, width=50))
hypo_BP + scale_y_discrete(labels=function(x) yulab.utils::str_wrap(x, width=50))
invisible(dev.off())


## Enrich plot: GO
DMgene_merge <- data.frame(REFSEQ = regiongeneall_significant$id, 
                           Methdiff = regiongeneall_significant$Methdiff1)
eg_merge <- merge(x = eg, y = DMgene_merge, by = "REFSEQ")
neweg <- data.frame(Entrez = eg_merge$ENTREZID, 
                    Methdiff = eg_merge$Methdiff)
neweg$group <- "Hyper-methylated"
neweg$group[neweg$Methdiff < 0] <- "Hypo-methylated"
formula_go <- compareCluster(Entrez ~ group, 
                             data = neweg, fun = "enrichGO", ont = "all", 
                             pvalueCutoff = 0.1, OrgDb = Dbannotation)
saveRDS(formula_go, file.path(res.tmp.path, "formula_go.rds"))
dotplot(formula_go, x = "group", 
        color = "p.adjust", showCategory = listnum, 
        split = NULL, font.size = 10, title = "Go term for significant gene")+ 
  scale_y_discrete(labels=function(x) yulab.utils::str_wrap(x, width=70))
ggsave(file.path(res.tmp.path, "GoTermforSignificantGene.pdf"), width = 9, height = 10)


# Circular Graph ----------------------------------------------------------

cat(as.character(Sys.time()), ": plot circular graph\n")
library(RCircos)

## circular graph: prepare plotdata (a list) from data
Circos_plot.prep <- function(inputcytofile = inputcytofile,
                             inputmethfile_QC = inputmethfile_QC,
                             inputrefseqfile = inputrefseqfile,
                             labelname = labelname,
                             inputcpgifeaturefile = inputcpgifeaturefile,
                             linecolor = NULL){
  windowfile <- Window_divide(inputcytofile, windowbp = 1e+06)
  windowfilealls <- Methmean_region(inputmethfile_QC, windowfile, 
                                    chrnum = "all")
  chrpos <- data.frame(chr = windowfilealls$chr, start = windowfilealls$start, 
                       end = windowfilealls$end)
  groupnum <- length(grep("Methgroup", colnames(windowfilealls)))
  grouppos <- grep("Methgroup", colnames(windowfilealls))
  inputcpgi <- filter(inputcpgifeaturefile, cpgfeature %in% 
                        "CpGisland")
  inputshore <- filter(inputcpgifeaturefile, cpgfeature %in% 
                         "Shores")
  out_density <- data.frame(chrpos, density = array(0, c(nrow(chrpos), 
                                                         3)))
  for (i in 1:nrow(chrpos)) {
    filter1 <- filter(inputrefseqfile, chr %in% as.vector(unlist(out_density[1]))[i])
    filter2 <- filter(inputcpgi, chr %in% as.vector(unlist(out_density[1]))[i])
    filter3 <- filter(inputshore, chr %in% as.vector(unlist(out_density[1]))[i])
    out_density[i, (ncol(chrpos) + 1)] <- sum(filter1$start > 
                                                out_density[i, 2] & filter1$end <= out_density[i, 
                                                                                               3])
    out_density[i, (ncol(chrpos) + 2)] <- sum(filter2$start > 
                                                out_density[i, 2] & filter2$end <= out_density[i, 
                                                                                               3])
    out_density[i, (ncol(chrpos) + 3)] <- sum(filter3$start > 
                                                out_density[i, 2] & filter3$end <= out_density[i, 
                                                                                               3])
  }
  labelnamefile <- data.frame(chr = labelname$chr, start = labelname$start, 
                              end = labelname$end, label = labelname$id)
  
  RC.RawData <- list("inputcytofile" = inputcytofile,
                     "out_density" = out_density,
                     "labelnamefile" = labelnamefile,
                     "windowfilealls" = windowfilealls,
                     "grouppos" = grouppos,
                     "groupnum" = groupnum,
                     "linecolor" = linecolor)
  return(RC.RawData)
}

labelname <- dplyr::arrange(regiongeneall_significant, Qvalue1)
labelname <- split(labelname, labelname$chr)
labelname <- lapply(labelname, function(x){x[1:10,]})
labelname <- do.call(rbind, labelname)
cir.plot.data <- Circos_plot.prep(inputcytofile = inputcytofile,
                                  inputmethfile_QC = inputmethfile_QC,
                                  inputrefseqfile = inputrefseqfile,
                                  labelname = labelname,
                                  inputcpgifeaturefile = inputcpgifeaturefile,
                                  linecolor = c("blue1", "green1"))

## circular graph: output pdf
pdf(file.path(res.tmp.path, "Circular_Graph.pdf"), width = 15, height = 15)
RCircos.Set.Core.Components(cir.plot.data$inputcytofile, chr.exclude = NULL, 
                            tracks.inside = 7 + cir.plot.data$groupnum, tracks.outside = 0)
RC.param <- RCircos.Get.Plot.Parameters()
RC.param["text.size"] <- 0.4
RC.param["scatter.color"] <- "purple"
RC.param["track.background"] <- "white"
RC.param["grid.line.color"] <- "white"
RCircos.Reset.Plot.Parameters(RC.param)
RCircos.Set.Plot.Area()
par(mai = c(0, 0.1, 0, 0.1))
plot.new()
plot.window(c(-2.5, 2.5), c(-2.5, 2.5))
RCircos.Chromosome.Ideogram.Plot()
RCircos.Gene.Connector.Plot(cir.plot.data$labelnamefile, track.num = 1, side = "in", outside.pos = 10)
RCircos.Gene.Name.Plot(cir.plot.data$labelnamefile, name.col = 4, track.num = 2, side = "in")
RCircos.Histogram.Plot(cir.plot.data$out_density, data.col = 4, track.num = 5, side = "in")
RCircos.Scatter.Plot(cir.plot.data$out_density, data.col = 5, track.num = 6, side = "in")
RC.param <- RCircos.Get.Plot.Parameters()
RC.param["scatter.color"] <- "purple1"
RCircos.Reset.Plot.Parameters(RC.param)
RCircos.Scatter.Plot(cir.plot.data$out_density, data.col = 6, track.num = 7, side = "in")
for (j in 1:cir.plot.data$groupnum) {
  RC.param <- RCircos.Get.Plot.Parameters()
  RC.param["line.color"] <- cir.plot.data$linecolor[j]
  RCircos.Reset.Plot.Parameters(RC.param)
  RCircos.Line.Plot(cir.plot.data$windowfilealls[, -1], data.col = cir.plot.data$grouppos[j] - 1, track.num = 7 + j, side = "in")
}
invisible(dev.off())


# Annotation --------------------------------------------------------------

cat(as.character(Sys.time()), ": annotate DMRs\n")
## Prepare annotation database
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

## Read DMR table
Sheetnames <- openxlsx::getSheetNames(file.path(res.tmp.path, "GeneDMR.xlsx"))
tab <- lapply(Sheetnames, function(sheetname){
  openxlsx::read.xlsx(xlsxFile = file.path(res.tmp.path, "GeneDMR.xlsx"),
                      sheet = sheetname)
})
names(tab) <- Sheetnames

## Annotation for significantly differential methylation regions ----------------------------

gene.Grange <- makeGRangesFromDataFrame(tab$region_sig, keep.extra.columns = T)

### Annotation for significantly differential methylation regions using ChIPseeker
annotations.ChIPseeker <- annotatePeak(peak = gene.Grange,
                                       TxDb=txdb, 
                                       annoDb = annodb,
                                       tssRegion=c(-3000, 3000), 
                                       verbose=FALSE)
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
                     file = file.path(res.tmp.path, "ChIPseeker.GeneDMRs.sig.dmr.xlsx"),
                     row.names = F, col.names = T, overwrite = T)


### Annotation for significantly differential methylation regions using annotatr
annotations.annotatr <- annotate_regions(regions = gene.Grange,
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
                     file = file.path(res.tmp.path, "annotatr.GeneDMRs.sig.dmr.xlsx"),
                     row.names = F, col.names = T, overwrite = T)


## Annotation for all differential methylation regions ----------------------

gene.Grange <- makeGRangesFromDataFrame(tab$region_all, keep.extra.columns = T)

### Annotation for all differential methylation regions using ChIPseeker
annotations.ChIPseeker <- annotatePeak(peak = gene.Grange,
                                       TxDb=txdb, 
                                       annoDb = annodb,
                                       tssRegion=c(-3000, 3000), 
                                       verbose=FALSE)
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
                     file = file.path(res.tmp.path, "ChIPseeker.GeneDMRs.all.dmr.xlsx"),
                     row.names = F, col.names = T, overwrite = T)


### Annotation for all differential methylation regions using annotatr
annotations.annotatr <- annotate_regions(regions = gene.Grange,
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
                     file = file.path(res.tmp.path, "annotatr.GeneDMRs.all.dmr.xlsx"),
                     row.names = F, col.names = T, overwrite = T)
