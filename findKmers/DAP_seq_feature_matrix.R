
library(tidyverse)
library(ChIPseeker)
library(clusterProfiler)
library(GenomicFeatures)

rm(list = ls())
setwd('E:/Study/Courses/Projects/ResearchProjects/TimeSeriesDrought_PEG/ML_pipeline')
# Load data
samplefiles <- list.files("dap_seq_all_peaks/", pattern= ".narrowPeak", full.names=T)
samplefiles <- as.list(samplefiles)
names(samplefiles) <- sapply(samplefiles,function(x) substr(strsplit(x, split='/')[[1]][2],1,9))

# convert the reference genome data
txdb <- makeTxDbFromGFF('TAIR10_GFF3_genes.gff')
# annotate peaks
peakAnnoList <- lapply(samplefiles, annotatePeak, TxDb=txdb, 
                       tssRegion=c(-1000, 0), verbose=FALSE)

################## Merge peaks to make peak matrix######################################
peakMatrix <- data.frame(geneId=0)
for (name in names(peakAnnoList)){
  
  peak_anno <- data.frame(peakAnnoList[[name]]@anno)
  peak_anno <- peak_anno[peak_anno$annotation=='Promoter','geneId']
  peak_df <- data.frame(geneId=peak_anno, name_tmp=rep(1,length(peak_anno)))
  colnames(peak_df) <- c('geneId',name)
  peakMatrix <- merge(peakMatrix, peak_df, by='geneId',all=TRUE)
  peakMatrix <- distinct(peakMatrix)
}
peakMatrix <- peakMatrix[-1,]
peakMatrix[is.na(peakMatrix)] <- 0
write.csv(peakMatrix,'DAP_seq_PeakMatrix.csv')


peakMatrix <- read.csv('DAP_seq_PeakMatrix.csv',row.names = 1)
control_genes <- read.csv('df_controlGenes.csv',row.names=1)
DEG_genes <- read.csv('../kmers/dehydration_allDE.csv',row.names = 1)
c <- rownames(control_genes)
d <- rownames(DEG_genes)
dif <- setdiff(union(c,d),rownames(peakMatrix))
dif.dataframe <- data.frame(matrix(0, nrow = length(dif),ncol = length(colnames(peakMatrix))))
colnames(dif.dataframe) <- colnames(peakMatrix)
rownames(dif.dataframe) <- dif
final_peak_matrix <- rbind(peakMatrix[intersect(rownames(peakMatrix),union(c,d)),],dif.dataframe)
write.csv(final_peak_matrix,'allPeaks_matrix.csv')

# Atalias <- read.table('gene_aliases_20181231.txt',sep='\t',quote='',fill = T,row.names=NULL) 
# 
# dd <- colnames(peakMatrix)[colnames(peakMatrix) %in% Atalias[,1]]
# Atalias[which(Atalias$V1 %in% dd),2]

TFnames <- read.table('TFnames.txt',sep=',',fill = T,header = T)
