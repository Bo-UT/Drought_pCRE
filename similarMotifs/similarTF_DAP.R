library(tidyverse)
library(reshape2)
library(ComplexHeatmap)
library(MotifDb)
library(universalmotif)
library(ggseqlogo)

rm(list=ls())

setwd('E:/Study/Courses/Projects/ResearchProjects/labDehydrationProject/ML_pipeline/bamcompare_importantTF/scripts')

# load DAP motifs
dapMotifFiles <- list.files('../data/DAPmotifs',full.names = T,pattern = '.txt',recursive = T)
dapMotifs <- sapply(dapMotifFiles,read_meme)
names(dapMotifs) <- unlist(lapply(names(dapMotifs), function (x) strsplit(x,'/')[[1]][5]))
length(dapMotifs)
# extract top 10 motifs for each gene cluster
top10CREfiles <- grep(list.files('../../findSimilarMotif/results/similarTF/similarTF'),
                      pattern='heatmap.tiff',invert=T,value = T) # exclude heatmap files

extractTop10CRE <- function (cluster) {
  
  clusterfile = grep(top10CREfiles, pattern = paste0('cluster_', cluster, '_'), value = T)
  testMatrix = matrix(0, nrow=10, ncol=2)
  for (k in seq(10)) testMatrix[k,] = unlist(strsplit(strsplit(clusterfile[k], '[.]')[[1]][1],'_'))[c(3,4)]
  # write top10 motifs to the file
  writeLines(testMatrix[order(testMatrix[,2]),][,1], paste0('../data/top10motifs/cluster_',cluster,'_top10motifs.txt' ))
}

for (i in seq(40)) extractTop10CRE(i)


##compare with DAP motifs to find similar TF with PCC larger than 0.8

compare_motif <- function(cre){
  cre <- create_motif(cre, pseudocount = 1, nsites = 50)
  res <- compare_motifs(c(dapMotifs,cre),method = 'PCC',min.mean.ic = 0,tryRC=TRUE,
                        score.strat = "a.mean")
  TF_highPCC <- res[dim(res)[1],] %>% # take the last row that compares cre with all other TF binding sites
    as.data.frame() %>% 
    rownames_to_column('TF') %>% 
    # mutate(TF=sapply(TF,function(x) unlist(strsplit(x,split=' '))[1])) %>% # split TF names with [duplicated..]
    stats::setNames(c('TF','pcc')) %>% 
    mutate(pcc=round(pcc,4)) %>% 
    distinct() %>% # remove duplicated
    dplyr::filter(!grepl('motif',TF)) %>%  # filter TF with pcc score larger than 0.8
    dplyr::filter(pcc > 0.8)
    # arrange(desc(pcc)) %>% # sort in descending order
    # dplyr::slice(2:12) #extract top 4, the first one is the CRE itself
  return(TF_highPCC) 
}

for (afile in list.files('../data/top10motifs',full.names = T)) {
  cluster = strsplit(basename(afile),'_')[[1]][2]
  data = read.table(afile,header = F)
  results = apply(data, 1, function(x) compare_motif(x))
  df = do.call(rbind, results)
  df = df %>% distinct(TF, .keep_all=T)
  write.table(df, file = paste0("../data/highpccTF/cluster_",cluster,"_highpccTF.txt"), sep = "\t",
              row.names = F, col.names = F,quote = F)
  
}


####match motifs with SraRun 
motifName1 = names(dapMotifs)
motifName2 = unlist(lapply(names(dapMotifs),function(x) dapMotifs[[x]]['name']))           
df = data.frame(motifName1=motifName1,motifName2=motifName2)
dim(df2)
# df[duplicated(df$motifName1),]
dapSra = read.csv('DAP_seq_SraRun2.csv',row.names = 1)
df2 = data.frame(gene_id = dapSra$gene_id, protein_family = dapSra$protein_family,
                 motifName1 = paste0(dapSra$Protein,'_',dapSra$dna_source,'_',dapSra$Replicate),
                 Run = dapSra$Run)
df2 %>% inner_join(df,by='motifName1') %>% select(-motifName1) %>% dim

# df2[df2$motifName1 %in% setdiff(df2$motifName1,df$motifName1),]
# df[df$motifName1 %in% setdiff(df$motifName1,df2$motifName1),]
# df[!(df$motifName1 %in% df2$motifName1),]

# select TF with biologica repeats
df2 %>% mutate(motifName3=paste0(dapSra$Protein,'_',dapSra$dna_source)) %>% 
  group_by(motifName3) %>% filter(n()>1) %>% ungroup %>% 
  select(-motifName1) %>% rename(motifName3 = 'motifName1') %>% inner_join(df,by='motifName1') %>% 
  # filter(str_detect(motifName1, "At1g22810")) %>% 
  select(-motifName1) %>% 
  bind_rows(df2 %>% inner_join(df,by='motifName1') %>% select(-motifName1)) %>% 
  # count(motifName2) %>% filter(n>1) #detect duplicated rows
  rename(motifName2='motifName') %>% 
  # test to select Run for TF in cluster 10, unname to convert a row/column to a vector without name
  # filter(motifName %in% unname(unlist(read.table('../data/highpccTF/cluster_10_highpccTF.txt') %>% select(V1))))
  write.table('DAP_motifName_Run.txt',sep='\t', row.names = F,col.names = F,quote = F)

Sra_Run = read.table('DAP_motifName_Run.txt')

for (afile in list.files("../data/highpccTF/",full.names = T)){
  cluster = unlist(strsplit(basename(afile),'_'))[2]
  Sra_Run %>% filter(V4 %in% unname(unlist(read.table(afile) %>% select(V1)))) %>%
    mutate(V5=paste0(V3, '_',sapply(.$V4,function(x) unlist(strsplit(x, '_'))[3]))) %>% select(c(V1,V2,V5,V4)) %>%   
    write.table(paste0("../data/highpccTF/cluster_",cluster,"Run.txt"), sep='\t', row.names = F,col.names = F,quote = F)
}


Sra_Run %>% filter(V4 %in% unname(unlist(read.table('../data/highpccTF/cluster_1_highpccTF.txt') %>% select(V1)))) %>%
  mutate(V5=paste0(V3, '_',sapply(.$V4,function(x) unlist(strsplit(x, '_'))[3]))) %>% select(c(V1,V2,V5,V4)) %>% 
  write.table(paste0("../data/highpccTF/cluster_1_Run.txt"), sep='\t', row.names = F,col.names = F,quote = F)


df2[duplicated(df2$gene_id),]
df[duplicated(df$motifName1),]
df2[grepl('SOL1',df2$motifName1),]
df2[df2$gene_id=='AT1G69560',]
df2[grepl('ERF15_col_',df2$motifName1),]
df[grepl('ERF15_col$',df$motifName1),]
df2
dim(df2)
df[which(sapply(strsplit(df$motifName1,'_'),function(x) length(x))=='2'),]
df2[which(sapply(strsplit(df2$motifName1,'_'),function(x) length(x))=='2'),]
df[grepl('ANAC017',df$motifName1),]
dapSra[grepl('ANAC017',dapSra$Protein),]

##extract Run(SRR file) and library(Single or paired end)

srarun = read.csv('DAP_seq_sra.csv',row.names = 1)
View(srarun)
read.csv('../data/highpccTF/cluster_10Run.txt',sep='\t') %>% select(3) %>% 
  apply(.,1,function (x) strsplit(x, split='_')[[1]][1])

colnames(srarun)
srarun[srarun$LibraryLayout=='PAIRED',]
