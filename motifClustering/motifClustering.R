# motif clustering 
# rsat: matrix-clustering

library(MASS)
library(universalmotif)
library(stats)
library(dendextend) 
library(tidyverse)
library(MotifDb)
library(ggseqlogo)
# dendextend to change parameters by
# set(object, what, value)
library(colorspace) # get nice colors
setwd('E:/Study/Courses/Projects/ResearchProjects/labDehydrationProject/ML_pipeline/motifsClustering/')

rm(list=ls())

################################################
############Generate jaspar. format file########
################################################

afile <- file('cluster_2_motif_jaspar.txt')
# input a motif matrix from universalmotif create_motif
motif_matrix <- function(cis){
  motif <- create_motif(cis, name=cis)
  motif <- motif['motif']
  write(paste0('>', cis), 'cluster_2_motif_jaspar.txt', append = T)
  for (i in 1:dim(motif)[1]){
    currow=paste(rownames(motif)[i], '[', sep=' ')
    for (j in 1:dim(motif)[2]){
      currow <- paste(currow,motif[i,j], sep=' ')
      
    }
    currow = paste(currow, ']', sep=' ')
    write(currow,'cluster_2_motif_jaspar.txt',append=T)
  }
  
}

kmer_files <- list.files('data/kmersResults/')[grep('df_cluster_',list.files('data/kmersResults/'))]
test <- read.csv(paste0('data/kmersResults/',kmer_files[2]), row.names = 1)

for (k in 1:length(colnames(test))){
  motif_matrix(colnames(test)[k])
  
}

##############################################
############Get PCC matrix####################
##############################################
################# input: kmers list
###################### output: dendrogram, cluster table
# compute the PCC matrix
PCCmatrix <- function(kmerList){
  kmers <- c()
  for (kmer in kmerList){
    kmers <- c(kmers,create_motif(kmer,name=kmer,pseudocount = 1, nsites = 50))
  }
  return(compare_motifs(kmers, method='PCC',min.mean.ic = 0,tryRC=TRUE,
                        score.strat = "a.mean"))
}

# dendrogram plot and return cluster table
main <- function(k){
  maindir <- paste0('cluster_', k)
  # create a folder
  # dir.create(file.path(mainDir, subDir))
  dir.create(maindir)
  # load data
  curCluster <- read.csv(paste0('data/kmersResults/',kmer_files[k]), row.names = 1)
  pcc <- PCCmatrix(colnames(curCluster))
  dend <- as.dist(1-pcc) %>% hclust(method = 'average') %>% as.dendrogram  
  length(labels(dend))
  # cutree with threshold 
  cluster <- cutree(dend, h=threshold)
  cluster_num <- length(unique(cluster))
  # dend <- color_branches(dend, k=5) # Color the branches based on the clusters
  # color labels based on clusters
  labels_colors(dend) <-
    rainbow_hcl(cluster_num)[sort_levels_values(
      cluster[order.dendrogram(dend)]
    )]
  tiff(paste0(maindir,'/cluster_',k,'.tiff'),res=300,width = 2000,height = 3000)
  maxcex <- 0.65
  cex <- maxcex-length(labels(dend))/1000
  par(mar=c(5,6,4,8),cex=cex)
  plot(dend, main='',
       horiz =  TRUE,  nodePar = list(cex = .007))
  abline(v=0.39, col=2,lty=2)
  par(cex=1.2)
  title(main = paste0('pCREs of cluster ', k))
  dev.off()
  
  #########Adjust dendrogram label and title size
  # par(cex=0.3, mar=c(5, 8, 4, 1))
  # plot(hc, xlab="", ylab="", main="", sub="", axes=FALSE)
  # par(cex=1)
  # title(xlab="xlab", ylab="ylab", main="main")
  # axis(2)
  #########
  
  # write out the clusters
  write.table(as.data.frame(cluster),paste0(maindir,'/cluster_',k,'.txt'),sep='\t')
  
}

kmer_files <- list.files('data/kmersResults/')[grep('df_cluster_',list.files('data/kmersResults/'))]
threshold <- 0.39
for (i in 1:length(kmer_files)) main(i)

######################pCRE merge and find similar TF
motifFiles <- list.files('results/', pattern='.txt', recursive = T,full.names = T)
# cluster 21
clu21 <- read.table(motifFiles[grep('cluster_21',motifFiles)],sep = '\t') %>% 
  rownames_to_column('motif')

# input a dataframe with two columns: motif, cluster
mergemotifs <- function(df){
  clus_num <- length(unique(df$cluster))
  motif_names <- paste0('cluster ',seq(clus_num))
  motifs <- c()
  for (k in seq(clus_num)){
    cur_df <- df[df$cluster==k,]
    cur_cluster_motif_num <- dim(cur_df)[1]
    if (cur_cluster_motif_num > 1){
      motif_to_merge <- c()
      for (mo in cur_df$motif) motif_to_merge <- c(motif_to_merge,create_motif(mo))
      motif_merged <- merge_motifs(motif_to_merge,method = "PCC")
      motif_merged['name'] <- motif_names[k]
      motifs <- c(motifs,motif_merged)
    } else {
      # in case there is only one motif
      singleMotif <- create_motif(cur_df$motif)
      singleMotif['name'] <- motif_names[k]
      motifs <- c(motifs,singleMotif)
    }
  }
  return(motifs)
}

clu21_motifs <- mergemotifs(clu21)
view_motifs(clu21_motifs)
dev.off()


logo_plot <- function(k, merged_motifs){
  ggplot() + geom_logo(merged_motifs[[k]]['motif'])+theme_logo()+ 
    theme(plot.title = element_text(hjust = 0.5))+
    ggtitle(merged_motifs[[k]]['name'])
}

tiff('cluster21_merged_motifs.tiff', res=300,height = 3500, width = 1000)
plot_list <- lapply(seq(length(clu21_motifs)),function(x) logo_plot(x,clu21_motifs))
do.call(gridExtra::grid.arrange, c(plot_list, ncol=2))
dev.off()
