library(tidyverse)
library(reshape2)
library(ComplexHeatmap)
library(MotifDb)
library(universalmotif)
library(ggseqlogo)
library(cowplot)
library(ade4)

rm(list=ls())

setwd('E:/Study/Courses/Projects/ResearchProjects/labDehydrationProject/ML_pipeline/DAP_pCRE_RF/Scripts/')

files <- list.files("../Results/CRE_DAP_RF_results",pattern = '.csv', full.names = T)
CRE_files <- files[grep('df_CRE_RF|df_CRE_Kbest',files)]

plot_rank <- function(afile){
  # read data
  df <- read.csv(afile,row.name=1)
  clu_num <- strsplit(strsplit(basename(afile), split='[.]')[[1]][1], split='_') %>% unlist()
  clu_num <- clu_num[length(clu_num)]
  methodname <- strsplit(basename(afile),'_')[[1]][3]
  
  # print(clu_num)
  # print(methodname)
  # get rank for each pCRE
  getRank <- function(x){rank(x,ties.method='average')}
  # make the heatmap
  heatmap_cluster <- df %>% apply(.,1,getRank) %>% t() %>% 
    Heatmap(name = "Rank", 
            show_row_names = FALSE,
            row_title='Repeat',
            column_title = paste('Cluster', clu_num,methodname,collapse = '_'),
            column_names_rot=90,
            column_names_gp = gpar(fontsize = 4))
  
  # save the plot
  tiff(paste0('../Results/cluster_',clu_num,'_',methodname,'_rankheatmap.tiff'),
       units="in", width=20, height=10, res=400)
  print(heatmap_cluster)
  dev.off()
  
  # extract left most 10 pCRE (top 10)
  top10 <- colnames(df)[column_order(heatmap_cluster)][1:10] # row_order(heatmap_clu1)
  tmpcolname <- paste0(methodname,'_',clu_num)
  df <- data.frame(row.names = seq(10))
  df[`tmpcolname`] <- top10
  return(df)
}

top10motifs <- data.frame(row.names = seq(10))
for (afile in CRE_files) top10motifs <- cbind(top10motifs,plot_rank(afile))

orderednames <- c()
for (i in seq(40)) orderednames <- c(orderednames,paste0(c('Kbest_','RF_'),i))
top10motifs <- top10motifs[,orderednames]
top10motifs_clu21 <- top10motifs[,colnames(top10motifs)[grep('_21', colnames(top10motifs))]]

kbest_motifs <- c()
for (motif in top10motifs_clu21$Kbest_21) kbest_motifs <- c(kbest_motifs,create_motif(motif,name = motif))
RF_motifs <- c()
for (motif in top10motifs_clu21$RF_21) RF_motifs <- c(RF_motifs,create_motif(motif,name = motif))

comp <- compare_motifs(c(kbest_motifs,RF_motifs),method = "PCC", min.mean.ic = 0,
                 score.strat = "a.mean")
comp <- comp[top10motifs_clu21$Kbest_21,top10motifs_clu21$RF_21]

pheatmap::pheatmap(comp,
                      cluster_cols=F,
                      cluster_rows=F,
                   display_numbers=T,
                   main = "Cluster21_top10motifs_kbest_RF", 
                   legend_labels = c("0", "0.2", "0.4", "0.6","0.8","1.0","PCC\n"),
                   cellwidth=25, cellheight=25,
                   number_color='black',
                   filename = '../Results/cluster21_top10motifs_pcc.tiff')
print(p)


##############################
library(motifStack)

############example from motifstack
############https://www.bioconductor.org/packages/release/bioc/vignettes/motifStack/inst/doc/motifStack_HTML.html
## plot the logo stack with heatmap.
library("MotifDb")
matrix.fly <- query(MotifDb, "Dmelanogaster")
motifs2 <- as.list(matrix.fly)
## use data from FlyFactorSurvey
motifs2 <- motifs2[grepl("Dmelanogaster\\-FlyFactorSurvey\\-",
                         names(motifs2))]
## format the names
names(motifs2) <- gsub("Dmelanogaster_FlyFactorSurvey_", "",
                       gsub("_FBgn\\d+$", "",
                            gsub("[^a-zA-Z0-9]","_",
                                 gsub("(_\\d+)+$", "", names(motifs2)))))
motifs2 <- motifs2[unique(names(motifs2))]
pfms <- sample(motifs2, 30)
pfms <- mapply(names(pfms), pfms, FUN=function(.ele, .pfm){
  new("pfm",mat=.pfm, name=.ele)}
  ,SIMPLIFY = FALSE)

df <- data.frame(A=runif(n = 30), B=runif(n = 30), C=runif(n = 30), D=runif(n = 30))
map2col <- function(x, pal){
  rg <- range(x) # return min and max
  pal[findInterval(x, seq(rg[1], rg[2], length.out = length(pal)+1), 
                   all.inside = TRUE)]
}
dl <- lapply(df, map2col, pal=heat.colors(10))
## alignment of the pfms, this step will make the motif logos occupy 
## more space. Users can skip this alignment to see the difference.
pfmsAligned <- DNAmotifAlignment(pfms)
library(ade4)
hc <- clusterMotifs(pfms) 
phylog <- ade4::hclust2phylog(hc)

library(RColorBrewer)
color <- brewer.pal(10, "Set3")

motifSig <- motifSignature(pfms, phylog, cutoffPval=0.0001, min.freq=1)
gpCol <- sigColor(motifSig)
## plot motifs
motifPiles(phylog=phylog, 
           pfms=pfms
           # col.tree=rep(color, each=5),
           # col.leaves=rep(rev(color), each=5),
           # col.pfms2=gpCol,
           # r.anno=rep(0.02, length(dl)),
           # col.anno=dl,
           # motifScale="logarithmic",
           # plotIndex=TRUE,
           # groupDistance=10
           )

###############################
# Kbest 
kbest_names <- sapply(kbest_motifs, function(x) x['name'])
kbestMotifs <- sapply(kbest_motifs, function(x) x['motif'])
# make the list with two vectors
kbest_list <- stats::setNames(as.list(kbestMotifs), kbest_names)
kbest_list2 <- mapply(names(kbest_list), kbest_list, FUN=function(.ele, .pfm){
  new("pfm",mat=.pfm, name=.ele)}
  ,SIMPLIFY = FALSE)

Kbest_hc <- clusterMotifs(kbest_list) 
kbest_phylog <- ade4::hclust2phylog(Kbest_hc)

# Randomfroest 
rf_names <- sapply(RF_motifs, function(x) x['name'])
rfMotifs <- sapply(RF_motifs, function(x) x['motif'])
# make the list with two vectors
rf_list <- stats::setNames(as.list(rfMotifs), rf_names)
rf_list2 <- mapply(names(rf_list), rf_list, FUN=function(.ele, .pfm){
  new("pfm",mat=.pfm, name=.ele)}
  ,SIMPLIFY = FALSE)

rf_hc <- clusterMotifs(rf_list) 
rf_phylog <- ade4::hclust2phylog(rf_hc)


### motiPiles to plot
motifPiles(phylog=kbest_phylog,
           pfms=kbest_list2
           # col.tree=rep(color, each=5),
           # col.leaves=rep(rev(color), each=5),
           # col.pfms2=gpCol,
           # r.anno=rep(0.02, length(dl)),
           # col.anno=dl,
           # motifScale="logarithmic",
           # plotIndex=TRUE,
           # groupDistance=10
)

motifPiles(phylog=rf_phylog, 
           pfms=rf_list2
           # col.tree=rep(color, each=5),
           # col.leaves=rep(rev(color), each=5),
           # col.pfms2=gpCol,
           # r.anno=rep(0.02, length(dl)),
           # col.anno=dl,
           # motifScale="logarithmic",
           # plotIndex=TRUE,
           # groupDistance=10
)
############use PCC to plot
# to compute the PCC between different motifs
PCCmatrix <- function(kmerList){
  kmers <- c()
  for (kmer in kmerList){
    kmers <- c(kmers,create_motif(kmer,name=kmer,pseudocount = 1, nsites = 50))
  }
  return(compare_motifs(kmers, method='PCC',min.mean.ic = 0,tryRC=TRUE,
                        score.strat = "a.mean"))
}


kbest_pcc <- PCCmatrix(kbest_names)
tiff('pcc.tiff')
pheatmap::pheatmap(kbest_pcc,cluster_cols=F,
                   cluster_rows=F)
dev.off()
# library(ade4)
kbest_hc <- as.dist(1-kbest_pcc) %>% hclust(method = 'average') %>% hclust2phylog # ade4::hclust2phylog
tiff('../Results/kbest_motifs_pcc_cluster.tiff', res=300, height = 3000, width = 2000)
motifPiles(phylog=kbest_hc, 
           pfms=kbest_list2,
           labels.nodes=""
           # col.tree=rep(color, each=5),
           # col.leaves=rep(rev(color), each=5),
           # col.pfms2=gpCol,
           # r.anno=rep(0.02, length(dl)),
           # col.anno=dl,
           # motifScale="logarithmic"
           # plotIndex=TRUE,
           # groupDistance=10
)
dev.off()

##RF
rf_pcc <- PCCmatrix(rf_names)
# library(ade4)
rf_hc <- as.dist(1-rf_pcc) %>% hclust(method = 'average') %>% hclust2phylog # ade4::hclust2phylog
tiff('../Results/rf_motifs_pcc_cluster.tiff', res=300, height = 3000, width = 2000)
motifPiles(phylog=rf_hc, 
           pfms=rf_list2,
           # labels.nodes=""
           # col.tree=rep(color, each=5),
           # col.leaves=rep(rev(color), each=5),
           # col.pfms2=gpCol,
           # r.anno=rep(0.02, length(dl)),
           # col.anno=dl,
           # motifScale="logarithmic"
           # plotIndex=TRUE,
           # groupDistance=10
)
dev.off()

################ Find motif patterns

# install.packages('tcR')
library(tcR)
get.all.substrings('ATCGCT', .min.len=4)

kbest_names
rf_names
kbest_patterns <- data.frame(Substring=0,Start=0,End=0)
rf_patterns <- data.frame(Substring=0,Start=0,End=0)
for (i in kbest_names) kbest_patterns <- rbind(kbest_patterns,tcR::get.all.substrings(i,.min.len = 4))
for (i in rf_names) rf_patterns <- rbind(rf_patterns,tcR::get.all.substrings(i,.min.len = 4))

kbest_patterns$method <- 'kbest'
rf_patterns$method <- 'rf'
allpatterns <- rbind(kbest_patterns[-1,],rf_patterns[-1,])
tiff('kbest_motif_substring.tiff')
kbest_patterns %>% dplyr::slice(2:n()) %>% ggplot(aes(x=Substring))+geom_bar() + 
  coord_flip()+
  theme(plot.title = element_text(hjust = 0.5))+
  ggtitle('kbest_substrings')
dev.off()

tiff('RF_motif_substring.tiff')
rf_patterns %>% dplyr::slice(2:n()) %>% ggplot(aes(x=Substring))+geom_bar() + 
  coord_flip()+
  theme(plot.title = element_text(hjust = 0.5))+
  ggtitle('RandomForest_substrings')
dev.off()