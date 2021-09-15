library(tidyverse)
library(reshape2)
library(ComplexHeatmap)
library(MotifDb)
library(universalmotif)
library(ggseqlogo)

rm(list=ls())

setwd('E:/Study/Courses/Projects/ResearchProjects/TimeSeriesDrought_PEG/ML_pipeline')
setwd('D:/Projects/ENAP1/Seed germination/Document')
##################################################
########### Heatmap of pCRE ranks########################
##################################################
plot_rank <- function(afile){
  # read data
  cluster <- read.csv(paste0('RandomForestResults/',afile),row.name=1)
  clu_num <- strsplit(strsplit(afile,split='_rf')[[1]][1],split='_')[[1]][2]
# get rank for each pCRE
getRank <- function(x){rank(x,ties.method='average')}
# make the heatmap
heatmap_cluster <- cluster %>% apply(.,1,getRank) %>% t() %>% 
  Heatmap(name = "Rank", 
          show_row_names = FALSE,
          row_title='Repeat',
          column_title = paste0('Cluster ', clu_num),
          column_names_rot=90,
          column_names_gp = gpar(fontsize = 4))
# save the plot
tiff(paste0('similarTF/',strsplit(afile,split='_rf')[[1]][1],'_heatmap.tiff'),
     units="in", width=20, height=10, res=400)
print(heatmap_cluster)
dev.off()
# extract left most 10 pCRE (top 10)
top10 <- colnames(cluster)[column_order(heatmap_cluster)][1:10] # row_order(heatmap_clu1)
return(top10)
}


##################################################
########### Compare motifs########################
##################################################
# iniput a pCRE, output top
compare_motif <- function(cre){
  cre <- create_motif(cre, pseudocount = 1, nsites = 50)
  res <- compare_motifs(c(motifs_ara,cre),method = 'PCC',min.mean.ic = 0,tryRC=TRUE,
                        score.strat = "a.mean")
  TF_top3 <- res[dim(res)[1],] %>% # take the last row that compares cre with all other TF binding sites
    as.data.frame() %>% 
    rownames_to_column('TF') %>% 
    mutate(TF=sapply(TF,function(x) unlist(strsplit(x,split=' '))[1])) %>% # split TF names with [duplicated..]
    stats::setNames(c('TF','pcc')) %>% 
    mutate(pcc=round(pcc,4)) %>% 
    distinct() %>% # remove duplicated
    arrange(desc(pcc)) %>% # sort in descending order
    dplyr::slice(1:12) #extract top 4, the first one is the CRE itself
  return(TF_top3) 
}

# test <- compare_motif('ATGCGT')
# test[test$TF=='motif',1] <- 'ATGCGT'
# plot_logo(test)
##################################################
########### Plot top 3 similar TF########################
##################################################
# input dataframe from compare_motif function
plot_logo <- function(adata.frame){ 
  TF_name=adata.frame$TF
  TF_pcc=adata.frame$pcc
  
  p <- function(name,pcc) {
    
    motif <- filter_motifs(motifs_ara,name = name)
    if (length(motif)==0) {
      motif_ppm <- create_motif(name)['motif'] # the first is the pCRE itself
    } else{
      motif_ppm <- motif[[1]]['motif']
    }
    # plot the motif
    pp <- ggplot()+geom_logo( motif_ppm, method = 'bits' ) + 
      annotate("text", x = 1.5, y = 1.6, label = paste('PCC=',pcc),size=3)+
      theme_logo()+theme(plot.title = element_text(hjust = 0.5))+
      ggtitle(name)
    return(pp)}
  len <- seq(length(TF_name))
  p_list <- lapply(len,function(x) p(TF_name[x],TF_pcc[x]))
  do.call(gridExtra::grid.arrange, c(p_list, ncol=3))
  # p1 <- p(TF_name[1],TF_pcc[1])
  # p2 <- p(TF_name[2],TF_pcc[2])
  # p3 <- p(TF_name[3],TF_pcc[3])
  # p4 <- p(TF_name[4],TF_pcc[4])
  # p5 <- p(TF_name[5],TF_pcc[5])
  # p6 <- p(TF_name[6],TF_pcc[6])
  
# place plots in multiple rows
  # ppp <- gridExtra::grid.arrange(p1,p2,p3,p4,p5,p6,nrow=4)
  
  # return(ppp)
}
##################################################
########### Main Function ########################
##################################################

main <- function(afile){

  top10CRE <- plot_rank(afile)
  # assign plot results
  for (k in 1:length(top10CRE)){
    top3SimilarTF <- compare_motif(top10CRE[k])
    top3SimilarTF[top3SimilarTF['TF']=='motif',1]=top10CRE[k]
    print(top3SimilarTF)
    # save plot
    plot_name <- paste0(strsplit(file,split='_rf')[[1]][1],'_',top10CRE[k],'_',k,'.tiff')
    tiff(paste0('similarTF/',plot_name),units="in", width=14, height=10, res=300)
    plot_logo(top3SimilarTF)
    dev.off()
}}
##################################################
########### Running ########################
##################################################
# Access Arabidopsis TF binding motifs 
motifs <- convert_motifs(MotifDb)
motifs_ara <- filter_motifs(motifs, organism = "Athaliana",
                            extrainfo = c("dataSource" = "jaspar2018"))

importanceFiles <- list.files('RandomForestResults')[grep('.importance',list.files('RandomForestResults'))]
for (file in importanceFiles) {
  main(file)
}




