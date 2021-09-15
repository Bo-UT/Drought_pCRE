library(tidyverse)
library(reshape2)

setwd('E:/Study/Courses/Projects/ResearchProjects/labDehydrationProject/ML_pipeline/geneClustering/scripts')
##load data
data <- read.csv('dehydration_allDE_cluster.csv',row.names=1)
data$d0h <- 0
colnames(data)
data <- data[,c('d0h',colnames(data)[1:7])]
# data[-7] to remove original cluster column
data2 <- t(scale(t(data[-8]))) %>% as.data.frame() %>% rownames_to_column(., "GeneID") # scale data and convert row name to firt column


## kmeans cluster
num_cluster <- 40 # set cluster number
set.seed(42)
kmeans_out <- kmeans(data2[-1],centers=num_cluster,iter.max = 500,nstart=50)
## add cluster info to orig matrix 
data_with_cust_info <- data2 %>% 
mutate(cluster = paste("cluster ", kmeans_out$cluster,sep = ""))

# data_with_cust_info = cbind(data2,data$cluster) %>% stats::setNames(c(colnames(data2),'cluster'))
# dim(data_with_cust_info)
## visualise  each cluster 
options(repr.plot.width = 300, repr.plot.height = 300) # set figure size
data_with_cust_info %>% 
  # filter(cluster %in% paste0('cluster ',c(12,21,27,30,4,13,5,1,2,10,9,12))) %>% 
  gather(key = "variable" , value = "value", -c(1,9)) %>%  ### 1 is the index of column 'geneName' and 8 is the index of column 'clust'
  group_by(variable) %>%  
  mutate(cluster=factor(cluster,levels = paste0('cluster ', seq(num_cluster)))) %>% # levlels name must be same with cluster
  ggplot(aes(x =  variable , y = value , group = GeneID)) +   
  geom_point(size=1.0) +  
  geom_line(alpha = 1 , aes(col = as.character(cluster))) + 
  theme_bw() +  
  theme(legend.position = "none" , axis.text.x = element_text(angle = 90 , vjust = 0.4)) +
  labs(x='Time courses',y='Normalized expression')+
  facet_wrap(~cluster,ncol=8)+
  theme(text = element_text(size = 50),axis.title.x = element_text(margin = margin(t = 20, r = 20, b = 0, l = 0)))+
  ggsave('E:/Study/Courses/Projects/ResearchProjects/labDehydrationProject/findKmers/KmersFinding2/genecluster.tiff',
         width = 100,height = 100,units = 'cm',limitsize = FALSE)


data_with_cust_info %>% filter(cluster %in% paste0('cluster ',c(3,5,14,17,20,22,26,32,37))) %>% 
  ggplot(.)+geom_bar(aes(cluster))

for (i in c(3,5,14,17,20,22,26,32,37)) {
  data_with_cust_info %>% filter(cluster==paste0('cluster ', i)) %>% select(1) %>% 
    write.table(file = paste0("E:/Study/Courses/Projects/ResearchProjects/labDehydrationProject/findKmers/KmersFinding2/geneClusters/cluster_",i,".txt"), 
                sep = "\t",row.names = F, col.names = F,quote = F)
}

data_with_cust_info %>% select(1)
