# %%
from scipy.stats import fisher_exact
import numpy as np
import pandas as pd
import re
import copy
import os
# %%
# calcualte the p-value with Fisher exact test
def enrichedTF(k, df_pos, df_neg):
    file = open('clusterSignificantTF/cluster_{}.txt'.format(k),'w')
    file.write('TF\tp-value')
    # try except ValueError
    fet_pvalue = {}
    for col in df_pos.columns:
        try:
            TP = df_pos[col].sum()                   #Positive Examples with kmer present
            FP = df_neg[col].sum()                   #Negative Examples with kmer present
            TN = df_neg[df_neg[col]==0].shape[0]     #Negative Examples without kmer
            FN = df_pos[df_pos[col]==0].shape[0]     #Positive Examples without kmer
            print([[TP,FN],[FP,TN]])
            oddsratio,pvalue = fisher_exact([[TP,FN],[FP,TN]],alternative='greater') # one tail Fisher exact test
            fet_pvalue[col] = pvalue
        except ValueError:
            fet_pvalue[col] = 1.0
        if pvalue<0.01:
            file.write('\n{}\t{}'.format(col, pvalue))
    file.close()
# %%
# load the peak matrix that includes control genes
peakMatrix = pd.read_csv('DAPpeakMatrix_inclControlGenes.csv', index_col=0)
# make control genes peak matrix
ctrl = pd.read_csv('df_controlGenes.csv', index_col=0)
peakMatrix_ctrl = peakMatrix.loc[ctrl.index,:]
# load differentially expressed genes 
deg = pd.read_csv('../kmers/dehydration_allDE_cluster.csv', index_col=0)
cluster_num = len(set(deg.cluster))
for i in range(1,cluster_num+1):
    posList = deg[deg.cluster==i].index
    df_pos = peakMatrix.loc[posList,:]
    enrichedTF(i, df_pos, peakMatrix_ctrl)