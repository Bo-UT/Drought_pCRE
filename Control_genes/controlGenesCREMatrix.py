# %%
import pandas as pd
import numpy as np
import os
from plotnine import *
from Bio import SeqIO
from Bio.Seq import Seq

# %%
# get all kmers from in positive genes
kmer_dir = '../../DAP_pCRE_RF/data/kmersResults'
allkmers = []
for file in os.listdir(kmer_dir)[:40]:
    filetoload = os.path.join(kmer_dir, file)
    df = pd.read_csv(filetoload, index_col=0)
    allkmers += df.columns.tolist()
allkmers = list(set(allkmers))
kmerstosave = open('../results/allkmers.txt','w')
for i in allkmers:
    kmerstosave.write('\n%s'%i)
kmerstosave.close()
# %%
# find significant kmers for cotrol genes
def find_kmers(afastafile, km):
    df = pd.DataFrame(columns=km,dtype=int)
    for seq_record in SeqIO.parse(afastafile, 'fasta'):
        seq = str(seq_record.seq)
        genes_array = []       
        # iterate each CRE
        ### Not consider occurence number
        for ki in km:
            kmer = Seq(ki)
            if str(kmer) in seq or str(kmer.reverse_complement()) in seq: # find all occurence of kmer in the promoter including reverse complement
                genes_array.append(1)
            else:
                genes_array.append(0)
        # add counting to the dataframe 
        df.loc[seq_record.id] = genes_array
    return df

# %%
with open ('../data/negativeFasta') as negativeFasta:
    df_controlGenesCRE = find_kmers(negativeFasta, allkmers)
    df_controlGenesCRE.to_csv('../results/df_ctrlGenesCREmatrix.csv')
# %%
