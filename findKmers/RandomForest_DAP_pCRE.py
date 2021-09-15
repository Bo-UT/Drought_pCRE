# %%
import numpy as np
import pandas as pd

from sklearn.model_selection import cross_val_score
from bayes_opt import BayesianOptimization
from sklearn.ensemble import RandomForestClassifier

from sklearn.feature_selection import SelectKBest
from sklearn.feature_selection import chi2

# %%
# RandomForest paremeter fine tuning
def RandomForestClassifierBayesianOptimization(X_dtrain, y_dtrain):
    def function(n_estimators, max_depth, min_samples_split, min_samples_leaf,max_features):
        return cross_val_score(
               RandomForestClassifier(
                   n_estimators=int(round(max(n_estimators,0))),                                                               
                   max_depth=int(round(max(max_depth,1))),
                   min_samples_split=int(round(max(min_samples_split,2))),
                   min_samples_leaf = int(round(max(min_samples_leaf,1))),
                   max_features = int(round(max(max_features,1))),
                   random_state=42),  
               X=X_dtrain, 
               y=y_dtrain, 
               cv=10,
               scoring='f1').mean()

    parameters = {"n_estimators": (10, 500),
                  "max_depth": (5, 100),
                  "min_samples_split": (2,10),
                 'min_samples_leaf':(1,10),
                 'max_features':(1,10)}
    #optimize with BayesianOptimization
    bayes = BayesianOptimization(function, parameters)
    bayes.maximize(init_points=3, n_iter=7)
    optimal_params = bayes.max
    return optimal_params['params']

# %%
def RF(X,y):
    parameters = RandomForestClassifierBayesianOptimization(X,y)
    # convert float parameter to int
    parameters = {key:round(value) for (key,value) in parameters.items()}
    # create a model
    model = RandomForestClassifier(**parameters)
    # cross validation to get f1, aucroc and importance score
    f1 = cross_val_score(model, X=X, y=y, cv=10, scoring = 'f1').mean()
    auroc = cross_val_score(model, X=X, y=y, cv=10, scoring = 'roc_auc').mean()
    importance = model.fit(X,y).feature_importances_
    
    return f1,auroc,importance.tolist()

# %% 
def feature_selection(X, y):
    bestfeatures = SelectKBest(score_func=chi2, k=X.shape[1])
    fit = bestfeatures.fit(X,y)
    return fit.scores_.tolist()

# %%
def main(k):
    # load cluster kth kmer matrix
    cur_kmer = kmer_dir+'df_cluster_{}.csv'.format(k)
    df_posCRE = pd.read_csv(cur_kmer, index_col=0)
    df_posCRE['target'] = 1

    # load cluster kth DAP matrix
    cur_DAP = DAP_dir+'cluster_{}.txt'.format(k)
    df_DAP_tmp = pd.read_table(cur_DAP)
    df_posDAP = peakMatrix.loc[df_posCRE.index.tolist(), df_DAP_tmp.TF]
    df_posDAP['target'] = 1

    CRE_scores = open('../results/metricScores/cluster_%s_CREscores.txt'%k,'w')
    CRE_scores.write('repeat\tf1\tauroc')
    df_CRE_RF_importance = pd.DataFrame(columns = df_posCRE.columns[:-1])
    df_CRE_Kbest_importance = pd.DataFrame(columns = df_posCRE.columns[:-1])

    DAP_scores = open('../results/metricScores/cluster_%s_DAPscores.txt'%k,'w')
    DAP_scores.write('repeat\tf1\tauroc')
    df_DAP_RF_importance = pd.DataFrame(columns = df_posDAP.columns[:-1])
    df_DAP_Kbest_importance = pd.DataFrame(columns = df_posDAP.columns[:-1])

    DAP_CRE_scores = open('../results/metricScores/cluster_%s_DAP_CREscores.txt'%k,'w')
    DAP_CRE_scores.write('repeat\tf1\tauroc')
    df_DAP_CRE_RF_importance = pd.DataFrame(columns = df_posCRE.columns[:-1].tolist()+df_posDAP.columns[:-1].tolist())
    df_DAP_CRE_Kbest_importance = pd.DataFrame(columns = df_posCRE.columns[:-1].tolist()+df_posDAP.columns[:-1].tolist())
    
    for i in range(rep_num):
        # sample control CRE data
        df_sampledCtrlCRE = df_ctrlCRE.sample(df_posCRE.shape[0])
        df_sampledCtrlCRE = df_sampledCtrlCRE.loc[:, df_posCRE.columns[:-1]]
        df_sampledCtrlCRE['target'] = 0

        # sample control DAP data 
        df_sampledCtrlDAP = peakMatrix.loc[df_sampledCtrlCRE.index.tolist(), df_posDAP.columns[:-1]]
        df_sampledCtrlDAP['target'] = 0

        # combine and shuffle the data for models training
        df_dataCRE = pd.concat([df_posCRE, df_sampledCtrlCRE], axis=0).sample(frac=1)
        df_dataDAP = pd.concat([df_posDAP, df_sampledCtrlDAP], axis=0).sample(frac=1)

        df_positive = pd.concat([df_posCRE.drop(columns='target'), df_posDAP], axis=1)
        df_control = pd.concat([df_sampledCtrlCRE.drop(columns='target'),df_sampledCtrlDAP], axis=1)
        df_dataCRE_DAP = pd.concat([df_positive, df_control], axis=0).sample(frac=1)

        # model I: DAP
        X_DAP = df_dataDAP.drop(columns='target')
        y_DAP = df_dataDAP['target']
        f1, auroc, importance_RF = RF(X_DAP, y_DAP)
        DAP_scores.write('\n{}\t{}\t{}'.format(i, f1, auroc))
        chi2_scores = feature_selection(X_DAP, y_DAP)
        df_DAP_RF_importance.loc[i] = importance_RF
        df_DAP_Kbest_importance.loc[i] = chi2_scores

        # model II: pCRE
        X_CRE = df_dataCRE.drop(columns='target')
        y_CRE = df_dataCRE['target']
        f1, auroc, importance_RF = RF(X_CRE, y_CRE)
        CRE_scores.write('\n{}\t{}\t{}'.format(i, f1, auroc))
        chi2_scores = feature_selection(X_CRE, y_CRE)
        df_CRE_RF_importance.loc[i] = importance_RF
        df_CRE_Kbest_importance.loc[i] = chi2_scores

        # model III: DAP+pCRE
        X_DAP_CRE = df_dataCRE_DAP.drop(columns='target')
        y_DAP_CRE = df_dataCRE_DAP['target']
        f1, auroc, importance_RF = RF(X_DAP_CRE, y_DAP_CRE)
        DAP_CRE_scores.write('\n{}\t{}\t{}'.format(i, f1, auroc))
        chi2_scores = feature_selection(X_DAP_CRE, y_DAP_CRE)
        df_DAP_CRE_RF_importance.loc[i] = importance_RF
        df_DAP_CRE_Kbest_importance.loc[i] = chi2_scores

        print('*'*50+'Iteration: '+str(i)+'*'*50)
    df_CRE_RF_importance.to_csv('../results/df_CRE_RF_importance.csv')
    df_CRE_Kbest_importance.to_csv('../results/df_CRE_Kbest_importance.csv')

    df_DAP_RF_importance.to_csv('../results/df_DAP_RF_importance.csv')
    df_DAP_Kbest_importance.to_csv('../results/df_DAP_Kbest_importance.csv')

    df_DAP_CRE_RF_importance.to_csv('../results/df_DAP_CRE_RF_importance.csv')
    df_DAP_CRE_Kbest_importance.to_csv('../results/df_DAP_CRE_Kbest_importance.csv')

    print('*'*70)
    print('Finished cluster {} '.format(k))
    print('*'*70)

# %%
rep_num = 10 # model repeat times
import os
kmer_dir = '../data/kmersResults/'
DAP_dir = '../results/clusterSignificantTF/'

peakMatrix = pd.read_csv('../data/DAPpeakMatrix_inclControlGenes.csv', index_col=0)
df_ctrlCRE = pd.read_csv('../../controlGenes/results/df_ctrlGenesCREmatrix.csv', index_col=0)
if __name__ == '__main__':
    for clu in range(1,41):
        main(clu)
    



