#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Feb 24 16:34:26 2023

@author: phuasx
"""

import pandas as pd
import numpy as np
from sklearn.decomposition import PCA
from sklearn import preprocessing
import matplotlib.pyplot as plt
import scipy.stats as stats
import seaborn as sns


def pca_df(df):
    """
    
    Parameters
    ----------
    df : pd.DataFrame
        p x n matrix, where columns n = samples and index p = protein 

    Returns
    -------
    pca_df : pd.DataFrame
        PCA matrix with n indices of samples with n columns of PCs.

    """
    from sklearn.decomposition import PCA
    from sklearn import preprocessing
    
    scaled_df = preprocessing.scale(df.T) 
    pca = PCA()
    pca.fit(scaled_df)
    pca_data = pca.transform(scaled_df)
    per_var = np.round(pca.explained_variance_ratio_* 100, decimals = 1)
    labels = ['PC' + str(x) for x in range(1, len(per_var)+1)]
    heads = list(df.columns)
    pca_df = pd.DataFrame(pca_data, index = heads, columns=labels)
    print('\t\t'.join([f'PC{i}' for i in range(1,11)]))
    print('\t\t'.join([str(i) for i in per_var[:10]]))
    return pca_df

def shannon_entropy(x):
    return -np.log2(x)

def pca_association_test(matrix, variable):
    """
    
    Parameters
    ----------
    matrix : pd.DataFrame
        p x n matrix, where columns n = samples and index p = protein 
    variable : pd.Series, list or iterable
        Vector with n number of variables to be tested

    Returns
    -------
    out : pd.Series
        Vector with n Kruskal Wallis H test p-values, indexed with each PC.

    """
    import scipy.stats as stats
    pcadf = pca_df(matrix)
    variable = list(variable)
    unique = list(set(variable))
    
    pval_all = []
    
    for pc in pcadf.columns:
        pc_series = pcadf.loc[:,pc]
        kw_input = [pcadf.iloc\
                    [[ind for (ind, j) in enumerate(variable) if j == i]]\
                        .loc[:,pc].tolist() for i in unique]
        pval = stats.kruskal(*kw_input).pvalue
        pval_all.append(pval)
        
    out = pd.Series(pval_all, index = pcadf.columns)
    return out


# columns = batch_sampleID_runID_class
# indices = log2 normalised protein intensity

file_list = ['ecl.csv', 'egl.csv', 'gcl.csv', 'ggl.csv', 
             'gr.csv', 'grp_jelena.csv', 'rcl.csv', 'rgl.csv']

ldict = {'gr.csv':'Uncorrected protein data',
         'grp_jelena.csv': 'Williams ' + '$\it{et}$' + ' ' + '$\it{al.}$'+ '(2021)',
         'rcl.csv': 'Class-specific protein level',
         'rgl.csv': 'Global protein level',
         'ecl.csv': 'Class-specific peptide level',
         'egl.csv': 'Global peptide level',
         'gcl.csv': 'Class-specific ambiguous peptide guided',
         'ggl.csv': 'Global ambiguous peptide guided'}


interbatch = []
interclass = []


for file in file_list:
    #loading data
    df = pd.read_csv(file,index_col = 0) # read file
    df = df.fillna(0) # zero impute if any missing value
    
    batch_variable = [i.split('_')[0] for i in df.columns] # batch variables
    class_variable = [i.split('_')[-1] for i in df.columns] # class variables
    
    batch_pval = pca_association_test(df, batch_variable) # get batch pvalue series
    class_pval = pca_association_test(df, class_variable) # get class pvalue series
    
    interbatch.append(batch_pval) # append to output variable
    interclass.append(class_pval) # append to output variable


ib_df = pd.concat(interbatch, axis = 1) # concatenate all outputs
ib_df.columns = file_list
ib_df = ib_df

    
ic_df = pd.concat(interclass, axis = 1)
ic_df.columns = file_list
ic_df = ic_df.apply(shannon_entropy)


def ib_plotter(unique, title, cm, ib_df = ib_df):
    sns.set('notebook',{'figure.figsize':(20,10),'figure.dpi':600})
    sns.set_style('white') # set style to white background without grid
    fig, ax = plt.subplots(figsize = (5,2), dpi = 600) # create subplot
    count = 0 # enumerate for colour
    plotter = ib_df.loc[ib_df.index[:10]] # to first ten PC
    plotter.index = [i.split('PC')[-1] for i in plotter.index] # change PC to number
    # for each correction method
    for file in unique: 
        # plot
        ax.plot(plotter[file],linewidth = 2,marker = 'o',c = cm[count],markersize = 6,label = ldict[file])
        count += 1
    # set title
    plt.title('Interbatch KW Test',fontsize = 10)
    # plot parameters start
    plt.xticks(fontsize = 10) 
    plt.yticks([i/4 for i in (range(0,5))]+[0.05],fontsize = 10)
    plt.axhline(0.05, linewidth = 2,linestyle = '--',c = 'grey',alpha = 0.5)
    pv = '$\it{p}$' + '-value'
    plt.ylabel(pv,fontsize = 10)
    plt.rcParams['legend.title_fontsize'] = 20
    plt.tight_layout()
    # plot parameters end
    plt.savefig(title,dpi=600) # save figure


def ic_plotter(unique, title, cm, ib_df = ic_df):
    sns.set('notebook',{'figure.figsize':(20,10),'figure.dpi':600})
    sns.set_style('white') # set style to white background without grid
    fig, ax = plt.subplots(figsize = (5,2), dpi = 600) 
    count = 0 # enumerate for color
    plotter = ic_df.loc[ic_df.index[:10]] # to first ten PC
    plotter.index = [i.split('PC')[-1] for i in plotter.index] # change PC to number
    # for each correction method
    for file in unique:
        ax.plot(plotter[file],linewidth = 2,marker = 'o',c = cm[count],markersize = 6,label = ldict[file])
        count += 1
    # set title
    plt.title('Interclass KW Test',fontsize = 10)
    # plot parameters start
    plt.xticks(fontsize = 10)
    plt.axhline(0.05, linewidth = 2,linestyle = '--',c = 'grey',alpha = 0.5)
    l10 = '$\mathregular{-log_2}$' + '('
    pv = l10+ '$\it{p}$' + '-value' + ')'
    plt.ylabel(pv,fontsize = 10)
    plt.ylim([0,150])
    plt.rcParams['legend.title_fontsize'] = 20
    plt.tight_layout()
    # plot parameters end
    plt.savefig(title,dpi=600) # save figure

title = 'interbatch_kw_test_gl.svg'
unique = ['gr.csv','grp_jelena.csv','ggl.csv','egl.csv','rgl.csv']
cm = ['#f5c1cf','#f0d39c','#c8e3ba','#8cad7b','#a8d392']
ib_plotter(title = title, unique = unique, cm = cm, ib_df = ib_df)

title = 'interbatch_kw_test_cl.svg'
unique = ['gr.csv','grp_jelena.csv','gcl.csv','ecl.csv','rcl.csv']
cm = ['#f5c1cf','#f0d39c','#e0bbe4','#957dad','#d291bc']
ib_plotter(title = title, unique = unique, cm = cm, ib_df = ib_df)

title = 'interclass_kw_test_gl'
unique = ['gr.csv','grp_jelena.csv','gcl.csv','ecl.csv','rcl.csv']
cm = ['#f5c1cf','#f0d39c','#e0bbe4','#957dad','#d291bc']
ic_plotter(unique, title, cm, ib_df = ic_df)

title = 'interclass_kw_test_cl'
unique = ['gr.csv','grp_jelena.csv','ggl.csv','egl.csv','rgl.csv']
cm = ['#f5c1cf','#f0d39c','#c8e3ba','#8cad7b','#a8d392']
ic_plotter(unique, title, cm, ib_df = ic_df)
