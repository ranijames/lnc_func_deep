#!/usr/bin/env python
# coding: utf-8

import pandas as pd
from scipy.stats import pearsonr
import argparse
from argparse import ArgumentParser, SUPPRESS
from matplotlib import pyplot as plt
import sys
from datetime import datetime
startTime = datetime.now()
#main()

__author__ = 'Alva James'
"""
useage example: python correlation_analysis.py -ag /home/alva/ukbio_phenotypes/derived/projects/autoencoder_lncRNAs/data/All_brca.pkl -a /home/alva/ukbio_phenotypes/derived/projects/autoencoder_lncRNAs/data/gene_annoataion.txt -o /home/alva/ukbio_phenotypes/derived/projects/autoencoder_lncRNAs/data/correlation_analysis
"""
parser = argparse.ArgumentParser()
parser.add_argument('-ag',help='the input dataframe', dest='all_genes')
parser.add_argument('-a', help='GTF edited annotation file', dest='annotation_file')
parser.add_argument('-o', help='Path you need the output', dest='output_path')
args = parser.parse_args()
 
if len(sys.argv) < 3:
    parser.print_help()
    sys.exit(0)
all_genes         = args.all_genes
annotation        = args.annotation_file
output            = args.output_path


# In[ ]:


def data_prep(all_genes,annotation ):
    """
    Data preparation for the correlation analysis
    """
    PC      = pd.DataFrame()
    lncRNAs = pd.DataFrame()
    all_genes           = pd.read_pickle(all_genes)
    all_genes.reset_index(inplace=True)
    all_genes['Gene_id'] = all_genes['Gene_id'].apply(lambda x:x.split('.')[0])
    annotation           = pd.read_csv(annotation,sep='\t',header=0)
    all_genes_annotation = pd.merge(all_genes, annotation[["Gene_id","Biotype"]], on ='Gene_id')
    all_genes_annotation["Biotype"]= all_genes_annotation["Biotype"].str.strip()
    PC                             = PC.append(all_genes_annotation.query('Biotype == "protein_coding"'))
    lncRNAs                        = lncRNAs.append(all_genes_annotation.query('Biotype == "lncRNA"'))
    lncRNAs['Gene_id']             = lncRNAs['Gene_id']  + '_lncRNAs'
    PC['Gene_id']                  = PC['Gene_id']  + '_PC'
    lncRNA_PC                      = pd.concat([PC,lncRNAs])
    del lncRNA_PC['Biotype']
    lncRNA_PC_T=lncRNA_PC.set_index("Gene_id").T
    #lncRNA_PC_T.rename_axis("Samples", inplace =True)
    return (lncRNA_PC_T)


# In[ ]:


def correlation_analysis(lncRNA_PC_T):
    """
    Function for correlation analysis
    """
    correlations = pd.DataFrame()
    for PC in [column for column in lncRNA_PC_T.columns if '_PC' in column]: 
        for lncRNA in [column for column in lncRNA_PC_T.columns if '_lncRNAs' in column]:
                    correlations = correlations.append(pd.Series(pearsonr(lncRNA_PC_T[PC],lncRNA_PC_T[lncRNA]),index=['PCC', 'p-value'],name=PC + '_' +lncRNA))
    correlations.reset_index(inplace=True)
    correlations.rename(columns={0:'name'},inplace=True)
    correlations['PC']= correlations['index'].apply(lambda x:x.split('PC')[0])
    correlations['lncRNAs']=correlations['index'].apply(lambda x:x.split('PC')[1])
    correlations['lncRNAs'] =correlations['lncRNAs'].apply(lambda x:x.split('_')[1])
    correlations['PC']=correlations.PC.str.strip('_')
    correlations.drop('index',axis=1,inplace=True)
    correlations=correlations.reindex(columns=['PC','lncRNAs','PCC','p-value'])         
    return(correlations)


#all_genes  ='/home/alva/ukbio_phenotypes/derived/projects/autoencoder_lncRNAs/data/All_brca.pkl'
#annotation ='/home/alva/ukbio_phenotypes/derived/projects/autoencoder_lncRNAs/data/gene_annoataion.txt'
#output     = '/home/alva/ukbio_phenotypes/derived/projects/autoencoder_lncRNAs/data/correlation_analysis'
lncRNA_PC_T  = data_prep(all_genes,annotation )
correlations = correlation_analysis(lncRNA_PC_T)
correlations.to_pickle(output+'/'+'BRCA_lncRNAs_PC_correlation.pkl')
print ("***Sucessfully finished correlation analysis ************")
print("--- Total time taken to finish the correlation analysis:", datetime.now() - startTime)


