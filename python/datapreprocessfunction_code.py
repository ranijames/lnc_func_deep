import pandas as pd
import numpy as np

"""
Global functions for data pre processing and data clean up
"""
def make_list(df,string):
    lists   =list(df[string].drop_duplicates())
    return lists

def annon(df1,string,df2):
    df_annon     =pd.merge(df1,df2,on=string)
    return df_annon

def add_column(df,lists,string,newcolname):
    df[newcolname]  = np.where(df['gene_sym'].isin(lists),'Yes','No')
    return df

def readfile(path,filename):
    df = pd.read_csv(path+filename,header=0,sep='\t',skipinitialspace=True)
    return df

def annotate_features(df,cancer_driver,drugtarget,triplex,promo_methy,SNPS):
    """
    Feature annotation, adding column to de lncRNAs dataframe based on the feature information
    """
    df =add_column(df,cancer_driver,"gene_sym","cancer_driver")
    df =add_column(df,drugtarget,"gene_sym","drug_target")
    df =add_column(df,triplex,"gene_sym","triplex")
    df =add_column(df,promo_methy,"gene_sym","promoter_methylated")
    df =add_column(df,SNPS,"gene_sym","SNPS")
    return df  

def rename_cols(df,oldname,newname):
    df.rename(columns={oldname:newname},inplace=True)
    return df

def outers_join(df1,df2,string):
    df_new = pd.merge(df1, df2, on=string,how="outer")
    return df_new

def label_paths (row):
    if row['Pathway'] == "unknown":
        return 0
    else:
        return 1
    
def percentage_calc(count_ofyes,countofno,question):
    """
    Percentage calculator for the number of pathways
    """
    pct_of_no_sub = countofno/(countofno+count_ofyes)
    pct_of_sub    = (count_ofyes)/(countofno+count_ofyes)
    return(print("percentage of lncRNAs without" + question, pct_of_no_sub*100),print("percentage of lncRNAs with"+question, pct_of_sub*100))

def function_process(df,string):
    """
    Spiliting the rows based on comma for each interactions/pathways
    """
    df1 = df.reindex(df.index.repeat(df[string].fillna("").str.split(',').apply(len)))
    df1 = df1.drop([string],1)
    s   = df[string].str.split(',', expand=True).stack().reset_index(level=1,drop=True)
    # Now grouping the series and df using cumcount.
    df1 = df1.set_index(df1.groupby(df1.index).cumcount(), append=True)
    s   = s.to_frame('intersections').set_index(s.groupby(s.index).cumcount(), append=True)
    df1 = df1.join(s, how='outer').reset_index(level=[0,1],drop=True)
    return df1
