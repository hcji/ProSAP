# -*- coding: utf-8 -*-
"""
Created on Mon Nov 29 08:43:58 2021

@author: jihon
"""

import numpy as np
import pandas as pd

R_isavailable = False

r_codes = '''

if (! 'BiocManager' %in% installed.packages()){
  suppressMessages(install.packages('BiocManager'))
}

if (!'TPP2D' %in% installed.packages()){
  suppressMessages(BiocManager::install('TPP2D'))
}

suppressMessages(library(TPP2D))

estimate_TPP2D <- function(data, maxit=200, B=20, alpha=0.1){
  suppressMessages({
    data <- resolveAmbiguousProteinNames(data)
    params_df <- getModelParamsDf(data, maxit = maxit)
    fstat_df <- computeFStatFromParams(params_df)
    
    null_df <- bootstrapNullAlternativeModel(
      df = data, params_df = params_df, 
      maxit = maxit, B = B,
      verbose = FALSE)
    
    fdr_df <- getFDR(df_out = fstat_df,
                     df_null = null_df,
                     squeezeDenominator = TRUE)
  })
  return(fdr_df)
}

findHits <- findHits

'''

try:
    from rpy2 import robjects
    from rpy2.robjects import numpy2ri, pandas2ri

    numpy2ri.activate()
    pandas2ri.activate()

    robjects.r(r_codes)

    estimate_TPP2D = robjects.globalenv['estimate_TPP2D']
    findHits = robjects.globalenv['findHits']
    R_isavailable = True

except: 
    pass



class TPP2D:
    
    def __init__(self, data):
        self.data = data
        self.fdr_df = None
        
    def check(self):
        X = self.data
        cols = ['representative', 'clustername', 'experiment', 'temperature', 'conc', 'raw_value', 'rel_value']
        for c in cols:
            if c not in X.columns:
                return False
        
        X.loc[:,'log_conc'] = np.log10(X.loc[:,'conc'])
        X.loc[:,'log2_value'] = np.log2(X.loc[:,'rel_value'])
        self.data = X
        return True
        
        
    def fit_data(self, maxit = 200, B = 20):
        
        if not R_isavailable:
            return None
        
        else:
            try:
                self.fdr_df = pd.DataFrame(estimate_TPP2D(self.data, maxit, B))
                self.fdr_df = self.fdr_df[self.fdr_df['dataset'] == 'true']
                return self.fdr_df
            except:
                return pd.DataFrame()
            
            
    def find_hits(self, alpha = 0.1):
        if self.fdr_df is None:
            return None
        else:
            return pd.DataFrame(findHits(self.fdr_df, alpha))
        
        