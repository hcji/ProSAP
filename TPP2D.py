# -*- coding: utf-8 -*-
"""
Created on Mon Nov 29 08:43:58 2021

@author: jihon
"""

import numpy as np
import pandas as pd

R_isavailable = False

# data = pd.read_csv('data/TPP2D/panobinostat_tpp2d_extract.csv')

r_codes = '''

if (! 'BiocManager' %in% installed.packages()){
  suppressMessages(install.packages('BiocManager'))
}

if (!'TPP2D' %in% installed.packages()){
  suppressMessages(BiocManager::install('TPP2D'))
}

suppressMessages(library(TPP2D))

estimate_TPP2D <- function(data, maxit=200, B=20){
  suppressMessages({
    data <- resolveAmbiguousProteinNames(data)
    data <- recomputeSignalFromRatios(data)
    params_df <- getModelParamsDf(data, maxit = 500)
    fstat_df <- computeFStatFromParams(params_df)
    
    null_df <- bootstrapNullAlternativeModel(
      df = data, params_df = params_df, 
      maxit = 500, B = 20,
      verbose = FALSE)
    
    fdr_df <- getFDR(df_out = fstat_df,
                     df_null = null_df,
                     squeezeDenominator = TRUE)
    
    hits_df <- findHits(fdr_df, alpha = 0.1)
  })
  return(fdr_df)
}

'''

try:
    from rpy2 import robjects
    from rpy2.robjects import numpy2ri, pandas2ri

    numpy2ri.activate()
    pandas2ri.activate()

    robjects.r(r_codes)

    estimate_TPP2D = robjects.globalenv['estimate_TPP2D']
    R_isavailable = True

except: 
    pass



class TPP2D:
    
    def __init__(self, data):
        self.data = data
        
        
    def check(self):
        X = self.data
        for col in ['log_conc', 'log2_value']:
            X.loc[:,col] = pd.to_numeric(X.loc[:,col], errors='coerce')
        self.data = X
        
        
    def fit_data(self, maxit = 200, B = 20):
        
        if not R_isavailable:
            return None
        
        else:
            try:
                return estimate_TPP2D(self.data, maxit, B)
            except:
                return np.nan
        
        
        