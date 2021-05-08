import numpy as np
import pandas as pd

from scipy.stats import ttest_ind
from rpy2 import robjects
from rpy2.robjects import numpy2ri, pandas2ri

numpy2ri.activate()
pandas2ri.activate()

r_codes = '''

if (! 'BiocManager' %in% installed.packages()){
  install.packages('BiocManager')
}

if (!'limma' %in% installed.packages()){
  BiocManager::install('limma')
}

library(limma)

p_value_adjust <- function(pvals, method = 'BH'){
  adjust_pvals = p.adjust(as.numeric(pvals), method = method)
  return(adjust_pvals)
}


do_limma <- function(X, y, names){
  X <- as.matrix(X)
  lbs <- unique(y)
  y[y==lbs[1]] <- 'G1'
  y[y==lbs[2]] <- 'G2'
  rownames(X) <- names
  design <- model.matrix(~0+factor(y))
  colnames(design) <- levels(factor(y))
  contrast.matrix <- makeContrasts(G2-G1,levels = c('G1', 'G2'))
  fit <- lmFit(X, design)
  fit2 <- contrasts.fit(fit, contrast.matrix)
  fit2 <- eBayes(fit2)
  DEG <- topTable(fit2, coef=1, n=Inf)
  return(DEG)
}
'''

robjects.r(r_codes)

do_limma = robjects.globalenv['do_limma']
p_value_adjust = robjects.globalenv['p_value_adjust']


def do_ttest(X, y):
    i = np.where(y == 'Case')[0]
    j = np.where(y == 'Control')[0]
    pval = []
    for k in range(X.shape[0]):
        pval.append(ttest_ind(X.iloc[k, i], X.iloc[k, j]).pvalue)
    return np.expand_dims(np.array(pval), 1)


class iTSA:
    
    def __init__(self, method = 't-Test'):
        self.method = method
        
        
    def fit_data(self, X, y, names):
        self.X = X
        self.y = np.array(y)
        self.names = names
        self.lbs = np.unique(y)
        
        if self.method == 'Limma':
            res = do_limma(self.X, self.y, self.names)
            colnames = np.array(res.colnames)
            res = pd.DataFrame(np.array(res).T)
            for i in range(1, res.shape[1]):
                if i == 1:
                    res.iloc[:,i] = np.round(res.iloc[:,i].astype(float), 4)
                else:
                    res.iloc[:,i] = -np.round(np.log10(res.iloc[:,i].astype(float)), 4)
            res.columns = colnames
            res = res[['ID', 'logFC', 'P.Value', 'adj.P.Val']]
            res.columns = ['Accession', 'logFC', '-logPval', '-logAdjPval']
        
        elif self.method == 't-Test':
            i = np.where(self.y == self.lbs[0])[0]
            j = np.where(self.y == self.lbs[1])[0]
            pval = []
            for k in range(X.shape[0]):
                pval.append(ttest_ind(X.iloc[k, i], X.iloc[k, j]).pvalue)
            logfc = self.fold_change(self.X, self.y, self.lbs)
            adjpval = np.array(p_value_adjust(pval))
            lgpval = -np.round(np.log10(np.array(pval)),4)
            lgapval = -np.round(np.log10(np.array(adjpval)),4)
            res = pd.DataFrame({'Accession': self.names, 'logFC': logfc, '-logPval': lgpval, '-logAdjPval': lgapval})
            
        else:
            raise IOError('{} is not a support method'.format(self.method))
        X['Accession'] = names
        res = pd.merge(res, X, on='Accession')
        res = res.sort_values(by = '-logAdjPval', ascending=False)
        return res
        
        
    def fold_change(self, X, y, lbs):
        case_val = np.nanmean(X.loc[:, y == lbs[0]], axis = 1)
        cont_val = np.nanmean(X.loc[:, y == lbs[1]], axis = 1)
        return np.round(case_val - cont_val, 4)
    