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

if (!'edgeR' %in% installed.packages()){
  BiocManager::install('edgeR')
}

if (!'DESeq2' %in% installed.packages()){
  BiocManager::install('DESeq2')
}

library(limma)
library(edgeR)
library(DESeq2)


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
  DEG$ID <- rownames(DEG)
  return(DEG)
}


do_edgeR <- function(X, y, names){
  lbs <- unique(y)
  y[y==lbs[1]] <- 'G1'
  y[y==lbs[2]] <- 'G2'
  
  DGElist <- DGEList( counts = X, group = y)
  design <- model.matrix(~0+y)
  
  rownames(design) = colnames(data)
  colnames(design) <- levels(y)
  
  DGElist <- calcNormFactors( DGElist )
  DGElist <- estimateGLMCommonDisp(DGElist, design)
  DGElist <- estimateGLMTrendedDisp(DGElist, design)
  DGElist <- estimateGLMTagwiseDisp(DGElist, design)
  
  fit <- glmFit(DGElist, design)
  results <- glmLRT(fit, contrast = c(-1, 1)) 
  nrDEG_edgeR <- topTags(results, n = nrow(DGElist))
  nrDEG_edgeR <- as.data.frame(nrDEG_edgeR)
  
  return(nrDEG_edgeR)
}


do_DESeq2 <- function(X, y, names){
  lbs <- unique(y)
  y[y==lbs[1]] <- 'G1'
  y[y==lbs[2]] <- 'G2'
  
  data <- X
  # DESeq2 need the max vaule lower than 2147483647
  if (max(X) > 2147483647){
    r <- round(max(X) / 2147483647 + 1)
  } else {
    r <- 1
  }
  data <- round(data / r)
  condition <- factor(y)
  coldata <- data.frame(row.names = colnames(X), condition)
  dds <- DESeqDataSetFromMatrix(countData = data,
                                colData = coldata,
                                design = ~condition)
  dds$condition<- relevel(dds$condition, ref = "G1")
  dds <- DESeq(dds)
  dds <- as.data.frame(results(dds))
  return(dds)
}
     
      
estimate_df <- function(rss1, rssDiff){
  rm_idx <- is.na(rssDiff) | (rssDiff <= 0)
  rss1 <- rss1[!rm_idx]
  rssDiff <- rssDiff[!rm_idx]

  M = median(rssDiff, na.rm = TRUE)
  V = mad(rssDiff, na.rm = TRUE)^2
  s0_sq = 1/2 * V/M
  rssDiff = rssDiff/s0_sq
  rss1 = rss1/s0_sq
  d1 = MASS::fitdistr(x = rssDiff, densfun = "chi-squared", start = list(df = 1), method = "Brent", lower = 0, upper = length(rssDiff))[["estimate"]]
  d2 = MASS::fitdistr(x = rss1, densfun = "chi-squared", start = list(df = 1),  method = "Brent", lower = 0, upper = length(rssDiff))[["estimate"]]

  out <- c(d1 = d1, d2 = d2, s0_sq = s0_sq)
  return(out)
}

'''


robjects.r(r_codes)

do_limma = robjects.globalenv['do_limma']
do_DESeq2 = robjects.globalenv['do_DESeq2']
do_edgeR = robjects.globalenv['do_edgeR']
estimate_df = robjects.globalenv['estimate_df']
p_value_adjust = robjects.globalenv['p_value_adjust']


class iTSA:
    
    def __init__(self, method = 't-Test'):
        self.method = method
        
        
    def fit_data(self, X, y, names):
        self.X = X
        self.y = np.array(y)
        self.names = np.array(names)
        self.lbs = np.unique(y)
        
        if self.method == 'Limma':
            res = do_limma(np.log2(self.X), self.y, self.names)
            res = pd.DataFrame(res)
            res = res[['ID', 'logFC', 'P.Value', 'adj.P.Val']]
            for i, c in enumerate(res.columns):
                if c == 'ID':
                    pass
                elif c == 'logFC':
                    res.iloc[:,i] = np.round(res.iloc[:,i].astype(float), 4)
                else:
                    res.iloc[:,i] = -np.round(np.log10(res.iloc[:,i].astype(float)), 4)
            
            res.columns = ['Accession', 'logFC', '-logPval', '-logAdjPval']         
            
        elif self.method == 'edgeR':
            res = do_edgeR(self.X, self.y, self.names)
            res = pd.DataFrame(res)
            res['ID'] = self.names
            res = res[['ID', 'logFC', 'PValue', 'FDR']]
            res.columns = ['Accession', 'logFC', '-logPval', '-logAdjPval'] 
            res['-logPval'] = -np.round(np.log10(res['-logPval']), 4)
            res['-logAdjPval'] = -np.round(np.log10(res['-logAdjPval']), 4)
            res['logFC'] = np.round(res['logFC'], 4)
            
        elif self.method == 'DESeq2':
            res = do_DESeq2(self.X, self.y, self.names)
            res = pd.DataFrame(res)
            res['ID'] = self.names
            res = res[['ID', 'log2FoldChange', 'pvalue', 'padj']]
            res.columns = ['Accession', 'logFC', '-logPval', '-logAdjPval'] 
            res['-logPval'] = -np.round(np.log10(res['-logPval']), 4)
            res['-logAdjPval'] = -np.round(np.log10(res['-logAdjPval']), 4)
            res['logFC'] = np.round(res['logFC'], 4)
        
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
        X_ = X.copy()
        X_['Accession'] = names
        res = pd.merge(res, X_)
        res = res.sort_values(by = '-logAdjPval', ascending=False)
        return res
        
        
    def fold_change(self, X, y, lbs):
        case_val = np.nanmean(X.loc[:, y == lbs[1]], axis = 1)
        cont_val = np.nanmean(X.loc[:, y == lbs[0]], axis = 1)
        return np.round(np.log2(case_val / cont_val), 4)
    