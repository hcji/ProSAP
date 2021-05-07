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


do_limma <- function(X, y){
  design <- model.matrix(~0+factor(y))
  colnames(design) <- levels(factor(y))
  rownames(design) <- colnames(X)
  contrast.matrix <- makeContrasts(Case-Control,levels = design)
  fit <- lmFit(X, design)
  fit2 <- contrasts.fit(fit, contrast.matrix)
  fit2 <- eBayes(fit2)
  p.value <- fit2$p.value
  return(p.value)
}