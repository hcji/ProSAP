if (!'BiocManager' %in% installed.packages()){
  install.packages('BiocManager')
}

if (!'TPP' %in% installed.packages()){
  BiocManager::install('TPP')
}


library(TPP)

# library(readr)
# 
# proteinData1 <- read_csv('data/TPCA_TableS14_DMSO.csv')
# proteinData2 <- read_csv('data/TPCA_TableS14_MTX.csv')
# 
# data_1 <- proteinData1[,c(1, 10:19)]
# data_2 <- proteinData2[,c(1, 10:19)]


runTPPTR <- function(data_1, data_2, temps, Pl = 0, a = 550, b = 10, minR2 = 0.8, maxPlateau = 0.3){
  n <- ncol(data_1)
  colName <- colnames(data_1)[2:ncol(data_1)]
  TPP_config <- as.data.frame(matrix(NA, nrow=2, ncol=n+3))
  colnames(TPP_config) <-c('Experiment','Condition','ComparisonVT1','ComparisonVT2',colName)
  
  TPP_config[,1] <- c('group_1', 'group_2')
  TPP_config[,2] <- c('Treatment', 'Vehicle')
  
  TPP_config[,3] <- c('x', '')
  TPP_config[,4] <- c('', 'x')
  TPP_config[1, 5:ncol(TPP_config)] <- temps
  TPP_config[2, 5:ncol(TPP_config)] <- temps
  
  TPP_data <- list()
  TPP_data$group_1 <- as.data.frame(data_1)
  TPP_data$group_2 <- as.data.frame(data_2)
  
  TRresults <- analyzeTPPTR(configTable = TPP_config,
                            normalize = FALSE,
                            methods = 'meltcurvefit',
                            data = TPP_data,
                            idVar = 'Accession',
                            fcStr = '',
                            nCores = detectCores(logical = F) - 1,
                            startPars = c(Pl = Pl, a = a, b = b),
                            pValFilter = list(minR2 = minR2, maxPlateau = maxPlateau),
                            xlsxExport = FALSE,
                            plotCurves = FALSE)
  return(TRresults)
}

