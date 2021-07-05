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

if (!'Inflect' %in% installed.packages()){
  BiocManager::install('Inflect')
}

if (!'writexl' %in% installed.packages()){
  BiocManager::install('writexl')
}

library(Inflect)

run_inflect <- function(Temperature, r1p1, r1p2, r2p1, r2p2, Rsq, NumSD){
    try(unlink('C:/inflect_tempdir', recursive =TRUE))
    dir.create('C:/inflect_tempdir')
    writexl::write_xlsx(r1p2, 'C:/inflect_tempdir/Condition 1.xlsx')
    writexl::write_xlsx(r2p2, 'C:/inflect_tempdir/Condition 2.xlsx')
    writexl::write_xlsx(r1p1, 'C:/inflect_tempdir/Control 1.xlsx')
    writexl::write_xlsx(r2p1, 'C:/inflect_tempdir/Control 2.xlsx')

    TRresults <- Inflect:::Inflect('C:/inflect_tempdir', Temperature, Rsq, NumSD, 2)
}
'''

robjects.r(r_codes)
run_inflect = robjects.globalenv['run_inflect']
