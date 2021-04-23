# -*- coding: utf-8 -*-
"""
Created on Wed Apr 14 10:54:02 2021

@author: hcji
"""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

from PyQt5 import QtCore, QtGui, QtWidgets
from PyQt5.QtCore import Qt

from rpy2 import robjects
from rpy2.robjects import numpy2ri, pandas2ri

from R_functions import R_functions

numpy2ri.activate()
pandas2ri.activate()

R = R_functions()

robjects.r(R['TPPTR_Rstring'])
runTPPTR = robjects.globalenv['runTPPTR']

class TableModel(QtCore.QAbstractTableModel):
    def __init__(self, data, showAllColumn=False):
        QtCore.QAbstractTableModel.__init__(self)
        self.showAllColumn = showAllColumn
        self._data = data

    def rowCount(self, parent=None):
        return self._data.shape[0]

    def columnCount(self, parent=None):
        return self._data.shape[1]

    def data(self, index, role=Qt.DisplayRole):
        if index.isValid():
            if role == Qt.DisplayRole:
                return str(self._data.iloc[index.row(), index.column()])
        return None

    def headerData(self,col,orientation,role):
        if orientation == Qt.Horizontal and role == Qt.DisplayRole:
            if type(self._data.columns[col]) == tuple:
                return self._data.columns[col][-1]
            else:
                return self._data.columns[col]
        elif orientation == Qt.Vertical and role == Qt.DisplayRole:
            return (self._data.axes[0][col])
        return None
    

def analTSA(proteinData1, proteinData2, columns, Pl = 0, a = 550, b = 10, minR2 = 0.8, maxPlateau = 0.6):
    '''
    columns = ['T37', 'T40', 'T43', 'T46', 'T49', 'T52', 'T55', 'T58', 'T61', 'T64']
    proteinData1 = pd.read_csv('data/TPCA_TableS14_MTX.csv')
    proteinData2 = pd.read_csv('data/TPCA_TableS14_DMSO.csv')
    '''
    temps = np.array([float(t.replace('T', '')) for t in columns])
    cols = ['Accession'] + columns
    data_1 = proteinData1.loc[:, cols]
    data_2 = proteinData2.loc[:, cols]
    
    result = runTPPTR(data_1, data_2, temps, Pl, a, b, minR2, maxPlateau)
    
    colName = np.array(result.colnames)
    output = pd.DataFrame(np.array(result).T)
    output.columns = colName
    
