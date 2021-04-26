# -*- coding: utf-8 -*-
"""
Created on Wed Apr 14 10:54:02 2021

@author: hcji
"""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

import os
from scipy import stats
from joblib import Parallel, delayed
from scipy.optimize import curve_fit, fsolve
from sklearn.metrics import r2_score

from PyQt5 import QtCore, QtGui, QtWidgets
from PyQt5.QtCore import Qt
'''
from rpy2 import robjects
from rpy2.robjects import numpy2ri, pandas2ri

from R_functions import R_functions

numpy2ri.activate()
pandas2ri.activate()

R = R_functions()

robjects.r(R['TPPTR_Rstring'])
runTPPTR = robjects.globalenv['runTPPTR']
'''

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
    
    
def meltCurve(T, a, b, pl):
    A = 1 - pl
    B = 1 + np.exp(b - a/T)
    return A / B + pl


def analTSA(proteinData1, proteinData2, columns, minR2 = 0.8, maxPlateau = 0.3, h_axis = 0.5):
    
    # columns = ['T37', 'T40', 'T43', 'T46', 'T49', 'T52', 'T55', 'T58', 'T61', 'T64']
    # proteinData1 = pd.read_csv('data/TPCA_TableS14_MTX.csv')
    # proteinData2 = pd.read_csv('data/TPCA_TableS14_DMSO.csv')
    
    temps = np.array([float(t.replace('T', '')) for t in columns])
    cols = ['Accession'] + columns
    data_1 = proteinData1.loc[:, cols]
    data_2 = proteinData2.loc[:, cols]

    prots = np.intersect1d(list(data_1.iloc[:,0]), list(data_2.iloc[:,0]))
    
    def fit_curve(p):
        x = temps
        y1 = np.array(data_1[data_1.iloc[:,0] == p].iloc[0,1:])
        y2 = np.array(data_2[data_2.iloc[:,0] == p].iloc[0,1:])
        
        paras1 = curve_fit(meltCurve, x, y1, bounds=(0, [float('inf'), float('inf'), maxPlateau]))[0]
        paras2 = curve_fit(meltCurve, x, y2, bounds=(0, [float('inf'), float('inf'), maxPlateau]))[0]
        
        yh1 = meltCurve(x, paras1[0], paras1[1], paras1[2])
        yh2 = meltCurve(x, paras2[0], paras2[1], paras2[2])
        
        r1 = max(r2_score(y1, yh1), 0)
        r2 = max(r2_score(y2, yh2), 0)
        
        i1 = x[ np.argmin(np.abs(y1 - h_axis))]
        i2 = x[ np.argmin(np.abs(y2 - h_axis))]
        
        Tm1 = fsolve(lambda x: meltCurve(x, paras1[0], paras1[1], paras1[2]) - h_axis, i1)[0]
        Tm2 = fsolve(lambda x: meltCurve(x, paras2[0], paras2[1], paras2[2]) - h_axis, i2)[0]
        
        Tm1 = min(Tm1, max(x))
        Tm1 = max(Tm1, min(x))
        Tm2 = min(Tm2, max(x))
        Tm2 = max(Tm2, min(x))
        
        if min(r1, r2) < minR2:
            deltaTm = 0
        else:
            deltaTm = abs(Tm1 - Tm2)
        
        return r1, r2, Tm1, Tm2, deltaTm
        
    n_core = os.cpu_count() - 1
    res = Parallel(n_jobs=n_core)(delayed(fit_curve)(p) for p in prots)
    res = pd.DataFrame(res)
    res.columns = ['Group1_R2', 'Group2_R2', 'Group1_Tm', 'Group2_Tm', 'delta_Tm']
    
    delta_Tm = res['delta_Tm']
    p_Val = []
    for i in range(len(res)):
        s = delta_Tm[i]
        pv = stats.t.sf((s - np.mean(delta_Tm)) / np.std(delta_Tm), len(delta_Tm)-1)
        p_Val.append(pv)
    score = -np.log10(np.array(p_Val)) * (res['Group1_R2'] * res['Group2_R2']) ** 2
    
    res['Accession'] = prots
    res['delta_Tm'] = delta_Tm
    res['p_Val (-log10)'] = -np.log10(p_Val)
    res['Score'] = score
    res = np.round(res, 3)
    
    res = res[['Accession', 'Score', 'p_Val (-log10)', 'delta_Tm', 'Group1_R2', 'Group2_R2', 'Group1_Tm', 'Group2_Tm']]
    res = res.sort_values(by = 'Score',axis = 0, ascending = False)
    
    '''
    result = runTPPTR(data_1, data_2, temps, Pl, a, b, minR2, maxPlateau)
    
    colName = np.array(result.colnames)
    output = pd.DataFrame(np.array(result).T)
    output.columns = colName
    '''
    return res