# -*- coding: utf-8 -*-
"""
Created on Wed Apr 14 10:54:02 2021

@author: hcji
"""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.stats import ttest_ind

import os
# import multiprocessing
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

# multiprocessing.freeze_support()

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


def fit_np(x, y1, y2):
    paras1 = curve_fit(meltCurve, x, y1, bounds=(0, [10**4, 10**4, 1.5]))[0]
    paras2 = curve_fit(meltCurve, x, y2, bounds=(0, [10**4, 10**4, 1.5]))[0]    
    yh1 = meltCurve(x, paras1[0], paras1[1], paras1[2])
    yh2 = meltCurve(x, paras2[0], paras2[1], paras2[2])
    rss1 = np.sum((y1 - yh1) ** 2) + np.sum((y2 - yh2) ** 2)
    
    x0 = np.concatenate([x, x])
    y0 = np.concatenate([y1, y2])
    paras0 = curve_fit(meltCurve, x0, y0, bounds=(0, [10**4, 10**4, 1.5]))[0]
    yh0 = meltCurve(x0, paras0[0], paras0[1], paras0[2])
    rss0 = np.sum((y0 - yh0) ** 2)
    
    diff = abs(rss1 - rss0)
    return rss0, rss1, diff


def fit_curve(x, y1, y2, minR2 = 0.8, maxPlateau = 0.3, h_axis = 0.5):
    
    paras1 = curve_fit(meltCurve, x, y1, bounds=(0, [10**4, 10**4, maxPlateau]))[0]
    paras2 = curve_fit(meltCurve, x, y2, bounds=(0, [10**4, 10**4, maxPlateau]))[0]
        
    yh1 = meltCurve(x, paras1[0], paras1[1], paras1[2])
    yh2 = meltCurve(x, paras2[0], paras2[1], paras2[2])
    
    r1 = max(r2_score(y1, yh1), 0)
    r2 = max(r2_score(y2, yh2), 0)
    '''
    i1 = x[np.argmin(np.abs(y1 - h_axis))]
    i2 = x[np.argmin(np.abs(y2 - h_axis))]
    Tm1 = fsolve(lambda x: meltCurve(x, paras1[0], paras1[1], paras1[2]) - h_axis, i1)[0]
    Tm2 = fsolve(lambda x: meltCurve(x, paras2[0], paras2[1], paras2[2]) - h_axis, i2)[0]
    Tm1 = min(Tm1, max(x))
    Tm1 = max(Tm1, min(x))
    Tm2 = min(Tm2, max(x))
    Tm2 = max(Tm2, min(x))
    '''
    x1 = np.arange(x[0], x[-1], 0.01)
    Tm1 = x1[np.argmin(np.abs(meltCurve(x1, paras1[0], paras1[1], paras1[2]) - h_axis))]
    Tm2 = x1[np.argmin(np.abs(meltCurve(x1, paras2[0], paras2[1], paras2[2]) - h_axis))]
        
    if min(r1, r2) < minR2:
        deltaTm = 0
    else:
        deltaTm = abs(Tm1 - Tm2)
        
    return r1, r2, Tm1, Tm2, deltaTm

