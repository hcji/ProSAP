# -*- coding: utf-8 -*-
"""
Created on Wed Apr 14 10:54:02 2021

@author: hcji
"""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

from scipy.stats import ttest_ind
from scipy.sparse import csc_matrix, eye, diags
from scipy.sparse.linalg import spsolve
from scipy.misc import derivative

import os
# from joblib import Parallel, delayed
from scipy.optimize import curve_fit
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
    try: 
        paras1 = curve_fit(meltCurve, x, y1, bounds=(0, [float('inf'), float('inf'), 1.5]))[0]
        paras2 = curve_fit(meltCurve, x, y2, bounds=(0, [float('inf'), float('inf'), 1.5]))[0]    
        yh1 = meltCurve(x, paras1[0], paras1[1], paras1[2])
        yh2 = meltCurve(x, paras2[0], paras2[1], paras2[2])
        rss1 = np.sum((y1 - yh1) ** 2) + np.sum((y2 - yh2) ** 2)
        
        x0 = np.concatenate([x, x])
        y0 = np.concatenate([y1, y2])
        paras0 = curve_fit(meltCurve, x0, y0, bounds=(0, [float('inf'), float('inf'), 1.5]))[0]
        yh0 = meltCurve(x0, paras0[0], paras0[1], paras0[2])
        rss0 = np.sum((y0 - yh0) ** 2)
        
        r1 = max(r2_score(y1, yh1), 0)
        r2 = max(r2_score(y2, yh2), 0)
        R = min(r1, r2)
    except:
        rss1, rss0, R = 0, 0, 0
    diff = abs(rss1 - rss0)
    return rss0, rss1, diff, R


def fit_curve(x, y1, y2, minR2 = 0.8, maxPlateau = 0.3, h_axis = 0.5):
    
    paras1 = curve_fit(meltCurve, x, y1, bounds=(0, [float('inf'), float('inf'), maxPlateau]))[0]
    paras2 = curve_fit(meltCurve, x, y2, bounds=(0, [float('inf'), float('inf'), maxPlateau]))[0]
        
    yh1 = meltCurve(x, paras1[0], paras1[1], paras1[2])
    yh2 = meltCurve(x, paras2[0], paras2[1], paras2[2])
    
    r1 = max(r2_score(y1, yh1), 0)
    r2 = max(r2_score(y2, yh2), 0)

    x1 = np.arange(x[0], x[-1], 0.01)
    yy1 = meltCurve(x1, paras1[0], paras1[1], paras1[2])
    yy2 = meltCurve(x1, paras2[0], paras2[1], paras2[2])
    Tm1 = x1[np.argmin(np.abs(yy1 - h_axis))]
    Tm2 = x1[np.argmin(np.abs(yy2 - h_axis))]
    
    sl1 = np.min((yy1[1:] - yy1[:-1]) / 0.01)
    sl2 = np.min((yy2[1:] - yy2[:-1]) / 0.01)
    sl = min(sl1, sl2)
    
    if min(r1, r2) < minR2:
        deltaTm = np.nan
    elif max(np.min(yh1), np.min(yh2)) > h_axis:
        deltaTm = np.nan
    else:
        deltaTm = Tm2 - Tm1
        
    return r1, r2, Tm1, Tm2, deltaTm, sl


def WhittakerSmooth(x, lambda_, differences=1):
    X = np.matrix(x)
    w = np.ones(x.shape[0])
    m = X.size
    E = eye(m, format='csc')
    D = E[1:] - E[:-1]
    W = diags(w, 0, shape=(m, m))
    A = csc_matrix(W + (lambda_ * D.T * D))
    B = csc_matrix(W * X.T)
    background = spsolve(A, B)
    return np.array(background)
