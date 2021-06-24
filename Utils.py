# -*- coding: utf-8 -*-
"""
Created on Wed Apr 14 10:54:02 2021

@author: hcji
"""

import numpy as np
from scipy.sparse import csc_matrix, eye, diags
from scipy.sparse.linalg import spsolve
from scipy.spatial.distance import pdist

# from joblib import Parallel, delayed
from scipy.optimize import curve_fit
from sklearn.metrics import r2_score

from PyQt5 import QtCore, QtGui, QtWidgets
from PyQt5.QtCore import Qt


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


def fit_NPARC(x, y11, y12, y21, y22, minR2_null = 0.8, minR2_alt = 0.8, maxPlateau = 0.3):
    x_null = np.concatenate([x, x])
    y1_null = np.concatenate([y11, y21])
    y2_null = np.concatenate([y12, y22])
    try:
        paras1_null = curve_fit(meltCurve, x_null, y1_null, bounds=(0, [15000, 250, maxPlateau]))[0]
        paras2_null = curve_fit(meltCurve, x_null, y2_null, bounds=(0, [15000, 250, maxPlateau]))[0]
        yh1_null = meltCurve(x_null, paras1_null[0], paras1_null[1], paras1_null[2])
        yh2_null = meltCurve(x_null, paras2_null[0], paras2_null[1], paras2_null[2])
        rss_null = np.sum((y1_null - yh1_null) ** 2) + np.sum((y2_null - yh2_null) ** 2)
        r1_null = max(r2_score(y1_null, yh1_null), 0)
        r2_null = max(r2_score(y2_null, yh2_null), 0)
        if min(r1_null, r2_null) < minR2_null:
            rss_null = np.nan
    except:
        rss_null = np.nan
        r1_null = np.nan
        r2_null = np.nan

    x_alt = np.concatenate([x, x, x, x])
    y_alt = np.concatenate([y11, y12, y21, y22])
    try:
        paras_alt = curve_fit(meltCurve, x_alt, y_alt, bounds=(0, [float('inf'), float('inf'), maxPlateau]))[0]
        yh_alt = meltCurve(x_alt, paras_alt[0], paras_alt[1], paras_alt[2])
        r_alt = max(r2_score(y_alt, yh_alt), 0)
        if r_alt < minR2_alt:
            rss_alt = np.nan
        else:
            rss_alt = np.sum((y_alt - yh_alt) ** 2)
    except:
        rss_alt = np.nan
    
    rss_diff = abs(rss_null - rss_alt)
    return r1_null, r2_null, rss_null, rss_alt, rss_diff
    

def fit_dist(x, y1, y2, method = 'cityblock', minR2 = 0.8, maxPlateau = 0.3):
    try: 
        paras1 = curve_fit(meltCurve, x, y1, bounds=(0, [15000, 250, maxPlateau]))[0]
        paras2 = curve_fit(meltCurve, x, y2, bounds=(0, [15000, 250, maxPlateau]))[0]
    except:
        return np.nan, np.nan, np.nan, np.nan, np.nan, np.nan
    
    yh1 = meltCurve(x, paras1[0], paras1[1], paras1[2])
    yh2 = meltCurve(x, paras2[0], paras2[1], paras2[2])
    rss1 = np.sum((y1 - yh1) ** 2) + np.sum((y2 - yh2) ** 2)
        
    x1 = np.arange(x[0], x[-1], 0.01)
    yy1 = meltCurve(x1, paras1[0], paras1[1], paras1[2])
    yy2 = meltCurve(x1, paras2[0], paras2[1], paras2[2])
    sl1 = np.min((yy1[1:] - yy1[:-1]) / 0.01)
    sl2 = np.min((yy2[1:] - yy2[:-1]) / 0.01)
    sl = min(sl1, sl2)
    
    rss0 = np.sum(yh1)
    rss1 = np.sum(yh2)
    diff = pdist(np.vstack([yh1, yh2]), metric = method)[0]
    
    r1 = max(r2_score(y1, yh1), 0)
    r2 = max(r2_score(y2, yh2), 0)
    if min(r1, r2) < minR2:
        diff = np.nan
        sl = np.nan
        rss0 = np.nan
        rss1 = np.nan
    return r1, r2, rss0, rss1, diff, sl 



def fit_curve(x, y1, y2, minR2 = 0.8, maxPlateau = 0.3, h_axis = 0.5):
    try:
        paras1 = curve_fit(meltCurve, x, y1, bounds=(0, [float('inf'), float('inf'), maxPlateau]))[0]
        paras2 = curve_fit(meltCurve, x, y2, bounds=(0, [float('inf'), float('inf'), maxPlateau]))[0]
    except:
        return 0, 0, 0, 0, 0, 0
    
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
    elif max(np.min(yh1), np.min(yh2)) > h_axis + 0.1:
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


def ReplicateCheck(tppTable):
    # tppTable = pd.read_csv('test_tpp.csv', index_col=0)
    if 'Rep2delta_Tm' not in tppTable.columns:
        return tppTable
    for i in tppTable.index:
        pval_1 = tppTable.loc[i, 'Rep1pVal (-log10)']
        pval_2 = tppTable.loc[i, 'Rep2pVal (-log10)']
        cond_1 = (max(pval_1, pval_2) > 1.302) and (min(pval_1, pval_2) > 1)
        
        delm_1 = tppTable.loc[i, 'Rep1delta_Tm']
        delm_2 = tppTable.loc[i, 'Rep2delta_Tm']
        cond_2 = delm_1 * delm_2 > 0
        
        tm = abs(tppTable.loc[i, 'Rep1Group1_Tm'] - tppTable.loc[i, 'Rep2Group1_Tm'])
        cond_3 = min(delm_1, delm_2) > tm

        mins_1 = tppTable.loc[i, 'Rep1min_Slope']
        mins_2 = tppTable.loc[i, 'Rep2min_Slope']
        cond_4 = max(mins_1, mins_2) < -0.06
        
        if cond_1 and cond_2 and cond_3 and cond_4:
            tppTable.loc[i, 'Score'] = 0
    return tppTable
            