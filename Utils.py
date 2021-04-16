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
    


def plot_protein_complex(proteinSubunit, proteinData, colNames):
    prots = proteinSubunit.split(',')
    prots = [p.replace('(', '') for p in prots]
    prots = [p.replace(')', '') for p in prots]

    try:
        temps = np.array([float(t.replace('T', '')) for t in colNames])
    except:
        raise ValueError('invalid column names')

    pltData = dict()
    for p in prots:
        if p in list(proteinData['Accession']):
            vals = proteinData.loc[proteinData.loc[:, 'Accession'] == p, colNames]
            pltData[p] = vals.values[0,:]

    plt.figure(dpi = 300)
    for p, vec in pltData.items():
        plt.plot(temps, vec, label = p)
    plt.xlabel('temperature (℃)')
    plt.ylabel('abundances')
    plt.legend()
    plt.show()
    
    
