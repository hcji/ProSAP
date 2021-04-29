# -*- coding: utf-8 -*-
"""
Created on Tue Apr 27 09:47:47 2021

@author: hcji
"""


import numpy as np

from PyQt5.QtCore import Qt
from PyQt5 import QtCore, QtGui, QtWidgets
from Utils import TableModel, fit_curve


class CurveFitThread(QtCore.QThread):

    _ind = QtCore.pyqtSignal(str)
    _res = QtCore.pyqtSignal(list)
 
    def __init__(self, prots, temps, data_1, data_2, minR2, maxPlateau, h_axis):
        super(CurveFitThread, self).__init__()
        self.prots = prots
        self.temps = temps
        self.data_1 = data_1
        self.data_2 = data_2
        self.minR2 = minR2
        self.maxPlateau = maxPlateau
        self.h_axis = h_axis
        self.working = True
 
    def __del__(self):
        self.wait()
        self.working = False
 
    def run(self):
        for i, p in enumerate(self.prots):
            x = self.temps
            y1 = np.array(self.data_1[self.data_1.iloc[:,0] == p].iloc[0,1:])
            y2 = np.array(self.data_2[self.data_2.iloc[:,0] == p].iloc[0,1:])
            rv = fit_curve(x, y1, y2, self.minR2, self.maxPlateau, self.h_axis)
            
            self._ind.emit(str(int(100 * (i+1) / len(self.prots))))
            self._res.emit(list(rv))
        # print('finished')


class ROCThread(QtCore.QThread):

    _ind = QtCore.pyqtSignal(str)
    _res = QtCore.pyqtSignal(float)
 
    def __init__(self, prot, data, dist, proteinPair):
        super(ROCThread, self).__init__()
        self.prot = prot
        self.dist = dist
        self.proteinPair = proteinPair
        self.data = data
        self.working = True
 
    def __del__(self):
        self.wait()
        self.working = False
 
    def run(self):
        for i in self.proteinPair.index:
            p1 = self.proteinPair['Protein A'][i]
            p2 = self.proteinPair['Protein B'][i]
            a = np.where(self.prot == p1)[0]
            b = np.where(self.prot == p2)[0]
            if (len(a) > 0) and (len(b) > 0):
                d = self.dist[a[0], b[0]]
                self._ind.emit(str(int( 100 * (i + 1) / len(self.proteinPair))))
                self._res.emit(d)
        self._ind.emit(str(int(100)))
        # print('finished')
        

class PairThread(QtCore.QThread):

    _ind = QtCore.pyqtSignal(str)
    _res = QtCore.pyqtSignal(list)
 
    def __init__(self, prot1, dist1, prot2, dist2, proteinPair, n):
        super(PairThread, self).__init__()
        self.prot1 = prot1
        self.dist1 = dist1
        self.prot2 = prot2
        self.dist2 = dist2
        self.proteinPair = proteinPair
        self.n = n
        self.working = True
 
    def __del__(self):
        self.wait()
        self.working = False
 
    def run(self):
        all_prot = np.intersect1d(self.prot1, self.prot2)
        w1, w2 = [], []
        for i, p in enumerate(all_prot):
            w1.append(np.where(self.prot1 == p)[0][0])
            w2.append(np.where(self.prot2 == p)[0][0])
            self._ind.emit(str(int(100 * (i + 1) / len(self.proteinPair))))
        
        a = np.random.choice(np.arange(len(w1)), self.n)
        b = np.random.choice(np.arange(len(w1)), self.n)
        m1, n1 = np.array(w1)[a], np.array(w1)[b]
        m2, n2 = np.array(w2)[a], np.array(w2)[b]
        negDist = np.abs(self.dist1[m1, n1] - self.dist2[m2, n2])
        
        for i in self.proteinPair.index:
            p1 = self.proteinPair['Protein A'][i]
            p2 = self.proteinPair['Protein B'][i]
            a1 = np.where(self.prot1 == p1)[0]
            b1 = np.where(self.prot1 == p2)[0]
            a2 = np.where(self.prot2 == p1)[0]
            b2 = np.where(self.prot2 == p2)[0]
            if (len(a1) > 0) and (len(b1) > 0):
                if (len(a2) > 0) and (len(b2) > 0):
                    d1 = round(self.dist1[a1[0], b1[0]], 3)
                    d2 = round(self.dist2[a2[0], b2[0]], 3)
                    d = round(abs(d1 - d2), 3)
                    p = round(1 - len(np.where(negDist < d)[0]) / len(negDist), 3)
                    self._ind.emit(str(50 + int( 50 * (i + 1) / len(self.proteinPair))))
                    self._res.emit([d, p, d1, d2])
        # print('finished')



class ComplexThread(QtCore.QThread):

    _ind = QtCore.pyqtSignal(str)
    _res = QtCore.pyqtSignal(list)
 
    def __init__(self, prot1, dist1, prot2, dist2, proteinComplex, n):
        super(ComplexThread, self).__init__()
        self.prot1 = prot1
        self.dist1 = dist1
        self.prot2 = prot2
        self.dist2 = dist2
        self.proteinComplex = proteinComplex
        self.n = n
        self.working = True
 
    def __del__(self):
        self.wait()
        self.working = False
 
    def run(self):
        pass
