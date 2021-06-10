# -*- coding: utf-8 -*-
"""
Created on Tue Apr 27 09:47:47 2021

@author: hcji
"""


import numpy as np
import pandas as pd
from scipy.stats import norm
from scipy.spatial.distance import pdist
from scipy.optimize import curve_fit
from sklearn.metrics import r2_score

from PyQt5.QtCore import Qt
from PyQt5 import QtCore, QtGui, QtWidgets
from Utils import TableModel, fit_curve, meltCurve, fit_np
from iTSA import estimate_df


class PreprocessThread(QtCore.QThread):

    _ind = QtCore.pyqtSignal(str)
    _val = QtCore.pyqtSignal(list)
    _rsd = QtCore.pyqtSignal(float)
    _prot = QtCore.pyqtSignal(str)
 
    def __init__(self, data, psm_column, psm_thres, std_thres, columns, fun):
        super(PreprocessThread, self).__init__()
        self.data = data
        self.psm_column = psm_column
        self.psm_thres = psm_thres
        self.std_thres = std_thres
        self.columns = columns
        self.fun = fun
        self.working = True
 
    def __del__(self):
        self.wait()
        self.working = False
 
    def run(self):
        data = self.data
        for i in data.index:
            if self.psm_column == 'None':
                pass
            else:  
                wh = np.where([self.psm_column == s.split('--')[0] for s in data.columns])[0]
                psm = np.nanmean(data.iloc[i, wh].values.astype(float))
                if psm < self.psm_thres:
                    continue

            prot = data.loc[i, 'Accession']
            vals = []
            for c in self.columns:
                wh = np.where([c == s.split('--')[0] for s in data.columns])[0]
                v = data.iloc[i, wh].values.astype(float)               
                v = np.round(v, 4)
                std = np.nanstd(v) / np.nanmean(v)
                vals.append(self.fun(v))
            if std > self.std_thres:
                continue
            
            vals = np.array(vals)            
            self._ind.emit(str(int(100 * (i+1) / len(data.index))))
            self._val.emit(list(vals))
            self._prot.emit(prot)
            self._rsd.emit(std)
            
            # print(prot)
            # print(vals)
            
        self._ind.emit(str(int(100)))


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
        self._ind.emit(str(int(100)))


class NPTSAThread(QtCore.QThread):

    _ind = QtCore.pyqtSignal(str)
    _res = QtCore.pyqtSignal(list)
 
    def __init__(self, prots, temps, data_1, data_2, method):
        super(NPTSAThread, self).__init__()
        self.prots = prots
        self.temps = temps
        self.data_1 = data_1
        self.data_2 = data_2
        self.method = method
        self.working = True
 
    def __del__(self):
        self.wait()
        self.working = False
 
    def run(self):
        for i, p in enumerate(self.prots):
            x = self.temps
            y1 = np.array(self.data_1[self.data_1.iloc[:,0] == p].iloc[0,1:])
            y2 = np.array(self.data_2[self.data_2.iloc[:,0] == p].iloc[0,1:])       
            rv = fit_np(x, y1, y2, self.method)
            
            self._ind.emit(str(int(100 * (i+1) / len(self.prots))))
            self._res.emit(list(rv))
        self._ind.emit(str(int(100)))


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
            self._ind.emit(str(int(50 * (i + 1) / len(self.proteinPair))))
        
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
                else:
                    d1, d2, d, p = np.nan, np.nan, np.nan, np.nan
            else:
                d1, d2, d, p = np.nan, np.nan, np.nan, np.nan
            self._ind.emit(str(50 + int( 50 * (i + 1) / len(self.proteinPair))))
            self._res.emit([d, p, d1, d2])
        # print('finished')


class ComplexThread(QtCore.QThread):

    _ind = QtCore.pyqtSignal(str)
    _res = QtCore.pyqtSignal(list)
 
    def __init__(self, prot1, dist1, prot2, dist2, proteinComplex):
        super(ComplexThread, self).__init__()
        self.prot1 = prot1
        self.dist1 = dist1
        self.prot2 = prot2
        self.dist2 = dist2
        self.proteinComplex = proteinComplex
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
            self._ind.emit(str(int(50 * (i + 1) / len(all_prot))))
        w1 = np.array(w1)
        w2 = np.array(w2)
        
        prot_align = self.prot1[w1]
        dist1_align = self.dist1[w1,:][:, w1]
        dist2_align = self.dist2[w2,:][:, w2]
        dist_change = dist1_align - dist2_align
        
        dist1_align[dist1_align == 0] = np.nan
        dist2_align[dist2_align == 0] = np.nan
        dist_change[dist_change == 0] = np.nan
        
        dist1_mu, dist1_sigma = np.nanmean(dist1_align), np.nanstd(dist1_align)
        dist2_mu, dist2_sigma = np.nanmean(dist2_align), np.nanstd(dist2_align)
        dist_change_mu, dist_change_sigma = np.nanmean(dist_change), np.nanstd(dist_change)
        
        for i in self.proteinComplex.index:
            subunits = self.proteinComplex.loc[i, 'Subunits_UniProt_IDs']
            subunits = subunits.split(',')
            subunits = [p.replace('(', '') for p in subunits]
            subunits = [p.replace(')', '') for p in subunits]
            
            l = []
            for s in subunits:
                a = np.where(prot_align == s)[0]
                if len(a) > 0:
                    l.append(a[0])
            if len(l) <= 2:
                pm, dm, p1, d1, p2, d2 = np.nan, np.nan, np.nan, np.nan, np.nan, np.nan
            else:
                dist1_sub = dist1_align[l,:][:,l]
                dist2_sub = dist2_align[l,:][:,l]
                dist_change_sub = dist_change[l,:][:,l]
                
                d1 = round(np.nanmean(dist1_sub),3)
                d2 = round(np.nanmean(dist2_sub),3)
                dm = round(np.nanmean(dist_change_sub),3)
                
                n = len(l)
                z1 = (d1 - dist1_mu) * (n ** (1/2)) / dist1_sigma
                z2 = (d2 - dist2_mu) * (n ** (1/2)) / dist2_sigma
                zm = (dm - dist_change_mu) * (n ** (1/2)) / dist_change_sigma
        
                p1 = round(norm.sf(z1),3)
                p2 = round(norm.sf(z2),3)
                pm = round(norm.sf(zm),3)
                
            self._ind.emit(str(50 + int( 50 * (i + 1) / len(self.proteinComplex))))
            self._res.emit([n, pm, dm, p1, d1, p2, d2])        
