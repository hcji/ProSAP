# -*- coding: utf-8 -*-
"""
Created on Fri Apr 16 16:23:29 2021

@author: hcji
"""


import matplotlib
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler
from scipy.optimize import curve_fit

matplotlib.use("Qt5Agg")
from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg as FigureCanvas
from matplotlib.figure import Figure
from PyQt5 import QtCore, QtGui, QtWidgets
from Utils import meltCurve


class MakeFigure(FigureCanvas):
    def __init__(self,width=5, height=5, dpi=300):
        self.fig = Figure(figsize=(width, height), dpi=dpi)
        self.fig.subplots_adjust(top=0.95,bottom=0.2,left=0.15,right=0.85)
        super(MakeFigure,self).__init__(self.fig) 
        self.axes = self.fig.add_subplot(111)
        self.axes.spines['bottom'].set_linewidth(0.5)
        self.axes.spines['left'].set_linewidth(0.5)
        self.axes.spines['right'].set_linewidth(0.5)
        self.axes.spines['top'].set_linewidth(0.5)
        self.axes.tick_params(labelsize=5)
        FigureCanvas.setSizePolicy(self, QtWidgets.QSizePolicy.Expanding, QtWidgets.QSizePolicy.Expanding)
        FigureCanvas.updateGeometry(self)        

    
    def ProteinComplexFigure(self, proteinSubunit, proteinData, colNames):
        self.fig.subplots_adjust(top=0.95,bottom=0.2,left=0.15,right=0.8)
        prots = proteinSubunit.split(',')
        prots = [p.replace('(', '') for p in prots]
        prots = [p.replace(')', '') for p in prots]

        try:
            temps = np.array([float(t.replace('T', '')) for t in colNames])
            temps_ = np.arange(temps[0], temps[-1], 0.1)
        except:
            raise ValueError('invalid column names')

        pltData = dict()
        for p in prots:
            if p in list(proteinData['Accession']):
                vals = proteinData.loc[proteinData.loc[:, 'Accession'] == p, colNames]
                pltData[p] = vals.values[0,:]
        
        self.axes.cla()
        for p, vec in pltData.items():
            paras = curve_fit(meltCurve, temps, vec, bounds=(0, [float('inf'), float('inf'), 0.3]))[0]
            self.axes.scatter(temps, vec, marker='.', label = p, s = 15)
            self.axes.plot(temps_, meltCurve(temps_, paras[0], paras[1], paras[2]), lw = 1)
            
        self.axes.set_xlabel('Temperature (℃)', fontsize=5)
        self.axes.set_ylabel('Abundances', fontsize=5)
        self.axes.legend(fontsize=4, bbox_to_anchor=(1,1), loc="upper left")
        self.draw()
    
    
    def AverageTSAFigure(self, proteinData1, proteinData2, colNames):
        try:
            temps = np.array([float(t.replace('T', '')) for t in colNames])
        except:
            raise ValueError('invalid column names')
        
        vec_1 = np.mean(proteinData1.loc[:, colNames], axis = 0)
        vec_2 = np.mean(proteinData2.loc[:, colNames], axis = 0)
        
        self.axes.cla()
        self.axes.plot(temps, vec_1, label = 'Group 1')
        self.axes.plot(temps, vec_2, label = 'Group 2')
        
        self.axes.tick_params(labelsize=4)
        self.axes.set_xlabel('Temperature (℃)', fontsize=3)
        self.axes.set_ylabel('Abundances', fontsize=3)
        self.axes.legend(fontsize=3)    
        self.draw()
        
    
    def SingleTSAFigure(self, proteinData1, proteinData2, colNames, ProteinAccession):
        try:
            temps = np.array([float(t.replace('T', '')) for t in colNames])
            temps_ = np.arange(temps[0], temps[-1], 0.1)
        except:
            raise ValueError('invalid column names')
            
        vec_1 = proteinData1.loc[proteinData1.loc[:, 'Accession'] == ProteinAccession, colNames].values[0,:]
        vec_2 = proteinData2.loc[proteinData2.loc[:, 'Accession'] == ProteinAccession, colNames].values[0,:]
        
        paras1 = curve_fit(meltCurve, temps, vec_1, bounds=(0, [float('inf'), float('inf'), 0.3]))[0]
        paras2 = curve_fit(meltCurve, temps, vec_2, bounds=(0, [float('inf'), float('inf'), 0.3]))[0]
        
        self.axes.cla()
        self.axes.scatter(temps, vec_1, marker='.', label = 'Group 1_{}'.format(ProteinAccession), s = 10)
        self.axes.scatter(temps, vec_2, marker='.', label = 'Group 2_{}'.format(ProteinAccession), s = 10)
        self.axes.plot(temps_, meltCurve(temps_, paras1[0], paras1[1], paras1[2]))
        self.axes.plot(temps_, meltCurve(temps_, paras2[0], paras2[1], paras2[2]))
        
        self.axes.tick_params(labelsize=4)
        self.axes.set_xlabel('Temperature (℃)', fontsize=4)
        self.axes.set_ylabel('Abundances', fontsize=4)
        self.axes.legend(fontsize=3)
        self.draw()
    
    
    def RSDHistFigure(self, rsdList):
        self.axes.cla()
        self.axes.hist(rsdList, 100)
        self.axes.tick_params(labelsize=3)
        self.axes.set_xlabel('RSD', fontsize=3)
        self.axes.set_ylabel('Number', fontsize=3)
        self.draw()
    
    
    def ROCFigure(self, fpr, tpr, auroc):
        self.axes.cla()
        self.axes.plot(fpr, tpr, label='AUC = {}'.format(auroc), color = 'red', lw=0.7)
        self.axes.plot([0, 1], [0, 1], color='black', linestyle='--', lw=0.7)
        self.axes.set_xlabel('False Positive Rate', fontsize = 4)
        self.axes.set_ylabel('True Positive Rate', fontsize = 4)
        self.axes.tick_params(labelsize=4)
        self.axes.legend(fontsize=3)
        self.draw()
        
        
    def ProteinPairFigure(self, p1, p2, proteinData, colNames):
        prots = [p1, p2]
        try:
            temps = np.array([float(t.replace('T', '')) for t in colNames])
            temps_ = np.arange(temps[0], temps[-1], 0.1)
        except:
            raise ValueError('invalid column names')

        pltData = dict()
        for p in prots:
            if p in list(proteinData['Accession']):
                vals = proteinData.loc[proteinData.loc[:, 'Accession'] == p, colNames]
                pltData[p] = vals.values[0,:]
                
        self.axes.cla()
        for p, vec in pltData.items():
            paras = curve_fit(meltCurve, temps, vec, bounds=(0, [float('inf'), float('inf'), 0.3]))[0]
            self.axes.scatter(temps, vec, marker='.', label = p, s = 5)
            self.axes.plot(temps_, meltCurve(temps_, paras[0], paras[1], paras[2]), lw=1)
        
        self.axes.set_xlabel('Temperature (℃)', fontsize=5)
        self.axes.set_ylabel('Abundances', fontsize=5)
        self.axes.legend(fontsize=4)
        self.draw()
        
        
    def iTSA_Volcano(self, iTSA_result, fc_thres, pv_thres):        
        fc = iTSA_result['logFC']
        pv = iTSA_result['-logAdjPval']
        sig = np.where(np.logical_and(np.abs(fc) >= np.log2(fc_thres), pv >= -np.log10(pv_thres)))[0]
        
        self.axes.cla()
        self.axes.scatter(fc, pv, color = 'gray', marker = '.', s = 10)
        self.axes.scatter(fc[sig], pv[sig], color = 'red', marker = '.', s = 10)
        self.axes.axvline(x = np.log2(fc_thres),ls = '-', color = 'black', lw=1)
        self.axes.axvline(x = -np.log2(fc_thres),ls = '-', color = 'black', lw=1)
        self.axes.axhline(y = -np.log10(pv_thres), ls = '-', color = 'black', lw=1)
        self.axes.set_xlabel('Log FC', fontsize = 4)
        self.axes.set_ylabel('-Log Adj P', fontsize = 4)
        self.axes.tick_params(labelsize=4)
        self.draw()
    
    
    def PCAPlot(self, X, y):
        pca = PCA(n_components=2)
        X_s = StandardScaler().fit_transform(X.T)
        X_r = pca.fit(X_s).transform(X_s)
        label = np.unique(y)
        target_names = ['group_{}'.format(i) for i in label]
        
        self.axes.cla()
        for i in range(len(label)):
            self.axes.scatter(X_r[y == label[i], 0], X_r[y == label[i], 1], alpha=.8, lw=1, label=target_names[i], s=10)
        self.axes.set_xlabel('PC 1', fontsize = 4)
        self.axes.set_ylabel('PC 2', fontsize = 4)
        self.axes.tick_params(labelsize=4)
        self.draw()        
    
    
    def BarChart(self, X, y):
        self.axes.cla()
        data = [[], []]
        for i in range(X.shape[0]):
            for j in range(X.shape[1]):
                data[0].append(np.log2(X.iloc[i,j]))
                data[1].append(X.columns[j])
        data = pd.DataFrame(data).T
        data.columns = ['Intensity', 'Sample']
        sns.boxplot(ax=self.axes, x='Sample', y='Intensity', data=data)
    
    
    def HeatMap(self, X):
        self.axes.cla()
        sns.heatmap(ax=self.axes, data=np.log2(X), cmap="YlOrBr", yticklabels=False)
        
    