# -*- coding: utf-8 -*-
"""
Created on Fri Apr 16 16:23:29 2021

@author: hcji
"""


import matplotlib
import numpy as np
import matplotlib.pyplot as plt
matplotlib.use("Qt5Agg")
from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg as FigureCanvas
from matplotlib.figure import Figure
from PyQt5 import QtCore, QtGui, QtWidgets


class MakeFigure(FigureCanvas):
    def __init__(self,width=5, height=4, dpi=300):
        self.fig = Figure(figsize=(width, height), dpi=dpi)
        self.fig.subplots_adjust(top=0.95,bottom=0.2,left=0.18,right=0.95)
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

        for p, vec in pltData.items():
            self.axes.plot(temps, vec, label = p)
            
        self.axes.set_xlabel('Temperature (℃)', fontsize=5)
        self.axes.set_ylabel('Abundances', fontsize=5)
        self.axes.legend(fontsize=4)
    
    
    def AverageTSAFigure(self, proteinData1, proteinData2, colNames):
        try:
            temps = np.array([float(t.replace('T', '')) for t in colNames])
        except:
            raise ValueError('invalid column names')
        
        vec_1 = np.mean(proteinData1.loc[:, colNames], axis = 0)
        vec_2 = np.mean(proteinData2.loc[:, colNames], axis = 0)
        
        self.axes.plot(temps, vec_1, label = 'Group 1')
        self.axes.plot(temps, vec_2, label = 'Group 2')
        
        self.axes.tick_params(labelsize=4)
        self.axes.set_xlabel('Temperature (℃)', fontsize=3)
        self.axes.set_ylabel('Abundances', fontsize=3)
        self.axes.legend(fontsize=3, bbox_to_anchor=(1,1), loc="upper left")        
        
    
    def SingleTSAFigure(self, proteinData1, proteinData2, colNames, ProteinAccession):
        try:
            temps = np.array([float(t.replace('T', '')) for t in colNames])
        except:
            raise ValueError('invalid column names')
            
        vec_1 = proteinData1.loc[proteinData1.loc[:, 'Accession'] == ProteinAccession, colNames].values[0,:]
        vec_2 = proteinData2.loc[proteinData2.loc[:, 'Accession'] == ProteinAccession, colNames].values[0,:]
        
        self.axes.plot(temps, vec_1, label = 'Group 1_{}'.format(ProteinAccession))
        self.axes.plot(temps, vec_2, label = 'Group 2_{}'.format(ProteinAccession))
        
        self.axes.tick_params(labelsize=4)
        self.axes.set_xlabel('Temperature (℃)', fontsize=3)
        self.axes.set_ylabel('Abundances', fontsize=3)
        self.axes.legend(fontsize=3)
    
    
    def RSDHistFigure(self, rsdList):
        self.axes.hist(rsdList, 100)
        self.axes.tick_params(labelsize=3)
        self.axes.set_xlabel('RSD', fontsize=3)
        self.axes.set_ylabel('Number', fontsize=3)
    
    
    def ROCFigure(self, fpr, tpr, auroc):
        self.axes.plot(fpr, tpr, label='AUC = {}'.format(auroc), color = 'red', lw=0.7)
        self.axes.plot([0, 1], [0, 1], color='black', linestyle='--', lw=0.7)
        self.axes.set_xlabel('False Positive Rate', fontsize = 4)
        self.axes.set_ylabel('True Positive Rate', fontsize = 4)
        self.axes.tick_params(labelsize=4)
        self.axes.legend(fontsize=3)
        