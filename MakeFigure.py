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



class MakeFigure(FigureCanvas):
    def __init__(self,width=5, height=4, dpi=300):
        self.fig = Figure(figsize=(width, height), dpi=dpi)
        super(MakeFigure,self).__init__(self.fig) 
        self.axes = self.fig.add_subplot(111)
    
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
        self.axes.set_xlabel('temperature (℃)')
        self.axes.set_ylabel('abundances')
        self.axes.legend()
        