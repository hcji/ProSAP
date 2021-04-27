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

        



