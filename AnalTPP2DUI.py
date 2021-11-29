# -*- coding: utf-8 -*-
"""
Created on Mon Nov 29 10:29:37 2021

@author: jihon
"""

import numpy as np
import pandas as pd

from PyQt5 import QtCore, QtGui, QtWidgets
from PyQt5.QtWidgets import QApplication, QMainWindow, QGridLayout, QHBoxLayout

from AnalTPP2D import Ui_Form
from MakeFigure import MakeFigure
from matplotlib.backends.backend_qt5 import NavigationToolbar2QT as NavigationToolbar

from Utils import TableModel
from TPP2D import TPP2D

from RunningUI import Running_Win


'''
data = pd.read_csv('data/TPP2D/panobinostat_tpp2d_cell.csv')

'''


class AnalTPP2D(QtWidgets.QWidget, Ui_Form):
    
    def __init__(self, parent=None): 
        super(AnalTPP2D, self).__init__(parent)
        self.setupUi(self)
        self.setWindowTitle("TPP2D Analysis")
        self.setWindowIcon(QtGui.QIcon("img/TPCA.ico"))
        self.Running_Win = Running_Win()
        
        self.comboBoxMethod.addItems(['Kurzawa'])
        
        self.pushButtonData.clicked.connect(self.LoadProteinFile)
        self.pushButtonOK.clicked.connect(self.DoPropress)
        self.pushButtonClose.clicked.connect(self.close)
        
        self.data = None
        self.maxit = None
        self.B = None
        self.result = None
    
    
    def WarnMsg(self, Text):
        msg = QtWidgets.QMessageBox()
        msg.resize(550, 200)
        msg.setIcon(QtWidgets.QMessageBox.Warning)
        msg.setText(Text)
        msg.setWindowTitle("Warning")
        msg.exec_()    
    
    
    def ErrorMsg(self, Text):
        msg = QtWidgets.QMessageBox()
        msg.resize(550, 200)
        msg.setIcon(QtWidgets.QMessageBox.Critical)
        msg.setText(Text)
        msg.setWindowTitle("Error")
        msg.exec_()

    
    def LoadProteinFile(self):
        options = QtWidgets.QFileDialog.Options()
        options |= QtWidgets.QFileDialog.DontUseNativeDialog
        fileName, _ = QtWidgets.QFileDialog.getOpenFileName(self, "Load", "","All Files (*);;CSV Files (*.csv)", options=options)
        if fileName:
            if fileName.split('.')[-1] == 'csv':
                self.data = pd.read_csv(fileName)
            elif fileName.split('.')[-1] == 'xlsx':
                self.data = pd.read_excel(fileName)
            else:
                self.ErrorMsg("Invalid format")
                return None
        self.tableViewData.setModel(TableModel(self.data))
        
        
    def DoPropress(self):
        self.Running_Win.show()
        tpp2d = TPP2D(self.data)
        tpp2d.check()
        result = tpp2d.fit_data(self.maxit, self.B)
        if result is None:
            self.ErrorMsg('There is something wrong with R, please check')
            self.Running_Win.close()
            return
        elif result == np.nan:
            self.ErrorMsg('There is something wrong when processing your data, please check')
            self.Running_Win.close()            
        else:
            self.result = result
            self.Running_Win.close()
    
        

if __name__ == '__main__':
    import sys
    
    app = QApplication(sys.argv)
    ui = AnalTPP2D()
    ui.show()
    sys.exit(app.exec_())