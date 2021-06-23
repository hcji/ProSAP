# -*- coding: utf-8 -*-
"""
Created on Wed Jun 23 14:00:41 2021

@author: hcji
"""


import numpy as np
import pandas as pd

from PyQt5 import QtCore, QtGui, QtWidgets
from PyQt5.QtWidgets import QApplication, QMainWindow, QGridLayout
from matplotlib.backends.backend_qt5 import NavigationToolbar2QT as NavigationToolbar

from AnalQC import Ui_MainWindow
from Utils import TableModel
from MakeFigure import MakeFigure


def calcRSD(data, Missing=False):
    vals = data.iloc[:,1:]
    for c in vals.columns:
        vals[c] = pd.to_numeric(vals[c], errors='coerce')
    
    missingNum = np.sum(vals.isnull().values, axis = 1)
    outputRSD = np.repeat(np.nan, len(missingNum))
    
    if Missing:
        keep = np.arange(len(missingNum))
    else:
        keep = np.where(missingNum == 0)[0]
        vals = vals.iloc[keep,:]
    
    medians = np.nanmedian(vals, axis = 0)
    ratios = medians / np.mean(medians)
    vals = vals / ratios
    
    sd = np.std(vals, axis = 1)
    m = np.mean(vals, axis = 1)
    rsd = sd / m
    outputRSD[keep] = rsd
    
    return missingNum, outputRSD


class AnalQCUI(QtWidgets.QMainWindow, Ui_MainWindow):
    
    def __init__(self, parent=None): 
        super(AnalQCUI, self).__init__(parent)
        self.setupUi(self)
        self.setWindowTitle("QC Analysis")
        self.setWindowIcon(QtGui.QIcon("img/TPCA.ico"))
          
        self.allData = []
        self.allRSD = []
        
        self.figureRSD = MakeFigure(10, 10, dpi = 250)
        self.figureRSD_ntb = NavigationToolbar(self.figureRSD, self)
        self.gridlayoutRSD = QGridLayout(self.groupBox1)
        self.gridlayoutRSD.addWidget(self.figureRSD)
        self.gridlayoutRSD.addWidget(self.figureRSD_ntb)
        
        self.figureBox = MakeFigure(10, 10, dpi = 250)
        self.figureBox_ntb = NavigationToolbar(self.figureBox, self)
        self.gridlayoutBox = QGridLayout(self.groupBox2)
        self.gridlayoutBox.addWidget(self.figureBox)
        self.gridlayoutBox.addWidget(self.figureBox_ntb)
        
        self.ButtonOpen.clicked.connect(self.LoadProteinFile)
        
        
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
        fileNames, _ = QtWidgets.QFileDialog.getOpenFileNames(self,"Load", "","All Files (*);;CSV Files (*.csv)", options=options)
        
        if len(fileNames) == 0:
            pass
        
        else:
            self.tableWidgetFile.setRowCount(len(fileNames))
            self.tableWidgetFile.setColumnCount(2)
            self.tableWidgetFile.setHorizontalHeaderLabels(['FilePath', 'Label'])
            self.tableWidgetFile.setVerticalHeaderLabels(np.arange(len(fileNames)).astype(str))
            for i in range(len(fileNames)):
                item = QtWidgets.QTableWidgetItem(str(fileNames[i]))
                self.tableWidgetFile.setItem(i, 0, item)


    
    def LoadData(self):
        self.allData = []
        fileNames = [self.tableWidgetFile.item(i,0) for i in range(self.tableWidgetFile.rowCount())]
        for f in fileNames:
            if f.split('.')[1] == 'csv':
                data = pd.read_csv(f)             
            elif f.text().split('.')[1] == 'xlsx':
                data = pd.read_excel(f)
            else:
                self.ErrorMsg('Cannot be load the selected file')
                return None
            
            if 'Accession' in data.columns:
                self.ErrorMsg('Accession is not given in the data')
                return None
            else:
                self.allData.append(data)
    
    
    def ShowRSD(self):
        allData = self.allData
        if len(allData) == 0:
            return None
        
        fileNames = [self.tableWidgetFile.item(i,0) for i in range(self.tableWidgetFile.rowCount())]
        fileLabel = [self.tableWidgetFile.item(i,1) for i in range(self.tableWidgetFile.rowCount())]
        if None in fileLabel:
            fileLabel = ['file_{}'.format(i) for i in range(len(fileNames))]
            
        self.allRSD = []
        for i in range(len(fileLabel)):
            nMissing, rsd = calcRSD(allData[i])
    
        
        
if __name__ == '__main__':
    import sys
    
    app = QApplication(sys.argv)
    ui = AnalQCUI()
    ui.show()
    sys.exit(app.exec_())
    