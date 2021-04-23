# -*- coding: utf-8 -*-
"""
Created on Wed Apr 14 09:29:07 2021

@author: hcji
"""


import numpy as np
import pandas as pd

from PyQt5.QtCore import Qt
from PyQt5 import QtCore, QtGui, QtWidgets
from PyQt5.QtWidgets import QApplication, QMainWindow

from MainWindow import Ui_MainWindow
from ColumnSelectUI import ColumnSelectUI
from AnalPvalComplexUI import AnalPvalComplexUI
from AnalROCUI import AnalROCUI
from AnalTSAUI import AnalTSAUI
from MakeFigure import MakeFigure
from Utils import TableModel

class TCPA_Main(QMainWindow, Ui_MainWindow):
    
    def __init__(self, parent=None):
        super(TCPA_Main, self).__init__(parent)
        self.setupUi(self)
        
        # main window
        self.resize(1300, 800)
        self.setMinimumWidth(1150)
        self.setMinimumHeight(650)
        self.move(75, 50)
        self.setWindowTitle("TPCA -- Thermal proximity coaggregation analysis")
        
        # widgets
        self.ColumnSelectUI = ColumnSelectUI()
        self.AnalPvalComplexUI = AnalPvalComplexUI()
        self.AnalROCUI = AnalROCUI()
        self.AnalTSAUI = AnalTSAUI()
        
        # menu action
        self.actionProteomics.triggered.connect(self.LoadProteinFile)
        self.actionDatabase.triggered.connect(self.LoadProteinComplex)
        self.action_CETSA.triggered.connect(self.OpenAnalTSA)
        self.actionCalcROC.triggered.connect(self.OpenAnalROC)
        self.actionCalcPval.triggered.connect(self.OpenAnalPvalComplex)
        
        # button action
        self.ButtonGroup1.clicked.connect(self.SetProteinTable1)
        self.ButtonGroup2.clicked.connect(self.SetProteinTable2)
        self.ButtonClearFileList.clicked.connect(self.ClearProteinFile)
        self.ButtonDatabaseConfirm.clicked.connect(self.SetProteinComplex)
        self.ButtonDatabaseRemove.clicked.connect(self.RemoveProteinComplex)
        self.ButtonClearDatabase.clicked.connect(self.ClearProteinComplex)
        
        self.ButtonShowCurve.clicked.connect(self.ProteinComplexCurve)
        
        # server data
        self.UsedColumns = []
        
    
    def LoadProteinFile(self):
        options = QtWidgets.QFileDialog.Options()
        options |= QtWidgets.QFileDialog.DontUseNativeDialog
        fileName, _ = QtWidgets.QFileDialog.getOpenFileName(self,"QFileDialog.getOpenFileName()", "","All Files (*);;CSV Files (*.csv)", options=options)
        if fileName:
            if fileName.split('.')[1] == 'csv':
                self.ListFile.addItem(fileName)
            else:
                msg = QtWidgets.QMessageBox()
                msg.resize(550, 200)
                msg.setIcon(QtWidgets.QMessageBox.Critical)
                msg.setText("Invalid format")
                msg.setWindowTitle("Error")
                msg.exec_()
        else:
            pass


    def ClearProteinFile(self):
        self.ListFile.clear()
        
    
    def LoadProteinComplex(self):
        options = QtWidgets.QFileDialog.Options()
        options |= QtWidgets.QFileDialog.DontUseNativeDialog
        fileName, _ = QtWidgets.QFileDialog.getOpenFileName(self,"QFileDialog.getOpenFileName()", "","All Files (*);;CSV Files (*.csv)", options=options)
        if fileName:
            if fileName.split('.')[1] == 'csv':
                self.ListDatabase.addItem(fileName)
            else:
                msg = QtWidgets.QMessageBox()
                msg.resize(550, 200)
                msg.setIcon(QtWidgets.QMessageBox.Critical)
                msg.setText("Invalid format")
                msg.setWindowTitle("Error")
                msg.exec_()
        else:
            pass
        
    
    def RemoveProteinComplex(self):
        self.ListDatabase.takeItem(self.ListDatabase.currentItem())
    
    
    def ClearProteinComplex(self):
        self.ListDatabase.clear()
    
    
    def SelectProteinTable(self):
        selectItem = self.ListFile.currentItem()
        selectData = pd.read_csv(selectItem.text())
        return selectData


    def SetProteinTable1(self):
        data = self.SelectProteinTable()
        self.tableProtein1.setModel(TableModel(data))


    def SetProteinTable2(self):
        data = self.SelectProteinTable()
        self.tableProtein2.setModel(TableModel(data))


    def SetProteinComplex(self):
        selectItem = self.ListDatabase.currentItem()
        selectData = pd.read_csv(selectItem.text())
        if 'Subunits_UniProt_IDs' not in selectData.columns:
            msg = QtWidgets.QMessageBox()
            msg.resize(550, 200)
            msg.setIcon(QtWidgets.QMessageBox.Critical)
            msg.setText("'Subunits_UniProt_IDs' not in columns")
            msg.setWindowTitle("Error")
            msg.exec_()
        else:
            self.tableProteinComplex.setRowCount(selectData.shape[0])
            self.tableProteinComplex.setColumnCount(selectData.shape[1])
            self.tableProteinComplex.setHorizontalHeaderLabels(selectData.columns)
            self.tableProteinComplex.setVerticalHeaderLabels(selectData.index.astype(str))
            for i in range(selectData.shape[0]):
                for j in range(selectData.shape[1]):
                    item = QtWidgets.QTableWidgetItem(str(selectData.iloc[i,j]))
                    self.tableProteinComplex.setItem(i, j, item)
    
    
    def ProteinComplexCurve(self):
        all_cols = self.tableProtein1.model()._data.columns
        for c in all_cols:
            self.ColumnSelectUI.listWidget.addItem(c)
        self.ColumnSelectUI.show()
        self.ColumnSelectUI.ButtonColumnSelect.clicked.connect(self.PlotProteinComplex)
        self.ColumnSelectUI.ButtonColumnCancel.clicked.connect(self.ColumnSelectUI.close)
        
    
    def OpenAnalROC(self):
        self.AnalROCUI.show()


    def OpenAnalPvalComplex(self):
        self.AnalPvalComplexUI.show()
        
    
    def OpenAnalTSA(self):
        self.AnalTSAUI.show()
        if self.tableProtein1.model() is None or (self.tableProtein2.model() is None):
            pass
        else:
            proteinData1 = self.tableProtein1.model()._data
            proteinData2 = self.tableProtein2.model()._data
        
        
        
        
        
    
    
    def PlotProteinComplex(self):
        colNames = [i.text() for i in self.ColumnSelectUI.listWidget.selectedItems()] 
        header = [self.tableProteinComplex.horizontalHeaderItem(i).text() for i in range(self.tableProteinComplex.columnCount())]
        # print(header)
        i = self.tableProteinComplex.selectedItems()[0].row()
        j = header.index('Subunits_UniProt_IDs')
        proteinSubunit = self.tableProteinComplex.item(i, j).text()
        # print(proteinSubunit)
        proteinData1 = self.tableProtein1.model()._data
        proteinData2 = self.tableProtein2.model()._data
        # print(proteinData)
        F_gr1 = MakeFigure(2.2, 1.6)
        F_gr2 = MakeFigure(2.2, 1.6)
        F_gr1.axes.cla()
        F_gr2.axes.cla()
        
        F_gr1.ProteinComplexFigure(proteinSubunit, proteinData1, colNames)
        F_gr2.ProteinComplexFigure(proteinSubunit, proteinData2, colNames)
        F1 = QtWidgets.QGraphicsScene()
        F1.addWidget(F_gr1)
        F2 = QtWidgets.QGraphicsScene()
        F2.addWidget(F_gr2)
        
        self.GraphicThermShift1.setScene(F1)
        self.GraphicThermShift2.setScene(F2)
        
        self.ColumnSelectUI.close
        
    

if __name__ == '__main__':
    import sys
    
    app = QApplication(sys.argv)
    ui = TCPA_Main()
    ui.show()
    sys.exit(app.exec_())