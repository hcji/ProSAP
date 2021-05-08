# -*- coding: utf-8 -*-
"""
Created on Fri May  7 08:36:21 2021

@author: jihon
"""

import numpy as np
import pandas as pd

from PyQt5 import QtCore, QtGui, QtWidgets
from PyQt5.QtWidgets import QApplication, QMainWindow, QGridLayout, QHBoxLayout

from AnaliTSA import Ui_Form
from MakeFigure import MakeFigure
from matplotlib.backends.backend_qt5 import NavigationToolbar2QT as NavigationToolbar

from Utils import TableModel
from iTSA import iTSA
from Thread import PreprocessThread
from MakeFigure import MakeFigure
from ColumnSelectUI import ColumnSelectUI


'''
proteinData = pd.read_csv('data/iTSA_TableS1_Proteomics.csv')
columns = ['V_log2.i._TMT_1_iTSA52', 'V_log2.i._TMT_3_iTSA52',
       'V_log2.i._TMT_5_iTSA52', 'V_log2.i._TMT_7_iTSA52',
       'V_log2.i._TMT_9_iTSA52', 'D_log2.i._TMT_2_iTSA52',
       'D_log2.i._TMT_4_iTSA52', 'D_log2.i._TMT_6_iTSA52',
       'D_log2.i._TMT_8_iTSA52', 'D_log2.i._TMT_10_iTSA52']
labels = [0,0,0,0,0,1,1,1,1,1]
X = proteinData.loc[:,columns]
'''


class AnaliTSAUI(QtWidgets.QWidget, Ui_Form):
    
    def __init__(self, parent=None): 
        super(AnaliTSAUI, self).__init__(parent)
        self.setupUi(self)
        self.setWindowTitle("iTSA Analysis")
        self.setWindowIcon(QtGui.QIcon("img/TPCA.ico"))

        self.figureTSA = MakeFigure(10, 10, dpi = 250)
        self.figureTSA_ntb = NavigationToolbar(self.figureTSA, self)
        self.gridlayoutTSA = QGridLayout(self.groupBoxVolcano)
        self.gridlayoutTSA.addWidget(self.figureTSA)
        self.gridlayoutTSA.addWidget(self.figureTSA_ntb)
        self.ColumnSelectUI = ColumnSelectUI()
        self.comboBoxMethod.addItems(['t-Test', 'Limma'])
        self.comboBoxLog2.addItems(['True', 'False'])
        
        self.pushButtonData.clicked.connect(self.LoadProteinFile)
        self.pushButtonOK.clicked.connect(self.DoPropress)
        
        self.data = None
        self.columns = None
        self.label = None
    
    
    def LoadProteinFile(self):
        options = QtWidgets.QFileDialog.Options()
        options |= QtWidgets.QFileDialog.DontUseNativeDialog
        fileName, _ = QtWidgets.QFileDialog.getOpenFileName(self,"QFileDialog.getOpenFileName()", "","All Files (*);;CSV Files (*.csv)", options=options)
        if fileName:
            if fileName.split('.')[1] == 'csv':
                self.data = pd.read_csv(fileName)
                self.tableViewData.setModel(TableModel(self.data))
                
                self.columns = self.data.columns
                for c in self.columns:
                    self.ColumnSelectUI.listWidget.addItem(c)
                self.ColumnSelectUI.show()
                self.ColumnSelectUI.ButtonColumnSelect.clicked.connect(self.SetLabel)
                self.ColumnSelectUI.ButtonColumnCancel.clicked.connect(self.ColumnSelectUI.close)
            else:
                self.ErrorMsg("Invalid format")
        else:
            pass
    
    
    def SetLabel(self):
        columns = [i.text() for i in self.ColumnSelectUI.listWidget.selectedItems()]
        self.columns = columns
        self.ColumnSelectUI.close()
        
        self.tableWidgetLabel.setRowCount(len(columns))
        self.tableWidgetLabel.setColumnCount(2)
        self.tableWidgetLabel.setHorizontalHeaderLabels(['columnName', 'label'])
        self.tableWidgetLabel.setVerticalHeaderLabels(np.arange(len(columns)).astype(str))
        for i in range(len(columns)):
            item = QtWidgets.QTableWidgetItem(str(columns[i]))
            self.tableWidgetLabel.setItem(i, 0, item)
    
    
    def DoPropress(self):
        if None in [self.tableWidgetLabel.item(i,1) for i in range(self.tableWidgetLabel.rowCount())]:
            msg = QtWidgets.QMessageBox()
            msg.resize(550, 200)
            msg.setIcon(QtWidgets.QMessageBox.Critical)
            msg.setText("Please input labels")
            msg.setWindowTitle("Error")
            msg.exec_()
        elif len(np.unique([str(self.tableWidgetLabel.item(i,1).text()) for i in range(self.tableWidgetLabel.rowCount())])) != 2:
            msg = QtWidgets.QMessageBox()
            msg.resize(550, 200)
            msg.setIcon(QtWidgets.QMessageBox.Critical)
            msg.setText("Only two unique labels allowed")
            msg.setWindowTitle("Error")
            msg.exec_()               
        else:
            X = self.data.loc[:,self.columns]
            if self.comboBoxLog2.currentText() == 'True':
                pass
            else:
                X = np.log2(X)
            y = np.array([str(self.tableWidgetLabel.item(i,1).text()) for i in range(self.tableWidgetLabel.rowCount())])
            names = self.data.loc[:,'Accession']
            method = self.comboBoxMethod.currentText()
            worker = iTSA(method = method)
            result = worker.fit_data(X, y, names)
            result = result.reset_index(drop=True)
            
            fc_thres = self.doubleSpinBoxFCthres.value()
            pv_thres = self.doubleSpinBoxPthres.value()
            
            self.tableViewData.setModel(TableModel(result))
            self.figureTSA.iTSA_Volcano(result, fc_thres, pv_thres)
            
        
        

if __name__ == '__main__':
    import sys
    
    app = QApplication(sys.argv)
    ui = AnaliTSAUI()
    ui.show()
    sys.exit(app.exec_())