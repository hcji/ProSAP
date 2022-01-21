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
from iTSA import iTSA, data_balance
from MakeFigure import MakeFigure
from ColumnSelectUI import ColumnSelectUI
from RunningUI import Running_Win


'''
proteinData = pd.read_csv('D:/project/CETSA_Benchmark/Data/Ouabain/iTSA.csv')
columns = ['sumionarea_protein_128H_Ctrl_1', 'sumionarea_protein_128H_Ctrl_2',
       'sumionarea_protein_128H_Ctrl_3', 'sumionarea_protein_128H_Ouabain_1',
       'sumionarea_protein_128H_Ouabain_2', 'sumionarea_protein_128H_Ouabain_3']
y = np.array([0,0,0,1,1,1])
X = proteinData.loc[:,columns]
names = proteinData.loc[:,'Accession'].values
'''


class AnaliTSAUI(QtWidgets.QWidget, Ui_Form):
    
    def __init__(self, parent=None): 
        super(AnaliTSAUI, self).__init__(parent)
        self.setupUi(self)
        self.setWindowTitle("iTSA Analysis")
        self.setWindowIcon(QtGui.QIcon("img/TPCA.ico"))

        self.figureVolcano = MakeFigure(20, 20, dpi = 300)
        self.figureVolcano_ntb = NavigationToolbar(self.figureVolcano, self)
        self.gridlayoutVolcano = QGridLayout(self.groupBoxVolcano)
        self.gridlayoutVolcano.addWidget(self.figureVolcano)
        self.gridlayoutVolcano.addWidget(self.figureVolcano_ntb)
        
        self.figureHeatmap = MakeFigure(10, 10, dpi = 200)
        self.figureHeatmap_ntb = NavigationToolbar(self.figureHeatmap, self)
        self.gridlayoutHeatmap = QGridLayout(self.groupBoxHeatmap)
        self.gridlayoutHeatmap.addWidget(self.figureHeatmap)
        self.gridlayoutHeatmap.addWidget(self.figureHeatmap_ntb)        
        
        self.figureBarchart = MakeFigure(10, 10, dpi = 200)
        self.figureBarchart_ntb = NavigationToolbar(self.figureBarchart, self)
        self.gridlayoutBarchart = QGridLayout(self.groupBoxBarchart)
        self.gridlayoutBarchart.addWidget(self.figureBarchart)
        self.gridlayoutBarchart.addWidget(self.figureBarchart_ntb)      
        
        self.figurePCA = MakeFigure(10, 10, dpi = 250)
        self.figurePCA_ntb = NavigationToolbar(self.figurePCA, self)
        self.gridlayoutPCA = QGridLayout(self.groupBoxPCA)
        self.gridlayoutPCA.addWidget(self.figurePCA)
        self.gridlayoutPCA.addWidget(self.figurePCA_ntb)  
        
        self.ColumnSelectUI = ColumnSelectUI()
        self.Running_Win = Running_Win()
        
        self.comboBoxMethod.addItems(['t-Test', 'Limma', 'edgeR', 'DESeq2'])
        self.comboBoxLog2.addItems(['True', 'False'])
        self.comboBoxBalance.addItems(['False', 'True'])
        
        self.pushButtonData.clicked.connect(self.LoadProteinFile)
        self.pushButtonOK.clicked.connect(self.DoPropress)
        self.pushButtonSave.clicked.connect(self.SaveData)
        
        self.data = None
        self.columns = None
        self.label = None
        self.RunningThread = None
    

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
            
            if 'Accession' in self.data.columns:
                pass
            else:
                self.ErrorMsg('Accession is not given in the data')
                return None
            
            self.tableViewData.setModel(TableModel(self.data))
            self.ColumnSelectUI.listWidget.clear()
            columns = self.data.columns
            for c in columns:
                self.ColumnSelectUI.listWidget.addItem(c)
            self.ColumnSelectUI.show()
            self.ColumnSelectUI.ButtonColumnSelect.clicked.connect(self.SetLabel)
            self.ColumnSelectUI.ButtonColumnCancel.clicked.connect(self.ColumnSelectUI.close)
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
            y = np.array([str(self.tableWidgetLabel.item(i,1).text()) for i in range(self.tableWidgetLabel.rowCount())])            
                        
            names = self.data.loc[:,'Accession'].values
            method = self.comboBoxMethod.currentText()
            if method in ['edgeR', 'DESeq2']:
                self.ErrorMsg('Block! The selected method is designed for transcriptomics with negative binomial distribution assumption. Whether the proteomics data are warranted for the mathematical assumption needs further validation. Contact the author if you want to use this function.')
                return
            worker = iTSA(method = method)
            
            nan_row = []
            for row in range(X.shape[0]):
                X.iloc[row,:] = pd.to_numeric(X.iloc[row,:], errors='coerce')
            X = X.astype(float)
            X[X <= 0] = np.nan
            
            for row in range(X.shape[0]):
                nan_row.append(X.iloc[row,:].isnull().sum())
            
            keep = np.where(np.array(nan_row) == 0)[0]
            X = X.iloc[keep,:]
            X = X.reset_index(drop = True)
            names = names[keep].astype(str)

            if self.comboBoxBalance.currentText() == 'True':
                X, y = data_balance(X, y)

            if self.comboBoxLog2.currentText() == 'True':
                if np.max(np.max(X)) >= 999:
                    self.ErrorMsg('The data seems not Log2 transformed, please change parameter')
                    return None
                else:
                    X = 2 ** X
            else:
                pass
            
            self.Running_Win.show()
            result = worker.fit_data(X, y, names)
            if result is None:
                self.ErrorMsg('There is something wrong with R, only t-test is support')
                self.Running_Win.close()
                return
            
            result = result.reset_index(drop=True)
            
            fc_thres = self.doubleSpinBoxFCthres.value()
            pv_thres = self.doubleSpinBoxPthres.value()

            self.tableViewData.setModel(TableModel(result))
            self.figureVolcano.iTSA_Volcano(result, fc_thres, pv_thres)
            self.figureHeatmap.CorrHeatMap(X)
            self.figureBarchart.BarChart(X, y)
            try:
                self.figurePCA.PCAPlot(X, y)
            except:
                pass
            self.Running_Win.close()
    
    
    def SaveData(self):
        options = QtWidgets.QFileDialog.Options()
        options |= QtWidgets.QFileDialog.DontUseNativeDialog
        fileName, _ = QtWidgets.QFileDialog.getSaveFileName(self, "Save", ".csv","All Files (*);;CSV Files (*.csv)", options=options)
        if fileName:
            data = pd.DataFrame(self.tableViewData.model()._data)
            data.to_csv(fileName, index = False)
        

if __name__ == '__main__':
    import sys
    
    app = QApplication(sys.argv)
    ui = AnaliTSAUI()
    ui.show()
    sys.exit(app.exec_())