# -*- coding: utf-8 -*-
"""
Created on Thu Apr 15 15:57:14 2021

@author: hcji
"""

import os
import numpy as np
import pandas as pd

from PyQt5 import QtCore, QtGui, QtWidgets
from PyQt5.QtWidgets import QApplication, QMainWindow

from Preprocess import Ui_Form
from ColumnSelectUI import ColumnSelectUI
from Utils import TableModel
from MakeFigure import MakeFigure

class PreprocessUI(QtWidgets.QWidget, Ui_Form):
    
    def __init__(self, parent=None): 
        super(PreprocessUI, self).__init__(parent)
        self.setupUi(self)
        self.setWindowTitle("Preprocessing")
        
        self.comboBoxMerging.addItems(['Median', 'Mean'])
        
        self.ColumnSelectUI = ColumnSelectUI()
        
        self.pushButtonOpen.clicked.connect(self.LoadProteinFile)
        self.pushButtonConfirm.clicked.connect(self.DoPropress)
        self.pushButtonClear.clicked.connect(self.ClearProteinFile)
        self.pushButtonSave.clicked.connect(self.SaveData)
        
        self.columns = None
    
    
    def ClearProteinFile(self):
        self.ListFile.clear()
    
        
    def LoadProteinFile(self):
        self.ListFile.clear()
        self.comboBoxReference.clear()
        
        options = QtWidgets.QFileDialog.Options()
        options |= QtWidgets.QFileDialog.DontUseNativeDialog
        fileNames, _ = QtWidgets.QFileDialog.getOpenFileNames(self,"QFileDialog.getOpenFileName()", "","All Files (*);;CSV Files (*.csv)", options=options)
        
        if len(fileNames) == 0:
            pass
        
        else:
            for fileName in fileNames:
                if fileName.split('.')[1] == 'csv':
                    self.ListFile.addItem(fileName)
                else:
                    pass
                
            if self.ListFile.count() == 0:
                pass
            
            else:    
                data = pd.read_csv(fileName)
        
                all_cols = data.columns
                self.comboBoxPSM.addItems(all_cols)
                for c in all_cols:
                    self.ColumnSelectUI.listWidget.addItem(c)
        
                self.ColumnSelectUI.show()
                self.ColumnSelectUI.ButtonColumnSelect.clicked.connect(self.SetTemperture)
                self.ColumnSelectUI.ButtonColumnCancel.clicked.connect(self.ColumnSelectUI.close)
        
        
    def SetTemperture(self):
        columns = [i.text() for i in self.ColumnSelectUI.listWidget.selectedItems()]
        self.columns = columns
        self.tableWidgetTemp.setRowCount(len(columns))
        self.tableWidgetTemp.setColumnCount(2)
        self.tableWidgetTemp.setHorizontalHeaderLabels(['columnName', 'temperture'])
        self.tableWidgetTemp.setVerticalHeaderLabels(np.arange(len(columns)).astype(str))
        for i in range(len(columns)):
            item = QtWidgets.QTableWidgetItem(str(columns[i]))
            self.tableWidgetTemp.setItem(i, 0, item)
        self.comboBoxReference.addItems(columns)
        self.ColumnSelectUI.close
    
    
    def DoPropress(self):
        
        # d = 'D:/project/PDFiles'
        # fileNames = os.listdir(d)
        fileNames = [self.ListFile.item(i).text() for i in range(self.ListFile.count())]
        columns = [i.text() for i in self.ColumnSelectUI.listWidget.selectedItems()]
        reference = self.comboBoxReference.currentText()
        
        if None in [self.tableWidgetTemp.item(i,1) for i in range(self.tableWidgetTemp.rowCount())]:
            msg = QtWidgets.QMessageBox()
            msg.resize(550, 200)
            msg.setIcon(QtWidgets.QMessageBox.Critical)
            msg.setText("Please input temperature")
            msg.setWindowTitle("Error")
            msg.exec_()
        else:
            
            tempertures = [float(self.tableWidgetTemp.item(i,1).text()) for i in range(self.tableWidgetTemp.rowCount())]
        
            psm_column = self.comboBoxPSM.currentText()
            psm_thres = self.spinBoxPSMFilter.value()
            std_thres = self.doubleSpinBoxRSDFilter.value()
            
            all_data = [pd.read_csv(f) for f in fileNames]
            for i in range(len(all_data)):
                ref = all_data[i].loc[:,reference].copy()
                for c in columns:
                    all_data[i].loc[:,c] /= ref
        
            data = all_data[0]
            if len(all_data) == 1:
                pass
            else:
                for data_ in all_data[1:]:
                    data = data.merge(data_, 'outer', on = 'Accession')
        
            if self.comboBoxMerging.currentText() == 'Median':
                fun = np.nanmedian
            else:
                fun = np.nanmean
        
            prot_list, val_list, rsd_list = [], [], []
            for i in data.index:
                wh = np.where([psm_column == s.split('_')[0] for s in data.columns])[0]
                psm = np.nanmean(data.iloc[i, wh].values.astype(float))
                if psm < psm_thres:
                    continue
                else:
                    prot = data.loc[i, 'Accession']
            
                vals = []
                for c in columns:
                    wh = np.where([c == s.split('_')[0] for s in data.columns])[0]
                    v = data.iloc[i, wh].values.astype(float)
                    v = np.round(v, 4)
                    std = np.nanstd(v) / np.nanmean(v)
                    vals.append(fun(v))
                if std > std_thres:
                    continue
            
                vals = np.array(vals)
                vals[np.isnan(vals)] = 0
            
                rsd_list.append(std)
                prot_list.append(prot)
                val_list.append(vals)
        
            val_cols = ['T{}'.format(t) for t in tempertures]
            val_list = pd.DataFrame(val_list)
            val_list.columns = val_cols
            val_list['Accession'] = prot_list
            val_cols = ['Accession'] + val_cols

            F = MakeFigure(1.2, 1.2)
            F.axes.cla()
            F.RSDHistFigure(np.array(rsd_list))
            f = QtWidgets.QGraphicsScene()
            f.addWidget(F)
            self.graphicsView.setScene(f)
        
            result = val_list[val_cols]
            self.tableView.setModel(TableModel(result))
        
        
    def SaveData(self):
        options = QtWidgets.QFileDialog.Options()
        options |= QtWidgets.QFileDialog.DontUseNativeDialog
        fileName, _ = QtWidgets.QFileDialog.getSaveFileName(self,"QFileDialog.getOpenFileName()", "","All Files (*);;CSV Files (*.csv)", options=options)
        data = pd.DataFrame(self.tableView.model()._data, index = False)
        data.to_csv(fileName)
        
        

if __name__ == '__main__':
    import sys
    
    app = QApplication(sys.argv)
    ui = PreprocessUI()
    ui.show()
    sys.exit(app.exec_())