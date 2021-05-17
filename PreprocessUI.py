# -*- coding: utf-8 -*-
"""
Created on Thu Apr 15 15:57:14 2021

@author: hcji
"""

import os
import numpy as np
import pandas as pd

from sklearn.impute import KNNImputer
from PyQt5 import QtCore, QtGui, QtWidgets
from PyQt5.QtWidgets import QApplication, QMainWindow, QGridLayout
from matplotlib.backends.backend_qt5 import NavigationToolbar2QT as NavigationToolbar

from Preprocess import Ui_Form
from ColumnSelectUI import ColumnSelectUI
from Utils import TableModel
from Thread import PreprocessThread
from MakeFigure import MakeFigure


class PreprocessUI(QtWidgets.QWidget, Ui_Form):
    
    def __init__(self, parent=None): 
        super(PreprocessUI, self).__init__(parent)
        self.setupUi(self)
        self.setWindowTitle("Preprocessing")
        self.setWindowIcon(QtGui.QIcon("img/TPCA.ico"))
        
        self.comboBoxMerging.addItems(['Median', 'Mean'])
        self.comboBoxNorm.addItems(['Reference', 'Median', 'None'])
        self.comboBoxMV.addItems(['None', 'KNN', 'Zero'])

        self.ColumnSelectUI = ColumnSelectUI()
        self.PreprocessThread = None
        
        self.figureRSD = MakeFigure(5, 5)
        self.figureRSD_ntb = NavigationToolbar(self.figureRSD, self)
        self.gridlayoutRSD = QGridLayout(self.groupBox)
        self.gridlayoutRSD.addWidget(self.figureRSD)
        self.gridlayoutRSD.addWidget(self.figureRSD_ntb)
        
        self.pushButtonOpen.clicked.connect(self.LoadProteinFile)
        self.pushButtonConfirm.clicked.connect(self.DoPropress)
        self.pushButtonClear.clicked.connect(self.ClearProteinFile)
        self.pushButtonSave.clicked.connect(self.SaveData)
        
        self.columns = None
        self.tempertures = None
        self.valueData = []
        self.rsdData = []
        self.protData = []
    
    
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
                if fileName.split('.')[1] in ['csv', 'xls', 'xlsx']:
                    self.ListFile.addItem(fileName)
                else:
                    pass
                
            if self.ListFile.count() == 0:
                pass
            
            else:
                if fileName.split('.')[1] == 'csv':
                    data = pd.read_csv(fileName)
                else:
                    data = pd.read_excel(fileName)
        
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
        self.ColumnSelectUI.close()
        
        self.tableWidgetTemp.setRowCount(len(columns))
        self.tableWidgetTemp.setColumnCount(2)
        self.tableWidgetTemp.setHorizontalHeaderLabels(['columnName', 'temperture'])
        self.tableWidgetTemp.setVerticalHeaderLabels(np.arange(len(columns)).astype(str))
        for i in range(len(columns)):
            item = QtWidgets.QTableWidgetItem(str(columns[i]))
            self.tableWidgetTemp.setItem(i, 0, item)
        self.comboBoxReference.addItems(columns)
    
    
    def DoPropress(self):
        # d = 'D:/project/PDFiles'
        # fileNames = os.listdir(d)
        fileNames = [self.ListFile.item(i).text() for i in range(self.ListFile.count())]
        fileNames_ = []
        for f in fileNames:
            if f.split('.')[1] in ['csv', 'xlsx', 'xls']:
                fileNames_.append(f)
        fileNames = fileNames_
        
        columns = [i.text() for i in self.ColumnSelectUI.listWidget.selectedItems()]
        reference = self.comboBoxReference.currentText()
        
        if None in [self.tableWidgetTemp.item(i,1) for i in range(self.tableWidgetTemp.rowCount())]:
            msg = QtWidgets.QMessageBox()
            msg.resize(550, 200)
            msg.setIcon(QtWidgets.QMessageBox.Information)
            msg.setText("No temperature input, use original colnames")
            msg.setWindowTitle("Information")
            msg.exec_()
            self.tempertures = None
        else:
            try:
                tempertures = [float(self.tableWidgetTemp.item(i,1).text()) for i in range(self.tableWidgetTemp.rowCount())]
            except:
                tempertures = None    
            self.tempertures = tempertures    
        psm_column = self.comboBoxPSM.currentText()
        psm_thres = self.spinBoxPSMFilter.value()
        std_thres = self.doubleSpinBoxRSDFilter.value()
            
        all_data = []
        for f in fileNames:
            if f.split('.')[1]  == 'csv':
                all_data.append(pd.read_csv(f))
            else:
                all_data.append(pd.read_excel(f))
            
        # devided by reference
        if self.comboBoxNorm.currentText() == 'Reference':
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
            
        self.valueData, self.rsdData, self.protData = [], [], []
            
        # normalize by median
        if self.comboBoxNorm.currentText() == 'Median':
            whs = np.where([c.split('_')[0] in columns for c in list(data.columns)])[0]
            medians = np.nanmedian(data.iloc[:,whs], axis = 0)
            ratios = medians / np.mean(medians)
            data.iloc[:,whs] = data.iloc[:,whs] / ratios
            
        self.PreprocessThread = PreprocessThread(data, psm_column, psm_thres, std_thres, columns, fun)
        self.PreprocessThread._ind.connect(self.ProcessBar)
        self.PreprocessThread._val.connect(self.ValueData)
        self.PreprocessThread._rsd.connect(self.RSDData)
        self.PreprocessThread._prot.connect(self.ProtData)
        self.PreprocessThread.start()
        self.PreprocessThread.finished.connect(self.VisualizeProprocess)            
    
    
    def VisualizeProprocess(self):
        val_list = pd.DataFrame(self.valueData)
        if self.comboBoxMV.currentText() == 'KNN':
            knn_imputer = KNNImputer(n_neighbors = 3)
            val_list_ = pd.DataFrame(knn_imputer.fit_transform(val_list))
            val_list_.columns = val_list.columns
            val_list = val_list_
        elif self.comboBoxMV.currentText() == 'Zero':
            val_list = val_list.fillna(0)
        else:
            pass
        
        if self.tempertures is not None:
            val_cols = ['T{}'.format(t) for t in self.tempertures]
        else:
            val_cols = [self.tableWidgetTemp.item(i,0).text() for i in range(self.tableWidgetTemp.rowCount())]
        val_list.columns = val_cols
        val_list['Accession'] = self.protData
        val_list['RSD'] = np.round(self.rsdData, 4)
        val_cols = ['Accession'] + val_cols + ['RSD']

        self.figureRSD.RSDHistFigure(np.array(self.rsdData))
        result = val_list[val_cols]
        result = result[result['RSD'] <= self.doubleSpinBoxRSDFilter.value()]
        
        self.tableView.setModel(TableModel(result))
        
    
    def RSDData(self, msg):
        self.rsdData.append(msg)
        

    def ProtData(self, msg):
        self.protData.append(msg)
    
    
    def ValueData(self, msg):
        self.valueData.append(msg)
        

    def ProcessBar(self, msg):
        self.progressBar.setValue(int(msg))    
    
        
    def SaveData(self):
        options = QtWidgets.QFileDialog.Options()
        options |= QtWidgets.QFileDialog.DontUseNativeDialog
        fileName, _ = QtWidgets.QFileDialog.getSaveFileName(self,"QFileDialog.getOpenFileName()", "","All Files (*);;CSV Files (*.csv)", options=options)
        data = pd.DataFrame(self.tableView.model()._data)
        data.to_csv(fileName, index = False)
        
        

if __name__ == '__main__':
    import sys
    
    app = QApplication(sys.argv)
    ui = PreprocessUI()
    ui.show()
    sys.exit(app.exec_())