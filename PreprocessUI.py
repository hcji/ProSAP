# -*- coding: utf-8 -*-
"""
Created on Thu Apr 15 15:57:14 2021

@author: hcji
"""


import numpy as np
import pandas as pd

from scipy.optimize import curve_fit
from sklearn.impute import KNNImputer
from PyQt5 import QtCore, QtGui, QtWidgets
from PyQt5.QtWidgets import QApplication, QMainWindow, QGridLayout
from matplotlib.backends.backend_qt5 import NavigationToolbar2QT as NavigationToolbar

from Preprocess import Ui_Form
from ColumnSelectUI import ColumnSelectUI
from Utils import TableModel
from Thread import PreprocessThread
from MakeFigure import MakeFigure
from Utils import fit_curve, meltCurve


class PreprocessUI(QtWidgets.QWidget, Ui_Form):
    
    def __init__(self, parent=None): 
        super(PreprocessUI, self).__init__(parent)
        self.setupUi(self)
        self.setWindowTitle("Preprocessing")
        self.setWindowIcon(QtGui.QIcon("img/TPCA.ico"))
        
        self.comboBoxMerging.addItems(['Median', 'Mean'])
        self.comboBoxNorm.addItems(['Reference', 'Median', 'None'])
        self.comboBoxMV.addItems(['None', 'KNN', 'Zero'])
        self.comboBoxPSM.addItems(['None'])

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

    
    def ClearProteinFile(self):
        self.ListFile.clear()
    
        
    def LoadProteinFile(self):
        self.ListFile.clear()
        self.comboBoxReference.clear()
        
        options = QtWidgets.QFileDialog.Options()
        options |= QtWidgets.QFileDialog.DontUseNativeDialog
        fileNames, _ = QtWidgets.QFileDialog.getOpenFileNames(self,"Load", "","All Files (*);;CSV Files (*.csv)", options=options)
        
        if len(fileNames) == 0:
            pass
        
        else:
            for fileName in fileNames:
                if fileName.split('.')[-1] in ['csv', 'xls', 'xlsx']:
                    self.ListFile.addItem(fileName)
                else:
                    self.WarnMsg("Find file with invalid format, only csv and xlsx are support")
                
            if self.ListFile.count() == 0:
                pass
            
            else:
                if fileName.split('.')[1] == 'csv':
                    data = pd.read_csv(fileName)
                else:
                    data = pd.read_excel(fileName)
        
                all_cols = data.columns
                self.comboBoxPSM.addItems(all_cols)
                
                self.ColumnSelectUI.listWidget.clear()
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
        # import os
        # d = 'D:/project/PDFiles'
        # fileNames = os.listdir(d)
        # fileNames = [d+'/'+f for f in fileNames]
        fileNames = [self.ListFile.item(i).text() for i in range(self.ListFile.count())]
        fileNames_ = []
        if len(fileNames) == 0:
            return None
        
        for f in fileNames:
            if f.split('.')[1] in ['csv', 'xlsx', 'xls']:
                fileNames_.append(f)
        fileNames = fileNames_
        
        columns = [i.text() for i in self.ColumnSelectUI.listWidget.selectedItems()]
        for c in columns:
            if '--' in str(c):
                msg = QtWidgets.QMessageBox()
                msg.resize(550, 200)
                msg.setIcon(QtWidgets.QMessageBox.Information)
                msg.setText("'--' cannot in column names")
                msg.setWindowTitle("Information")
                msg.exec_()
                return None
        
        if len(columns) == 0:
            return None
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
        mv_thres = self.doubleSpinBoxMVFilter.value()
            
        all_data = []
        for f in fileNames:
            if f.split('.')[1]  == 'csv':
                data = pd.read_csv(f)
            else:
                data = pd.read_excel(f)
            
            if 'Accession' in data.columns:
                for j in columns:
                    data.loc[:,j] = pd.to_numeric(data.loc[:,j], errors='coerce')
                all_data.append(data)
                # print(data)
                
        if len(all_data) == 0:
            msg = QtWidgets.QMessageBox()
            msg.resize(550, 200)
            msg.setIcon(QtWidgets.QMessageBox.Critical)
            msg.setText("No valid file input")
            msg.setWindowTitle("Error")
            msg.exec_()
            self.ListFile.clear()
            return None         
        
        # devided by reference
        if self.comboBoxNorm.currentText() == 'Reference':
            for i in range(len(all_data)):
                ref = all_data[i].loc[:,reference].copy() + 0.0001
                for c in columns:
                    all_data[i].loc[:,c] /= ref
                if not self.tempertures:
                    meds = np.nanmedian(all_data[i], axis = 0)
                    paras = curve_fit(meltCurve, self.tempertures, meds, bounds=(0, [15000, 250, 0.3]))[0]
                    meds_ = meltCurve(self.tempertures, paras[0], paras[1], paras[2])
                    norm_ = meds_ - meds
                    all_data[i] = all_data[i] + norm_
                    
        
        data = all_data[0]
        if len(all_data) == 1:
            pass
        else:
            for data_ in all_data[1:]:
                data = data.merge(data_, 'outer', suffixes=('--x', '--y'), on = 'Accession')
        
        if self.comboBoxMerging.currentText() == 'Median':
            fun = np.nanmedian
        else:
            fun = np.nanmean
            
        self.valueData, self.rsdData, self.protData = [], [], []
            
        # normalize by median
        if self.comboBoxNorm.currentText() == 'Median':
            whs = np.where([c.split('--')[0] in columns for c in list(data.columns)])[0]
            medians = np.nanmedian(data.iloc[:,whs], axis = 0)
            ratios = medians / np.mean(medians)
            data.iloc[:,whs] = data.iloc[:,whs] / ratios
            
        self.PreprocessThread = PreprocessThread(data, psm_column, psm_thres, std_thres, columns, fun, mv_thres)
        self.PreprocessThread._ind.connect(self.ProcessBar)
        self.PreprocessThread._val.connect(self.ValueData)
        self.PreprocessThread._rsd.connect(self.RSDData)
        self.PreprocessThread._prot.connect(self.ProtData)
        self.PreprocessThread.start()
        self.PreprocessThread.finished.connect(self.VisualizeProprocess)            
    
    
    def VisualizeProprocess(self):
        val_list = pd.DataFrame(self.valueData)
        if self.comboBoxMV.currentText() == 'KNN':
            knn_imputer = KNNImputer(n_neighbors = 2)
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
        fileName, _ = QtWidgets.QFileDialog.getSaveFileName(self, "Save", ".csv","All Files (*);;CSV Files (*.csv)", options=options)
        if fileName:
            data = pd.DataFrame(self.tableView.model()._data)
            data.to_csv(fileName, index = False)
        
        

if __name__ == '__main__':
    import sys
    
    app = QApplication(sys.argv)
    ui = PreprocessUI()
    ui.show()
    sys.exit(app.exec_())