# -*- coding: utf-8 -*-
"""
Created on Thu Apr 15 15:57:14 2021

@author: hcji
"""

import numpy as np
import pandas as pd
from scipy import stats

from PyQt5.QtCore import Qt, QVariant
from PyQt5 import QtCore, QtGui, QtWidgets
from PyQt5.QtWidgets import QApplication, QMainWindow, QGridLayout
from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg as FigureCanvas
from matplotlib.backends.backend_qt5 import NavigationToolbar2QT as NavigationToolbar

from AnalTSA import Ui_MainWindow
from ColumnSelectUI import ColumnSelectUI
from MakeFigure import MakeFigure

from Thread import CurveFitThread, AnalDistThread
from MakeFigure import MakeFigure
from Utils import TableModel, fit_curve
from iTSA import estimate_df, p_value_adjust


class AnalTSAUI(QtWidgets.QMainWindow, Ui_MainWindow):
    
    def __init__(self, parent=None): 
        super(AnalTSAUI, self).__init__(parent)
        self.setupUi(self)
        self.setWindowTitle("TPP Analysis")
        self.setWindowIcon(QtGui.QIcon("img/TPCA.ico"))
        
        # main window
        self.resize(1300, 800)
        self.setMinimumWidth(1150)
        self.setMinimumHeight(650)
        self.move(75, 50)
        
        # figure box
        self.figureTSA1 = MakeFigure(10, 10, dpi = 250)
        self.figureTSA1_ntb = NavigationToolbar(self.figureTSA1, self)
        self.gridlayoutTSA = QGridLayout(self.groupBox1)
        self.gridlayoutTSA.addWidget(self.figureTSA1)
        self.gridlayoutTSA.addWidget(self.figureTSA1_ntb)
        
        self.figureTSA2 = MakeFigure(10, 10, dpi = 250)
        self.figureTSA2_ntb = NavigationToolbar(self.figureTSA2, self)
        self.gridlayoutTSA2 = QGridLayout(self.groupBox2)
        self.gridlayoutTSA2.addWidget(self.figureTSA2)
        self.gridlayoutTSA2.addWidget(self.figureTSA2_ntb)
        
        # threads
        self.CurveFitThread = None
        self.AnalDistThread = None
        
        # menu action
        self.actionProteomics.triggered.connect(self.LoadProteinFile)
        
        # button action
        self.ButtonR1P1.clicked.connect(self.SetR1P1)
        self.ButtonR1P2.clicked.connect(self.SetR1P2)
        self.ButtonR2P1.clicked.connect(self.SetR2P1)
        self.ButtonR2P2.clicked.connect(self.SetR2P2)
        
        self.ColumnSelectUI = ColumnSelectUI()
        self.ColumnSelectUI.ButtonColumnSelect.clicked.connect(self.SetProteinColumn)
        self.ColumnSelectUI.ButtonColumnCancel.clicked.connect(self.ColumnSelectUI.close)
        
        
        # server data
        self.columns = None
        self.prots = None
        self.TSA_table = pd.DataFrame()
        self.resultDataTSA = []
             

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
            for fileName in fileNames:
                if fileName:
                    if fileName.split('.')[1] in ['csv', 'xlsx']:
                        self.ListFile.addItem(fileName)
                    else:
                        pass


    def ClearProteinFile(self):
        self.ListFile.clear()

        
    def FillTable(self, TSA_table):
        self.tableWidgetProteinList.setRowCount(TSA_table.shape[0])
        self.tableWidgetProteinList.setColumnCount(TSA_table.shape[1])
        self.tableWidgetProteinList.setHorizontalHeaderLabels(TSA_table.columns)
        self.tableWidgetProteinList.setVerticalHeaderLabels(TSA_table.index.astype(str))
        for i in range(TSA_table.shape[0]):
            for j in range(TSA_table.shape[1]):
                if type(TSA_table.iloc[i,j]) == np.float64:
                    item = QtWidgets.QTableWidgetItem()
                    item.setData(Qt.EditRole, QVariant(float(TSA_table.iloc[i,j])))
                    # item = QtWidgets.QTableWidgetItem(str(TSA_table.iloc[i,j]))
                else:
                    item = QtWidgets.QTableWidgetItem(str(TSA_table.iloc[i,j]))
                self.tableWidgetProteinList.setItem(i, j, item)
                
                
    def ReplaceNonNumeric(self, data):
        for col in data.columns[1:]:
            data[col] = pd.to_numeric(data[col], errors='coerce')
        # data = data.astype(float)
        keep = np.logical_not(np.isnan(np.sum(np.array(data)[:,1:], axis = 1).astype(float)))
        data = data.iloc[keep,:]
        data = data.reset_index(drop = True)
        return data


    def SelectProteinTable(self):
        selectItem = self.ListFile.currentItem()
        if not selectItem:
            self.ErrorMsg('No item is selected')
            return None
        try:
            if selectItem.text().split('.')[1] == 'csv':
                selectData = pd.read_csv(selectItem.text())
            elif selectItem.text().split('.')[1] == 'xlsx':
                selectData = pd.read_excel(selectItem.text())
            else:
                return None
        except:
            self.ErrorMsg('Cannot be load the selected file')
        
        if 'Accession' in selectData.columns:
            return selectData
        else:
            self.ErrorMsg('Accession is not given in the data')
            return None


    def SetR1P1(self):
        data = self.SelectProteinTable()
        if data is None:
            return None
        self.ColumnSelectUI.listWidget.clear()
        all_cols = data.columns
        for c in all_cols:
            self.ColumnSelectUI.listWidget.addItem(c)
        self.ColumnSelectUI.show()


    def SetR1P2(self):
        if self.columns == None:
            self.ErrorMsg('Please set Group 1')
        else:
            data = self.SelectProteinTable()
            if data is None:
                return None
            columns = ['Accession'] + self.columns
            try:
                data = data.loc[:, columns]
                data = self.ReplaceNonNumeric(data)
                if np.nanmax(data.loc[:, self.columns]) > 10:
                    self.WarnMsg('The data seems not normalized')
                self.tableRep1Protein2.setModel(TableModel(data))
            except:
                self.ErrorMsg('No columns matched with Group 1')


    def SetR2P1(self):
        if self.columns == None:
            self.ErrorMsg('Please set Group 1')
        else:
            data = self.SelectProteinTable()
            if data is None:
                return None
            columns = ['Accession'] + self.columns
            try:
                data = data.loc[:, columns]
                data = self.ReplaceNonNumeric(data)
                if np.nanmax(data.loc[:, self.columns]) > 10:
                    self.WarnMsg('The data seems not normalized')
                self.tableRep2Protein1.setModel(TableModel(data))
            except:
                self.ErrorMsg('No columns matched with Replicate 1')


    def SetR2P2(self):
        if self.columns == None:
            self.ErrorMsg('Please set Group 1')
        else:
            data = self.SelectProteinTable()
            if data is None:
                return None
            columns = ['Accession'] + self.columns
            try:
                data = data.loc[:, columns]
                data = self.ReplaceNonNumeric(data)
                if np.nanmax(data.loc[:, self.columns]) > 10:
                    self.WarnMsg('The data seems not normalized')
                self.tableRep2Protein2.setModel(TableModel(data))
            except:
                self.ErrorMsg('No columns matched with Replicate 1')


    def SetProteinColumn(self):
        data = self.SelectProteinTable()
        self.columns = [i.text() for i in self.ColumnSelectUI.listWidget.selectedItems()]
        try:
            [float(t.replace('T', '')) for t in self.columns]
        except:
            self.columns = None
            self.ErrorMsg('Selected columns can only be Txx, where xx is a number representing temperature')
            return None
        columns = ['Accession'] + self.columns
        data = data.loc[:, columns]
        data = self.ReplaceNonNumeric(data)
        if np.nanmax(data.loc[:, self.columns] > 10):
            self.WarnMsg('The data seems not normalized')
        self.ColumnSelectUI.close()
        self.tableRep1Protein1.setModel(TableModel(data))


    # TPP Analysis  
    
    def ProcessBarTPP(self, msg):
        self.progressBar.setValue(int(msg))
        
    
    def ResultDataTPP(self, msg):
        self.resultDataTSA.append(msg)
        # print(msg)
    
    
    def ShowAnalTPP(self):
        columns = self.columns
        self.ColumnSelectUI.close()
        self.tableWidgetProteinList.clear()
        self.progressBar.setValue(0)
        
        proteinData1 = self.tableProtein1.model()._data
        proteinData2 = self.tableProtein2.model()._data

        h_axis = self.AnalTSAUI.Boxhaxis.value()
        minR2 = self.AnalTSAUI.BoxR2.value()
        maxPlateau = self.AnalTSAUI.BoxPlateau.value()
        
        temps = np.array([float(t.replace('T', '')) for t in columns])
        cols = ['Accession'] + columns
        data_1 = proteinData1.loc[:, cols]
        data_2 = proteinData2.loc[:, cols]

        self.prots = np.intersect1d(list(data_1.iloc[:,0]), list(data_2.iloc[:,0]))
        
        self.CurveFitThread = CurveFitThread(self.prots, temps, data_1, data_2, minR2, maxPlateau, h_axis)
        self.CurveFitThread._ind.connect(self.ProcessBarTSA)
        self.CurveFitThread._res.connect(self.ResultDataTSA)
        self.CurveFitThread.start()
        self.CurveFitThread.finished.connect(self.VisualizeTSA)
        
        '''
        res = []
        for i, p in enumerate(prots):
            x = temps
            y1 = np.array(data_1[data_1.iloc[:,0] == p].iloc[0,1:])
            y2 = np.array(data_2[data_2.iloc[:,0] == p].iloc[0,1:])
            res.append(fit_curve(x, y1, y2, minR2, maxPlateau, h_axis))
            self.AnalTSAUI.progressBar.setValue(int(i / len(prots)))
        '''
        # res = Parallel(n_jobs=n_core, backend = 'threading')(delayed(fit_curve)(p) for p in prots)
        # res = pd.DataFrame(res)
    
    
    def VisualizeTPP(self):   
        prots = self.prots
        res = pd.DataFrame(self.resultDataTSA)
        res.columns = ['Group1_R2', 'Group2_R2', 'Group1_Tm', 'Group2_Tm', 'delta_Tm', 'min_Slope']
    
        delta_Tm = res['delta_Tm']
        p_Val = []
        for i in range(len(res)):
            s = delta_Tm[i]
            pv = stats.t.sf(abs(s - np.nanmean(delta_Tm)) / np.nanstd(delta_Tm), len(delta_Tm)-1)
            p_Val.append(pv)
        score = -np.log10(np.array(p_Val)) * (res['Group1_R2'] * res['Group2_R2']) ** 2
    
        res['Accession'] = prots
        res['delta_Tm'] = delta_Tm
        res['p_Val (-log10)'] = -np.log10(p_Val)
        res['Score'] = score
        res = np.round(res, 3)
    
        res = res[['Accession', 'Score', 'p_Val (-log10)', 'delta_Tm', 'Group1_R2', 'Group2_R2', 'Group1_Tm', 'Group2_Tm', 'min_Slope']]
        TSA_table = res.sort_values(by = 'Score', axis = 0, ascending = False)
        
        self.resultDataTSA = []
        self.TSA_table = TSA_table
        self.AnalTSAUI.ButtonConfirm.setEnabled(True)
        self.AnalTSAUI.FillTable(TSA_table)

    
    
    
    
    
    def ShowMeltCurve(self):
        columns = self.columns
        proteinData1 = self.tableProtein1.model()._data
        proteinData2 = self.tableProtein2.model()._data
        
        header = [self.AnalTSAUI.tableWidgetProteinList.horizontalHeaderItem(i).text() for i in range(self.AnalTSAUI.tableWidgetProteinList.columnCount())]
        i = self.AnalTSAUI.tableWidgetProteinList.selectedItems()[0].row()
        j = header.index('Accession')
        ProteinAccession = self.AnalTSAUI.tableWidgetProteinList.item(i, j).text()

        self.AnalTSAUI.figureTSA.SingleTSAFigure(proteinData1, proteinData2, columns, ProteinAccession)
        

    def SaveData(self):
        options = QtWidgets.QFileDialog.Options()
        options |= QtWidgets.QFileDialog.DontUseNativeDialog
        fileName, _ = QtWidgets.QFileDialog.getSaveFileName(self, "Save", ".csv","CSV Files (*.csv)", options=options)
        if fileName:
            data = self.TSA_table
            data.to_csv(fileName)    
        

        

if __name__ == '__main__':
    import sys
    
    app = QApplication(sys.argv)
    ui = AnalTSAUI()
    ui.show()
    sys.exit(app.exec_())