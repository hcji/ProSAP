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
from ParamTSAUI import ParamTSAUI
from RunningUI import Running_Win

from Thread import TPPThread, NPAThread, DistThread
from MakeFigure import MakeFigure
from Utils import TableModel, ReplicateCheck
from Inflect import run_inflect

'''
r1p1Data = pd.read_csv('D:/project/CETSA_Benchmark/Data/Ball_STS/DMSO_1.csv')
r1p2Data = pd.read_csv('D:/project/CETSA_Benchmark/Data/Ball_STS/Stauro_1.csv')
r2p1Data = pd.read_csv('D:/project/CETSA_Benchmark/Data/Ball_STS/DMSO_2.csv')
r2p2Data = pd.read_csv('D:/project/CETSA_Benchmark/Data/Ball_STS/Stauro_2.csv')
columns = ['T37', 'T41.2', 'T44', 'T46.8', 'T50', 'T53.2', 'T54', 'T56.1', 'T59.1', 'T63.2', 'T66.9']
'''

class AnalTSAUI(QtWidgets.QMainWindow, Ui_MainWindow):
    
    def __init__(self, parent=None): 
        super(AnalTSAUI, self).__init__(parent)
        self.setupUi(self)
        self.setWindowTitle("TPP Analysis")
        self.setWindowIcon(QtGui.QIcon("img/TPCA.ico"))
        
        # main window
        self.resize(1500, 900)
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
        self.TPPThread = None
        self.NPARCThread = None
        self.DistThread = None
        
        # menu action
        self.actionProteomics.triggered.connect(self.LoadProteinFile)
        self.actionTPP.triggered.connect(self.ShowAnalTPP)
        self.actionNPARC.triggered.connect(self.ShowAnalNPARC)
        self.actionDistance.triggered.connect(self.ShowAnalDist)
        self.actionINFLECT.triggered.connect(self.ShowAnalInflect)
        
        # button action
        self.tableWidgetProteinList.setSortingEnabled(True)
        self.ButtonR1P1.clicked.connect(self.SetR1P1)
        self.ButtonR1P2.clicked.connect(self.SetR1P2)
        self.ButtonR2P1.clicked.connect(self.SetR2P1)
        self.ButtonR2P2.clicked.connect(self.SetR2P2)
        self.ButtonRemove.clicked.connect(self.RemoveProteinFile)
        self.ButtonClearFileList.clicked.connect(self.ClearProteinFile)
        
        self.ButtonParam.clicked.connect(self.OpenParams)
        self.ButtonShow.clicked.connect(self.ShowMeltCurve)
        self.ButtonSave.clicked.connect(self.SaveData)
        
        self.ColumnSelectUI = ColumnSelectUI()
        self.ColumnSelectUI.ButtonColumnSelect.clicked.connect(self.SetProteinColumn)
        self.ColumnSelectUI.ButtonColumnCancel.clicked.connect(self.ColumnSelectUI.close)

        self.ParamTSAUI = ParamTSAUI()
        self.ParamTSAUI.ButtonConfirm.clicked.connect(self.SetParams)
        self.ParamTSAUI.ButtonCancel.clicked.connect(self.ParamTSAUI.close)
        
        self.Running_Win = Running_Win()
        
        # server data
        self.columns = None
        self.prots = None
        self.resultTable = pd.DataFrame()
        self.resultData = []
        
        # default params
        self.haxis_TPP = 0.5
        self.minR2_TPP = 0.8
        self.maxPlateau_TPP = 0.3
        self.repCheck_TPP = 'True'
        self.minR2_NP_Null = 0.8
        self.minR2_NP_Alt = 0.8
        self.maxPlateau_NP = 0.3
        self.minR2_Infl = 0.8
        self.numSD_Infl = 2
        self.Metr_Dist = 'cityblock'
        self.minR2_Dist = 0.8
        self.maxPlateau_Dist = 0.3
        

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
    
    
    def OpenParams(self):
        self.ParamTSAUI.show()
    
    
    def SetParams(self):
        self.haxis_TPP = self.ParamTSAUI.BoxHAxis.value()
        self.minR2_TPP = self.ParamTSAUI.BoxR2.value()
        self.maxPlateau_TPP = self.ParamTSAUI.BoxPla.value()
        self.repCheck_TPP = self.ParamTSAUI.BoxCheck.currentText()
        self.minR2_NP_Null = self.ParamTSAUI.BoxR2_Null.value()
        self.minR2_NP_Alt = self.ParamTSAUI.BoxR2_Alt.value()
        self.maxPlateau_NP = self.ParamTSAUI.BoxPlaN.value()
        self.minR2_Infl = self.ParamTSAUI.BoxR2_Infl.value()
        self.numSD_Infl = self.ParamTSAUI.BoxNumSD.value()
        self.Metr_Dist = self.ParamTSAUI.BoxMetrics.currentText()
        self.minR2_Dist = self.ParamTSAUI.BoxR2_Dist.value()
        self.maxPlateau_Dist = self.ParamTSAUI.BoxPla_Dist.value()
        self.ParamTSAUI.close()
    
    
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
        
    
    def RemoveProteinFile(self):
          listItems = self.ListFile.selectedItems()
          if not listItems: return
          for item in listItems:
              self.ListFile.takeItem(self.ListFile.row(item))

        
    def FillTable(self, resultTable):
        self.tableWidgetProteinList.setRowCount(resultTable.shape[0])
        self.tableWidgetProteinList.setColumnCount(resultTable.shape[1])
        self.tableWidgetProteinList.setHorizontalHeaderLabels(resultTable.columns)
        self.tableWidgetProteinList.setVerticalHeaderLabels(resultTable.index.astype(str))
        for i in range(resultTable.shape[0]):
            for j in range(resultTable.shape[1]):
                if type(resultTable.iloc[i,j]) == np.float64:
                    item = QtWidgets.QTableWidgetItem()
                    item.setData(Qt.EditRole, QVariant(float(resultTable.iloc[i,j])))
                    # item = QtWidgets.QTableWidgetItem(str(resultTable.iloc[i,j]))
                else:
                    item = QtWidgets.QTableWidgetItem(str(resultTable.iloc[i,j]))
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


    def ProcessBar(self, msg):
        self.progressBar.setValue(int(msg))
        
    
    def ResultData(self, msg):
        self.resultData.append(msg)
        # print(msg)
    
    
    def DisableMenu(self):
        self.actionTPP.setEnabled(False)
        self.actionNPARC.setEnabled(False)
        self.actionDistance.setEnabled(False)
        self.actionINFLECT.setEnabled(False)
        
        
    def EnableMenu(self):
        self.actionTPP.setEnabled(True)
        self.actionNPARC.setEnabled(True)
        self.actionDistance.setEnabled(True)
        self.actionINFLECT.setEnabled(True)        
    
    
    # TPP Analysis 
    def ShowAnalTPP(self):
        self.DisableMenu()
        
        columns = self.columns
        self.tableWidgetProteinList.clear()
        self.progressBar.setValue(0)
        self.resultData = []
        
        r1p1Data = self.tableRep1Protein1.model()._data
        r1p2Data = self.tableRep1Protein2.model()._data
        try:
            r2p1Data = self.tableRep2Protein1.model()._data
            r2p2Data = self.tableRep2Protein2.model()._data
        except:
            r2p1Data = None
            r2p2Data = None

        haxis = self.haxis_TPP
        minR2 = self.minR2_TPP
        maxPlateau = self.maxPlateau_TPP
        
        temps = np.array([float(t.replace('T', '')) for t in columns])
        cols = ['Accession'] + columns
        
        r1p1 = r1p1Data.loc[:, cols]
        r1p2 = r1p2Data.loc[:, cols]
        if (r2p1Data is not None) and (r2p2Data is not None):
            r2p1 = r2p1Data.loc[:, cols]
            r2p2 = r2p2Data.loc[:, cols]
            prot_1 = np.intersect1d(list(r1p1.iloc[:,0]), list(r1p2.iloc[:,0]))
            prot_2 = np.intersect1d(list(r2p1.iloc[:,0]), list(r2p2.iloc[:,0]))
            self.prots = np.intersect1d(prot_1, prot_2)
        else:
            r2p1 = None
            r2p2 = None
            self.prots = np.intersect1d(list(r1p1.iloc[:,0]), list(r1p2.iloc[:,0]))
        
        self.TPPThread = TPPThread(self.prots, temps, r1p1, r1p2, r2p1, r2p2, minR2, maxPlateau, haxis)
        self.TPPThread._ind.connect(self.ProcessBar)
        self.TPPThread._res.connect(self.ResultData)
        self.TPPThread.start()
        self.TPPThread.finished.connect(self.VisualizeTPP)
    
    
    def VisualizeTPP(self):   
        prots = self.prots
        r2p1Data = self.tableRep2Protein1.model()
        r2p2Data = self.tableRep2Protein2.model()
        
        res = pd.DataFrame(self.resultData)
        if (r2p1Data is None) or (r2p2Data is None):
            res.columns = ['Rep1Group1_R2', 'Rep1Group2_R2', 'Rep1Group1_Tm', 'Rep1Group2_Tm', 'Rep1delta_Tm', 'Rep1min_Slope']
        else:
            res.columns = ['Rep1Group1_R2', 'Rep1Group2_R2', 'Rep1Group1_Tm', 'Rep1Group2_Tm', 'Rep1delta_Tm', 'Rep1min_Slope',
                           'Rep2Group1_R2', 'Rep2Group2_R2', 'Rep2Group1_Tm', 'Rep2Group2_Tm', 'Rep2delta_Tm', 'Rep2min_Slope']
    
        delta_Tm = res['Rep1delta_Tm']
        p_Val = []
        for i in range(len(res)):
            s = delta_Tm[i]
            pv = stats.t.sf(abs(s - np.nanmean(delta_Tm)) / np.nanstd(delta_Tm), len(delta_Tm)-1)
            p_Val.append(pv)
        res['Rep1pVal (-log10)'] = np.round(-np.log10(p_Val), 3)
        score = -np.log10(np.array(p_Val)) * (res['Rep1Group1_R2'] * res['Rep1Group2_R2']) ** 2
            
        if (r2p1Data is not None) and (r2p2Data is not None):
            delta_Tm = res['Rep2delta_Tm']
            p_Val = []
            for i in range(len(res)):
                s = delta_Tm[i]
                pv = stats.t.sf(abs(s - np.nanmean(delta_Tm)) / np.nanstd(delta_Tm), len(delta_Tm)-1)
                p_Val.append(pv)
            res['Rep2pVal (-log10)'] = np.round(-np.log10(p_Val), 3)
            score_2 = -np.log10(np.array(p_Val)) * (res['Rep2Group1_R2'] * res['Rep2Group2_R2']) ** 2
            score += score_2
        
        score[np.isnan(score)] = 0
        res['Score'] = score
        res = np.round(res, 3)
        res['Accession'] = prots
        if (r2p1Data is None) or (r2p2Data is None):
            res = res[['Accession', 'Score', 'Rep1pVal (-log10)', 'Rep1delta_Tm', 'Rep1Group1_R2', 'Rep1Group2_R2', 'Rep1Group1_Tm', 'Rep1Group2_Tm', 'Rep1min_Slope']]
        else:
            res = res[['Accession', 'Score', 'Rep1pVal (-log10)', 'Rep1delta_Tm', 'Rep1Group1_R2', 'Rep1Group2_R2', 'Rep1Group1_Tm', 'Rep1Group2_Tm', 'Rep1min_Slope',
                       'Rep2pVal (-log10)', 'Rep2delta_Tm', 'Rep2Group1_R2', 'Rep2Group2_R2', 'Rep2Group1_Tm', 'Rep2Group2_Tm', 'Rep2min_Slope']]
            if self.repCheck_TPP == 'True':
                res = ReplicateCheck(res)
                
        resultTable = res.sort_values(by = 'Score', axis = 0, ascending = False)
        self.resultData = []
        self.resultTable = resultTable
        self.FillTable(resultTable)
        self.EnableMenu()
        

    # NPARC Analysis 
    def ShowAnalNPARC(self):
        self.DisableMenu()
        
        columns = self.columns
        self.tableWidgetProteinList.clear()
        self.progressBar.setValue(0)
        self.resultData = []
        
        r1p1Data = self.tableRep1Protein1.model()._data
        r1p2Data = self.tableRep1Protein2.model()._data
        try:
            r2p1Data = self.tableRep2Protein1.model()._data
            r2p2Data = self.tableRep2Protein2.model()._data
        except:
            r2p1Data = None
            r2p2Data = None
        
        minR2_null = self.minR2_NP_Null
        minR2_alt = self.minR2_NP_Alt
        maxPlateau = self.maxPlateau_NP
        
        temps = np.array([float(t.replace('T', '')) for t in columns])
        cols = ['Accession'] + columns
        
        if (r2p1Data is None) or (r2p2Data is None):
            self.ErrorMsg('NPARC must run with two replicates')
            return None
        r1p1 = r1p1Data.loc[:, cols]
        r1p2 = r1p2Data.loc[:, cols]
        r2p1 = r2p1Data.loc[:, cols]
        r2p2 = r2p2Data.loc[:, cols]
        
        prot_1 = np.intersect1d(list(r1p1.iloc[:,0]), list(r1p2.iloc[:,0]))
        prot_2 = np.intersect1d(list(r2p1.iloc[:,0]), list(r2p2.iloc[:,0]))
        self.prots = np.intersect1d(prot_1, prot_2)
        
        self.NPARCThread = NPAThread(self.prots, temps, r1p1, r1p2, r2p1, r2p2, minR2_null, minR2_alt, maxPlateau)
        self.NPARCThread._ind.connect(self.ProcessBar)
        self.NPARCThread._res.connect(self.ResultData)
        self.NPARCThread.start()
        self.NPARCThread.finished.connect(self.VisualizeNPARC)


    def VisualizeNPARC(self):
        prots = self.prots
        res = pd.DataFrame(self.resultData)
        res.columns = ['Group1_R2', 'Group2_R2', 'RSS_Null', 'RSS_Alt', 'RSS_Diff']
        
        d1 = 3
        d2 = len(self.columns) * 4 - 6
        RSS_Diff = res['RSS_Diff'].values
        RSS_Alt = res['RSS_Alt'].values

        p_Val = []
        for i in range(len(res)):
            s = (RSS_Diff[i] / d1) / (RSS_Alt[i] / d2)
            pv = stats.f.sf(s, d1, d2)
            p_Val.append(pv)
        score = -np.log10(np.array(p_Val)) * (res['Group1_R2'] * res['Group2_R2']) ** 2
        
        score[np.isnan(score)] = 0
        res['Accession'] = prots
        res['p_Val (-log10)'] = -np.log10(p_Val)
        res['Score'] = score
        res = np.round(res, 3)        
        
        res = res[['Accession', 'Score', 'p_Val (-log10)', 'Group1_R2', 'Group2_R2', 'RSS_Null', 'RSS_Alt', 'RSS_Diff']]
        resultTable = res.sort_values(by = 'Score', axis = 0, ascending = False)
        self.resultData = []
        self.resultTable = resultTable
        self.FillTable(resultTable)        
        self.EnableMenu()
        
    
    # Distance methods
    def ShowAnalDist(self):
        self.DisableMenu()
        
        columns = self.columns
        self.tableWidgetProteinList.clear()
        self.progressBar.setValue(0)
        self.resultData = []
        
        r1p1Data = self.tableRep1Protein1.model()._data
        r1p2Data = self.tableRep1Protein2.model()._data
        try:
            r2p1Data = self.tableRep2Protein1.model()._data
            r2p2Data = self.tableRep2Protein2.model()._data
        except:
            r2p1Data = None
            r2p2Data = None

        method = self.Metr_Dist
        minR2 = self.minR2_Dist
        maxPlateau = self.maxPlateau_Dist
        
        temps = np.array([float(t.replace('T', '')) for t in columns])
        cols = ['Accession'] + columns
        
        r1p1 = r1p1Data.loc[:, cols]
        r1p2 = r1p2Data.loc[:, cols]
        if (r2p1Data is not None) and (r2p2Data is not None):
            r2p1 = r2p1Data.loc[:, cols]
            r2p2 = r2p2Data.loc[:, cols]
            prot_1 = np.intersect1d(list(r1p1.iloc[:,0]), list(r1p2.iloc[:,0]))
            prot_2 = np.intersect1d(list(r2p1.iloc[:,0]), list(r2p2.iloc[:,0]))
            self.prots = np.intersect1d(prot_1, prot_2)
        else:
            r2p1 = None
            r2p2 = None
            self.prots = np.intersect1d(list(r1p1.iloc[:,0]), list(r1p2.iloc[:,0]))
        
        self.DistThread = DistThread(self.prots, temps, r1p1, r1p2, r2p1, r2p2, method, minR2, maxPlateau)
        self.DistThread._ind.connect(self.ProcessBar)
        self.DistThread._res.connect(self.ResultData)
        self.DistThread.start()
        self.DistThread.finished.connect(self.VisualizeDist)
    
    
    def VisualizeDist(self):
        self.DisableMenu()
        
        prots = self.prots
        r2p1Data = self.tableRep2Protein1.model()
        r2p2Data = self.tableRep2Protein2.model()
        
        res = pd.DataFrame(self.resultData)
        if (r2p1Data is None) or (r2p2Data is None):
            res.columns = ['Rep1Group1_R2', 'Rep1Group2_R2', 'Rep1Group1_SH', 'Rep1Group2_SH', 'Rep1delta_SH', 'Rep1min_Slope']
        else:
            res.columns = ['Rep1Group1_R2', 'Rep1Group2_R2', 'Rep1Group1_SH', 'Rep1Group2_SH', 'Rep1delta_SH', 'Rep1min_Slope',
                           'Rep2Group1_R2', 'Rep2Group2_R2', 'Rep2Group1_SH', 'Rep2Group2_SH', 'Rep2delta_SH', 'Rep2min_Slope']
        
        if 'Rep2delta_SH' in res.columns:
            sign_1 = np.sign(res['Rep1Group2_SH'] - res['Rep1Group1_SH'])
            sign_2 = np.sign(res['Rep2Group2_SH'] - res['Rep2Group1_SH'])
            delta_SH = np.abs(res['Rep1delta_SH'] * sign_1 + res['Rep2delta_SH'] * sign_2)
        else:
            delta_SH = res['Rep1delta_SH']
        
        p_Val = []
        for i in range(len(res)):
            s = delta_SH[i]
            pv = stats.t.sf(abs(s - np.nanmean(delta_SH)) / np.nanstd(delta_SH), len(delta_SH)-1)
            p_Val.append(pv)
        res['pVal (-log10)'] = np.round(-np.log10(p_Val), 3)
        score = -np.log10(np.array(p_Val)) * (res['Rep1Group1_R2'] * res['Rep1Group2_R2']) ** 2
        
        score[np.isnan(score)] = 0
        res['Score'] = score
        res = np.round(res, 3)
        res['Accession'] = prots
        if (r2p1Data is None) or (r2p2Data is None):
            res = res[['Accession', 'Score', 'pVal (-log10)', 'Rep1delta_SH', 'Rep1Group1_R2', 'Rep1Group2_R2', 'Rep1Group1_SH', 'Rep1Group2_SH', 'Rep1min_Slope']]
        else:
            res = res[['Accession', 'Score', 'pVal (-log10)', 'Rep1delta_SH', 'Rep2delta_SH', 'Rep1Group1_R2', 'Rep1Group2_R2', 'Rep1Group1_SH', 'Rep1Group2_SH', 'Rep1min_Slope',
                       'Rep2Group1_R2', 'Rep2Group2_R2', 'Rep2Group1_SH', 'Rep2Group2_SH', 'Rep2min_Slope']]
            if self.repCheck_TPP == 'True':
                res = ReplicateCheck(res)  
        resultTable = res.sort_values(by = 'Score', axis = 0, ascending = False)
        self.resultData = []
        self.resultTable = resultTable
        self.FillTable(resultTable)
        self.EnableMenu()
        
        
    # Inflect method
    def ShowAnalInflect(self):
        columns = self.columns
        self.tableWidgetProteinList.clear()
        self.progressBar.setValue(0)
        self.resultData = []
        
        r1p1Data = self.tableRep1Protein1.model()._data
        r1p2Data = self.tableRep1Protein2.model()._data
        try:
            r2p1Data = self.tableRep2Protein1.model()._data
            r2p2Data = self.tableRep2Protein2.model()._data
        except:
            r2p1Data = None
            r2p2Data = None
        
        Rsq = self.minR2_Infl
        NumSD = self.numSD_Infl
        
        temps = np.array([float(t.replace('T', '')) for t in columns])
        cols = ['Accession'] + columns
        
        if (r2p1Data is None) or (r2p2Data is None):
            self.ErrorMsg('NPARC must run with two replicates')
            return None
        r1p1 = r1p1Data.loc[:, cols]
        r1p2 = r1p2Data.loc[:, cols]
        r2p1 = r2p1Data.loc[:, cols]
        r2p2 = r2p2Data.loc[:, cols]
        
        self.Running_Win.show()
        run_inflect(temps, r1p1, r1p2, r2p1, r2p2, Rsq, NumSD)
        self.Running_Win.close()
        self.VisualizeInflect()
    
    
    def VisualizeInflect(self):
        try:
            Rep1Result = pd.read_excel('C:/inflect_tempdir/Rep 1/Results.xlsx')
            Rep2Result = pd.read_excel('C:/inflect_tempdir/Rep 2/Results.xlsx')
            Rep1Sig = pd.read_excel('C:/inflect_tempdir/Rep 1/SignificantResults.xlsx')
            Rep2Sig = pd.read_excel('C:/inflect_tempdir/Rep 2/SignificantResults.xlsx')
        except:
            self.ErrorMsg('Unexpected error, try other method')
            return None
        
        Rep1Result = Rep1Result[['Accession', 'R2..Control', 'R2..Condition', 'Tm..Control', 'Tm..Condition']]
        Rep1Result.columns = ['Accession', 'Rep1Group1_R2', 'Rep1Group2_R2', 'Rep1Group1_Tm', 'Rep1Group2_Tm']
        
        Rep2Result = Rep2Result[['Accession', 'R2..Control', 'R2..Condition', 'Tm..Control', 'Tm..Condition']]
        Rep2Result.columns = ['Accession', 'Rep2Group1_R2', 'Rep2Group2_R2', 'Rep2Group1_Tm', 'Rep2Group2_Tm']
        
        all_sig = set(Rep1Sig['Accession']) & set(Rep2Sig['Accession'])
        
        res = Rep1Result.merge(Rep2Result, on = 'Accession')
        Significant = [str(x in all_sig) for x in res['Accession']]
        res['Significant'] = Significant
        self.resultTable = res.sort_values(by = 'Significant', axis = 0, ascending = False)
        
        self.FillTable(self.resultTable)
        self.EnableMenu()
    

    # Common functions
    def ShowMeltCurve(self):
        columns = self.columns
        r1p1Data = self.tableRep1Protein1.model()._data
        r1p2Data = self.tableRep1Protein2.model()._data
        try:
            r2p1Data = self.tableRep2Protein1.model()._data
            r2p2Data = self.tableRep2Protein2.model()._data
        except:
            r2p1Data = None
            r2p2Data = None
        header = [self.tableWidgetProteinList.horizontalHeaderItem(i).text() for i in range(self.tableWidgetProteinList.columnCount())]
        i = self.tableWidgetProteinList.selectedItems()[0].row()
        j = header.index('Accession')
        ProteinAccession = self.tableWidgetProteinList.item(i, j).text()

        self.figureTSA1.SingleTSAFigure(r1p1Data, r1p2Data, columns, ProteinAccession)
        if (r2p1Data is not None) and (r2p2Data is not None):
            self.figureTSA2.SingleTSAFigure(r2p1Data, r2p2Data, columns, ProteinAccession)
    

    def SaveData(self):
        options = QtWidgets.QFileDialog.Options()
        options |= QtWidgets.QFileDialog.DontUseNativeDialog
        fileName, _ = QtWidgets.QFileDialog.getSaveFileName(self, "Save", ".csv","CSV Files (*.csv)", options=options)
        if fileName:
            data = self.resultTable
            data.to_csv(fileName)    
        

        

if __name__ == '__main__':
    import sys
    
    app = QApplication(sys.argv)
    ui = AnalTSAUI()
    ui.show()
    sys.exit(app.exec_())