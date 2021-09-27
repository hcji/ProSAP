# -*- coding: utf-8 -*-
"""
Created on Wed Apr 14 09:29:07 2021

@author: hcji
"""

# import ctypes
# ctypes.windll.shell32.SetCurrentProcessExplicitAppUserModelID("TPCA")

import numpy as np
import pandas as pd
from platform import system

from scipy import stats
from sklearn import metrics

from PyQt5.QtCore import Qt, QVariant
from PyQt5 import QtCore, QtGui, QtWidgets
from PyQt5.QtWidgets import QApplication, QMainWindow, QGridLayout
from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg as FigureCanvas
from matplotlib.backends.backend_qt5 import NavigationToolbar2QT as NavigationToolbar

from AnalTPCA import Ui_MainWindow
from ColumnSelectUI import ColumnSelectUI
from AnalROCUI import AnalROCUI
from Thread import ROCThread, PairThread, ComplexThread
from MakeFigure import MakeFigure
from Utils import TableModel, fit_curve

# proteinData1 = pd.read_csv('data/TPCA_TableS14_DMSO.csv')
# proteinData2 = pd.read_csv('data/TPCA_TableS14_MTX.csv')
# columns = ['T37', 'T40', 'T43', 'T46', 'T49', 'T52', 'T55', 'T58', 'T61', 'T64']
# proteinPair = pd.read_csv('data/TPCA_TableS2_Protein_Pairs.csv')
# proteinComplex = pd.read_csv('data/TPCA_TableS3_Protein_Complex.csv')


class AnalTPCAUI(QMainWindow, Ui_MainWindow):
    
    def __init__(self, parent=None):
        super(AnalTPCAUI, self).__init__(parent)
        self.setupUi(self)
        
        # main window
        self.resize(1300, 800)
        self.setMinimumWidth(1150)
        self.setMinimumHeight(650)
        self.move(75, 50)
        self.setWindowTitle("TPCA Analysis")
        self.setWindowIcon(QtGui.QIcon("img/TPCA.ico"))
        
        # fix macos menu
        if system() == 'Darwin':
            menubar = self.menuBar()
            menubar.setNativeMenuBar(False)
        
        # threads
        self.CurveFitThread = None
        self.ROCThread = None
        self.PairThread = None
        self.ComplexThread = None
        self.AnalDistThread = None
        
        # groupbox
        self.figureG1 = MakeFigure(10, 10)
        self.figureG1_ntb = NavigationToolbar(self.figureG1, self)
        self.gridlayoutG1 = QGridLayout(self.groupBox)
        self.gridlayoutG1.addWidget(self.figureG1)
        self.gridlayoutG1.addWidget(self.figureG1_ntb)
        
        self.figureG2 = MakeFigure(10, 10)
        self.figureG2_ntb = NavigationToolbar(self.figureG2, self)
        self.gridlayoutG2 = QGridLayout(self.groupBox_2)
        self.gridlayoutG2.addWidget(self.figureG2)
        self.gridlayoutG2.addWidget(self.figureG2_ntb)
        
        # widgets
        self.ColumnSelectUI = ColumnSelectUI()
        self.AnalROCUI = AnalROCUI()
        
        # menu action
        self.actionProteomics.triggered.connect(self.LoadProteinFile)
        self.actionDatabase.triggered.connect(self.LoadProteinComplex)
        self.actionCalcROC.triggered.connect(self.OpenAnalROC)
        
        # button action
        self.ButtonGroup1.clicked.connect(self.SetProteinTable1)
        self.ButtonGroup2.clicked.connect(self.SetProteinTable2)
        self.ButtonClearFileList.clicked.connect(self.ClearProteinFile)
        self.ButtonDatabaseConfirm.clicked.connect(self.SetProteinComplex)
        self.ButtonDatabaseRemove.clicked.connect(self.RemoveProteinComplex)
        self.ButtonClearDatabase.clicked.connect(self.ClearProteinComplex)
        self.ButtonCalcComplex.clicked.connect(self.CalcProteinComplexChange)
        self.ButtonShowCurve.clicked.connect(self.PlotProteinComplex)
        self.ButtonSaveComp.clicked.connect(self.SaveProteinComplex)
        
        self.ColumnSelectUI.ButtonColumnSelect.clicked.connect(self.SetProteinColumn)
        self.ColumnSelectUI.ButtonColumnCancel.clicked.connect(self.ColumnSelectUI.close)
        
        self.AnalROCUI.pushButtonDatabase.clicked.connect(self.LoadProteinPair)
        self.AnalROCUI.pushButtonConfirm.clicked.connect(self.ShowAnalROC)
        self.AnalROCUI.pushButtonCancel.clicked.connect(self.AnalROCUI.close)
        self.AnalROCUI.pushButtonPval.clicked.connect(self.CalcProteinPairChange)
        self.AnalROCUI.pushButtonCurve.clicked.connect(self.PlotProteinPairCurve)
        
        # table sort
        self.tableProteinComplex.setSortingEnabled(True)
        
        # server data
        self.columns = None
        self.prots = None
        self.resultDataComplex = []
        self.resultDataROC = []
        self.resultProtPair = []
        self.ROCNegValues = []  
    
    
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
        
    
    def LoadProteinComplex(self):
        options = QtWidgets.QFileDialog.Options()
        options |= QtWidgets.QFileDialog.DontUseNativeDialog
        fileName, _ = QtWidgets.QFileDialog.getOpenFileName(self,"Load", "","All Files (*);;CSV Files (*.csv)", options=options)
        if fileName:
            if fileName.split('.')[1] in ['csv', 'xlsx']:
                self.ListDatabase.addItem(fileName)
            else:
                self.ErrorMsg("Invalid format")
        else:
            pass
        
    
    def RemoveProteinComplex(self):
        self.ListDatabase.takeItem(self.ListDatabase.currentRow())
    
    
    def ClearProteinComplex(self):
        self.tableProteinComplex.setRowCount(0)
        self.ListDatabase.clear()
    
    
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


    def SetProteinTable1(self):
        data = self.SelectProteinTable()
        if data is None:
            return None
        self.ColumnSelectUI.listWidget.clear()
        all_cols = data.columns
        for c in all_cols:
            self.ColumnSelectUI.listWidget.addItem(c)
        self.ColumnSelectUI.show()


    def SetProteinTable2(self):
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
                self.tableProtein2.setModel(TableModel(data))
            except:
                self.ErrorMsg('No columns matched with Group 1')


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
        self.tableProtein1.setModel(TableModel(data))


    def FillProteinComplex(self, proteinComplex):
        self.tableProteinComplex.setRowCount(proteinComplex.shape[0])
        self.tableProteinComplex.setColumnCount(proteinComplex.shape[1])
        self.tableProteinComplex.setHorizontalHeaderLabels(proteinComplex.columns)
        self.tableProteinComplex.setVerticalHeaderLabels(proteinComplex.index.astype(str))
        for i in range(proteinComplex.shape[0]):
            for j in range(proteinComplex.shape[1]):
                if type(proteinComplex.iloc[i,j]) == np.float64:
                    item = QtWidgets.QTableWidgetItem()
                    item.setData(Qt.EditRole, QVariant(float(proteinComplex.iloc[i,j])))
                else:
                    item = QtWidgets.QTableWidgetItem(str(proteinComplex.iloc[i,j]))
                self.tableProteinComplex.setItem(i, j, item)        
     
     
    def TakeProteinComplex(self):
        ncol = self.tableProteinComplex.columnCount()
        nrow = self.tableProteinComplex.rowCount()
        header = [self.tableProteinComplex.horizontalHeaderItem(i).text() for i in range(ncol)]
        output = pd.DataFrame(np.zeros((nrow, ncol)))
        output.columns = header
        for i in range(nrow):
            for j in range(ncol):        
                v = self.tableProteinComplex.item(i, j).text()
                try:
                    v = float(v)
                except:
                    pass
                output.iloc[i,j] = v
        return output
     

    def SetProteinComplex(self):
        selectItem = self.ListDatabase.currentItem()
        try:
            if selectItem.text().split('.')[1] == 'csv':
                selectData = pd.read_csv(selectItem.text())
            elif selectItem.text().split('.')[1] == 'xlsx':
                selectData = pd.read_excel(selectItem.text())
            else:
                pass
        
            if 'Subunits_UniProt_IDs' not in selectData.columns:
                self.ErrorMsg("Subunits_UniProt_IDs' not in columns")
            else:
                self.FillProteinComplex(selectData)
        except:
            self.ErrorMsg("Cannot load the selected file")
        
    
    def CalcProteinComplexChange(self):
        self.resultDataComplex = []
        self.progressBar.setValue(0)
        self.ButtonCalcComplex.setEnabled(False)
        
        columns = self.columns
        proteinComplex = self.TakeProteinComplex()
        if len(proteinComplex) == 0:
            self.ErrorMsg('Not set protein complex')
        elif (self.tableProtein1.model() is None) or (self.tableProtein2.model() is None):
            self.ErrorMsg('Not set proteomics data')
        else:       
            proteinData1 = self.tableProtein1.model()._data
            proteinData2 = self.tableProtein2.model()._data

            data1 = proteinData1.loc[:, columns]
            data2 = proteinData2.loc[:, columns]
            dist1 = metrics.pairwise_distances(data1, metric = 'cityblock')
            dist2 = metrics.pairwise_distances(data2, metric = 'cityblock')
            prot1 = proteinData1.loc[:, 'Accession']
            prot2 = proteinData2.loc[:, 'Accession']

            self.ComplexThread = ComplexThread(prot1, dist1, prot2, dist2, proteinComplex)
            self.ComplexThread._ind.connect(self.ProcessBarComplex)
            self.ComplexThread._res.connect(self.ResultDataComplex)
            self.ComplexThread.start()
            self.ComplexThread.finished.connect(self.VisualizeComplex)


    def ResultDataComplex(self, msg):
        self.resultDataComplex.append(msg)
        

    def ProcessBarComplex(self, msg):
        self.progressBar.setValue(int(msg))
    
    
    def VisualizeComplex(self):
        proteinComplex = self.TakeProteinComplex()
        keep = [c for c in proteinComplex.columns 
                    if c not in ['Num subunit found', 'p-value (change)', 'Avg distance (change)', 'TPCA Sig 1', 'Avg distance 1', 'TPCA Sig 2', 'Avg distance 2']]
        proteinComplex = proteinComplex.loc[:,keep]
        
        resultDataComplex = pd.DataFrame(self.resultDataComplex)
        resultDataComplex.columns = ['Num subunit found', 'p-value (change)', 'Avg distance (change)', 'TPCA Sig 1', 'Avg distance 1', 'TPCA Sig 2', 'Avg distance 2']
        proteinComplex = pd.concat([proteinComplex, resultDataComplex], axis=1)
        proteinComplex = proteinComplex.sort_values(by = 'p-value (change)')
        proteinComplex = proteinComplex.reset_index(drop = True)
        self.FillProteinComplex(proteinComplex)
        self.ButtonCalcComplex.setEnabled(True)
        
    
    def PlotProteinComplex(self):
        colNames = self.columns
        header = [self.tableProteinComplex.horizontalHeaderItem(i).text() for i in range(self.tableProteinComplex.columnCount())]
        # print(header)
        try:
            i = self.tableProteinComplex.selectedIndexes()[0].row()
            j = list(header).index('Subunits_UniProt_IDs')
            proteinSubunit = self.tableProteinComplex.item(i, j).text()
        except:
            self.ErrorMsg('Can not plot protein curves')
            return None
        # print(proteinSubunit)
        proteinData1 = self.tableProtein1.model()._data
        proteinData2 = self.tableProtein2.model()._data
        # print(proteinData)

        self.figureG1.ProteinComplexFigure(proteinSubunit, proteinData1, colNames)
        self.figureG2.ProteinComplexFigure(proteinSubunit, proteinData2, colNames)
    

    def SaveProteinComplex(self):
        options = QtWidgets.QFileDialog.Options()
        options |= QtWidgets.QFileDialog.DontUseNativeDialog
        fileName, _ = QtWidgets.QFileDialog.getSaveFileName(self,"Save", ".csv","CSV Files (*.csv)", options=options)
        if fileName:
            data = self.TakeProteinComplex()
            data.to_csv(fileName)


    def LoadProteinPair(self):
        options = QtWidgets.QFileDialog.Options()
        options |= QtWidgets.QFileDialog.DontUseNativeDialog
        fileName, _ = QtWidgets.QFileDialog.getOpenFileName(self,"Load", "","All Files (*);;CSV Files (*.csv)", options=options)
        if fileName:
            if fileName.split('.')[1] == 'csv':
                proteinPair = pd.read_csv(fileName)
                header = [c for c in proteinPair.columns if c in ['Protein A', 'Protein B', 'Publications']]
                proteinPair = proteinPair.loc[:, header]
                self.AnalROCUI.tableView.setModel(TableModel(proteinPair))
            elif fileName.split('.')[1] == 'xlsx':
                proteinPair = pd.read_excel(fileName)
                header = [c for c in proteinPair.columns if c in ['Protein A', 'Protein B', 'Publications']]
                proteinPair = proteinPair.loc[:, header]
                self.AnalROCUI.tableView.setModel(TableModel(proteinPair))
            else:
                self.ErrorMsg("Invalid format")
            self.AnalROCUI.spinBoxRandom.setProperty("value", len(proteinPair))
        
    
    def OpenAnalROC(self):
        self.AnalROCUI.show()
        self.AnalROCUI.progressBar.setValue(0)
        self.AnalROCUI.comboBoxDataset.clear()
        self.AnalROCUI.comboBoxDataset.addItems(['Group1', 'Group2'])
        self.AnalROCUI.comboBoxDistance.clear()
        self.AnalROCUI.comboBoxDistance.addItems(['manhattan', 'cityblock', 'cosine', 'euclidean', 'l1', 'l2'])

        if self.tableProtein1.model() is None or (self.tableProtein2.model() is None):
            self.ErrorMsg('Please set proteomics data')
            self.AnalROCUI.close()
        else:
            pass
        
    
    def ShowAnalROC(self):
        pub_thres = self.AnalROCUI.spinBoxPub.value()
        columns = self.columns
        
        self.resultDataROC = []
        if self.tableProtein1.model() is None or (self.tableProtein2.model() is None):
            self.ErrorMsg("Protein matrix is not available")
            self.AnalROCUI.close()
        else:
            proteinData1 = self.tableProtein1.model()._data
            proteinData2 = self.tableProtein2.model()._data
                
        if self.AnalROCUI.tableView.model() is None:
            self.ErrorMsg("Protein pairs is not available")
            self.AnalROCUI.close()
        
        else:
            if self.AnalROCUI.comboBoxDataset.currentText() == 'Group1':
                proteinData = proteinData1
            else:
                proteinData = proteinData2
            
            proteinPair = self.AnalROCUI.tableView.model()._data
            header = [c for c in proteinPair.columns if c in ['Protein A', 'Protein B', 'Publications']]
            proteinPair = proteinPair.loc[:, header]
            proteinPair = proteinPair.reset_index(drop = True)
                
            if ('Protein A' not in proteinPair.columns) or ('Protein B' not in proteinPair.columns):
                self.ErrorMsg("Protein pairs is not available")
            
            if 'Publications' in proteinPair.columns:
                proteinPair = proteinPair[proteinPair['Publications'] >= pub_thres]
            
            prot = proteinData.loc[:, 'Accession']
            data = proteinData.loc[:, columns]
            dist = metrics.pairwise_distances(data, metric = self.AnalROCUI.comboBoxDistance.currentText())
            neg_values = np.triu(dist, k = 0).flatten()
            neg_values = neg_values[neg_values > 0]
            
            self.ROCNegValues = np.array(np.random.choice(list(neg_values), self.AnalROCUI.spinBoxRandom.value(), replace=False))
            self.ROCThread = ROCThread(prot, data, dist, proteinPair)
            self.ROCThread._ind.connect(self.ProcessBarROC)
            self.ROCThread._res.connect(self.ResultDataROC)
            self.ROCThread.start()
            self.ROCThread.finished.connect(self.VisualizeROC)

    
    def VisualizeROC(self):
        pos_values = 1 - np.array(self.resultDataROC)
        neg_values = 1 - np.array(self.ROCNegValues)
        values = list(pos_values) + list(neg_values)
        labels = [1] * len(pos_values) + [0] * len(neg_values)
        
        fpr, tpr, threshold = metrics.roc_curve(labels, values, pos_label = 1)
        auroc = np.round(metrics.roc_auc_score(labels, values), 4)
        self.AnalROCUI.figureROC.ROCFigure(fpr, tpr, auroc)
        

    def PlotProteinPairCurve(self):
        columns = self.columns
        proteinPair = self.AnalROCUI.tableView.model()._data
        proteinPair = proteinPair.reset_index(drop = True)
        
        i = self.AnalROCUI.tableView.selectedIndexes()[0].row()
        p1 = proteinPair.loc[i, 'Protein A']
        p2 = proteinPair.loc[i, 'Protein B']
        
        proteinData1 = self.tableProtein1.model()._data
        proteinData2 = self.tableProtein2.model()._data
        
        self.AnalROCUI.figureGroup1.ProteinPairFigure(p1, p2, proteinData1, columns)
        self.AnalROCUI.figureGroup2.ProteinPairFigure(p1, p2, proteinData2, columns)


    def ResultDataROC(self, msg):
        self.resultDataROC.append(msg)
        

    def ProcessBarROC(self, msg):
        self.AnalROCUI.progressBar.setValue(int(msg))
        
        
    def CalcProteinPairChange(self):
        self.resultProtPair = []
        self.AnalROCUI.progressBar.setValue(0)
        self.AnalROCUI.pushButtonPval.setEnabled(False)
        columns = self.columns
        pub_thres = self.AnalROCUI.spinBoxPub.value()
        
        if self.AnalROCUI.tableView.model() is None:
            self.ErrorMsg('Not set protein complex')
        elif (self.tableProtein1.model() is None) or (self.tableProtein2.model() is None):
            self.ErrorMsg('Not set proteomics data')
        else:      
            proteinData1 = self.tableProtein1.model()._data
            proteinData2 = self.tableProtein2.model()._data    
            proteinPair = self.AnalROCUI.tableView.model()._data
            
            header = [c for c in proteinPair.columns if c in ['Protein A', 'Protein B', 'Publications']]
            proteinPair = proteinPair.loc[:, header]
            proteinPair = proteinPair.reset_index(drop = True)
        
            if 'Publications' in proteinPair.columns:
                proteinPair = proteinPair[proteinPair['Publications'] >= pub_thres]
            data1 = proteinData1.loc[:, columns]
            data2 = proteinData2.loc[:, columns]
            dist1 = metrics.pairwise_distances(data1, metric = self.AnalROCUI.comboBoxDistance.currentText())
            dist2 = metrics.pairwise_distances(data2, metric = self.AnalROCUI.comboBoxDistance.currentText())
            prot1 = proteinData1.loc[:, 'Accession']
            prot2 = proteinData2.loc[:, 'Accession']
        
            n = self.AnalROCUI.spinBoxRandom.value()
            self.PairThread = PairThread(prot1, dist1, prot2, dist2, proteinPair, n)
            self.PairThread._ind.connect(self.ProcessBarROC)
            self.PairThread._res.connect(self.ResultProtPair)
            self.PairThread.start()
            self.PairThread.finished.connect(self.VisualizeProtPair)        


    def VisualizeProtPair(self):
        pub_thres = self.AnalROCUI.spinBoxPub.value()
        proteinPair = self.AnalROCUI.tableView.model()._data
        header = [c for c in proteinPair.columns if c in ['Protein A', 'Protein B', 'Publications']]
        proteinPair = proteinPair.loc[:, header]
        proteinPair = proteinPair.reset_index(drop = True)
        
        if 'Publications' in proteinPair.columns:
            proteinPair = proteinPair[proteinPair['Publications'] >= pub_thres]
        proteinPairDist = pd.DataFrame(self.resultProtPair)
        proteinPairDist.columns = ['Distance change', 'p-value', 'Distance Group1', 'Distance Group2']
        proteinPair = pd.concat([proteinPair, proteinPairDist], axis=1)
        proteinPair = proteinPair.sort_values(by = 'p-value')
        self.AnalROCUI.tableView.setModel(TableModel(proteinPair))
        self.AnalROCUI.pushButtonPval.setEnabled(True)
        

    def ResultProtPair(self, msg):
        self.resultProtPair.append(msg)
    

if __name__ == '__main__':
    import sys
    
    app = QApplication(sys.argv)
    ui = AnalTPCAUI()
    ui.show()
    sys.exit(app.exec_())