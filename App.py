# -*- coding: utf-8 -*-
"""
Created on Wed Apr 14 09:29:07 2021

@author: hcji
"""


import numpy as np
import pandas as pd

from scipy import stats
from sklearn import metrics

from PyQt5.QtCore import Qt
from PyQt5 import QtCore, QtGui, QtWidgets
from PyQt5.QtWidgets import QApplication, QMainWindow, QGridLayout
from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg as FigureCanvas
from matplotlib.backends.backend_qt5 import NavigationToolbar2QT as NavigationToolbar

from MainWindow import Ui_MainWindow
from ColumnSelectUI import ColumnSelectUI
from AnalROCUI import AnalROCUI
from AnalTSAUI import AnalTSAUI
from PreprocessUI import PreprocessUI
from Thread import CurveFitThread, ROCThread, PairThread, ComplexThread
from MakeFigure import MakeFigure
from Utils import TableModel, fit_curve


# proteinData1 = pd.read_csv('D:/project/PDFiles/MTX_CETSA.csv')
# proteinData2 = pd.read_csv('D:/project/PDFiles/PAN_CETSA.csv')
# columns = ['T37.0', 'T40.0', 'T43.0', 'T46.0', 'T49.0', 'T52.0', 'T55.0', 'T58.0', 'T61.0', 'T64.0']
# proteinPair = pd.read_csv('data/TPCA_TableS2_Protein_Pairs.csv')
# proteinComplex = pd.read_csv('data/TPCA_TableS3_Protein_Complex.csv')


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
        
        # Threads
        self.CurveFitThread = None
        self.ROCThread = None
        self.PairThread = None
        self.ComplexThread = None
        
        # widgets
        self.ColumnSelectUI = ColumnSelectUI()
        self.AnalROCUI = AnalROCUI()
        self.AnalTSAUI = AnalTSAUI()
        self.PreprocessUI = PreprocessUI()
        
        # menu action
        self.actionProteomics.triggered.connect(self.LoadProteinFile)
        self.actionDatabase.triggered.connect(self.LoadProteinComplex)
        self.actionPreprocessing.triggered.connect(self.OpenPreprocessing)
        self.action_CETSA.triggered.connect(self.OpenAnalTSA)
        self.actionCalcROC.triggered.connect(self.OpenAnalROC)
        
        # button action
        self.ButtonGroup1.clicked.connect(self.SetProteinTable1)
        self.ButtonGroup2.clicked.connect(self.SetProteinTable2)
        self.ButtonClearFileList.clicked.connect(self.ClearProteinFile)
        self.ButtonDatabaseConfirm.clicked.connect(self.SetProteinComplex)
        self.ButtonDatabaseRemove.clicked.connect(self.RemoveProteinComplex)
        self.ButtonClearDatabase.clicked.connect(self.ClearProteinComplex)
        self.ButtonShowCurve.clicked.connect(self.PlotProteinComplex)
        
        # server data
        self.columns = None
        self.prots = None
        self.proteinPair = None
        self.proteinComplex = None
        self.ProteinTable1 = None
        self.ProteinTable2 = None
        self.TSA_table = pd.DataFrame()
        self.resultDataComplex = []
        self.resultDataTSA = []
        self.resultDataROC = []
        self.resultProtPair = []
        self.ROCNegValues = []
        
        
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
        fileName, _ = QtWidgets.QFileDialog.getOpenFileName(self,"QFileDialog.getOpenFileName()", "","All Files (*);;CSV Files (*.csv)", options=options)
        if fileName:
            if fileName.split('.')[1] == 'csv':
                self.ListFile.addItem(fileName)
            else:
                self.ErrorMsg("Invalid format")
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
                self.ErrorMsg("Invalid format")
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
        self.ProteinTable1 = data
        self.ColumnSelectUI.listWidget.clear()
        all_cols = data.columns
        for c in all_cols:
            self.ColumnSelectUI.listWidget.addItem(c)
        self.ColumnSelectUI.show()
        self.ColumnSelectUI.ButtonColumnSelect.clicked.connect(self.SetProteinColumn)
        self.ColumnSelectUI.ButtonColumnCancel.clicked.connect(self.ColumnSelectUI.close)


    def SetProteinTable2(self):
        if self.columns == None:
            self.ErrorMsg('Please set Group 1')
        else:
            columns = ['Accession'] + self.columns
            data = self.SelectProteinTable()
            data = data.loc[:, columns]
            self.tableProtein2.setModel(TableModel(data))
            self.ProteinTable2 = data


    def SetProteinColumn(self):
        self.columns = [i.text() for i in self.ColumnSelectUI.listWidget.selectedItems()]
        columns = ['Accession'] + self.columns
        self.ColumnSelectUI.close()
        data = self.ProteinTable1.loc[:, columns]
        self.tableProtein1.setModel(TableModel(data))
        

    def SetProteinComplex(self):
        selectItem = self.ListDatabase.currentItem()
        selectData = pd.read_csv(selectItem.text())
        if 'Subunits_UniProt_IDs' not in selectData.columns:
            self.ErrorMsg("Subunits_UniProt_IDs' not in columns")
        else:
            self.tableProteinComplex.setRowCount(selectData.shape[0])
            self.tableProteinComplex.setColumnCount(selectData.shape[1])
            self.tableProteinComplex.setHorizontalHeaderLabels(selectData.columns)
            self.tableProteinComplex.setVerticalHeaderLabels(selectData.index.astype(str))
            for i in range(selectData.shape[0]):
                for j in range(selectData.shape[1]):
                    item = QtWidgets.QTableWidgetItem(str(selectData.iloc[i,j]))
                    self.tableProteinComplex.setItem(i, j, item)
            self.proteinComplex = selectData
        self.ButtonCalcComplex.clicked.connect(self.CalcProteinComplexChange)
        
    
    def CalcProteinComplexChange(self):
        columns = self.columns
        proteinComplex = self.proteinComplex
        proteinData1 = self.tableProtein1.model()._data
        proteinData2 = self.tableProtein2.model()._data
        
        data1 = proteinData1.loc[:, columns]
        data2 = proteinData2.loc[:, columns]
        dist1 = metrics.pairwise_distances(data1, metric = 'manhattan')
        dist2 = metrics.pairwise_distances(data2, metric = 'manhattan')
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
        proteinComplex = self.proteinComplex
        resultDataComplex = pd.DataFrame(self.resultDataComplex)
        resultDataComplex.columns = ['Num subunit found', 'p-value (change)', 'Avg distance (change)', 'TPCA Sig 1', 'Avg distance 1', 'TPCA Sig 2', 'Avg distance 2']
        proteinComplex = pd.concat([proteinComplex, resultDataComplex], axis=1)
        proteinComplex = proteinComplex.sort_values(by = 'p-value (change)')
        
        self.tableProteinComplex.setRowCount(proteinComplex.shape[0])
        self.tableProteinComplex.setColumnCount(proteinComplex.shape[1])
        self.tableProteinComplex.setHorizontalHeaderLabels(proteinComplex.columns)
        self.tableProteinComplex.setVerticalHeaderLabels(proteinComplex.index.astype(str))
        for i in range(proteinComplex.shape[0]):
            for j in range(proteinComplex.shape[1]):
                item = QtWidgets.QTableWidgetItem(str(proteinComplex.iloc[i,j]))
                self.tableProteinComplex.setItem(i, j, item)       
    
    
    def PlotProteinComplex(self):       
        colNames = self.columns        
        header = [self.tableProteinComplex.horizontalHeaderItem(i).text() for i in range(self.tableProteinComplex.columnCount())]
        # print(header)
        i = self.tableProteinComplex.selectedItems()[0].row()
        j = header.index('Subunits_UniProt_IDs')
        proteinSubunit = self.tableProteinComplex.item(i, j).text()
        # print(proteinSubunit)
        proteinData1 = self.tableProtein1.model()._data
        proteinData2 = self.tableProtein2.model()._data
        # print(proteinData)
        F_gr1 = MakeFigure(1.5, 1.5)
        F_gr2 = MakeFigure(1.5, 1.5)
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
        

    def LoadProteinPair(self):
        options = QtWidgets.QFileDialog.Options()
        options |= QtWidgets.QFileDialog.DontUseNativeDialog
        fileName, _ = QtWidgets.QFileDialog.getOpenFileName(self,"QFileDialog.getOpenFileName()", "","All Files (*);;CSV Files (*.csv)", options=options)
        if fileName:
            if fileName.split('.')[1] == 'csv':
                self.proteinPair = pd.read_csv(fileName)
                self.AnalROCUI.tableView.setModel(TableModel(self.proteinPair))
            else:
                self.ErrorMsg("Invalid format")
        else:
            self.AnalROCUI.spinBoxRandom.setProperty("value", len(self.proteinPair))
        
    
    def OpenAnalROC(self):
        self.AnalROCUI.show()
        self.AnalROCUI.progressBar.setValue(0)
        self.AnalROCUI.comboBoxDataset.clear()
        self.AnalROCUI.comboBoxDataset.addItems(['Group1', 'Group2'])
        self.AnalROCUI.comboBoxDistance.clear()
        self.AnalROCUI.comboBoxDistance.addItems(['manhattan', 'cityblock', 'cosine', 'euclidean', 'l1', 'l2'])

        if self.tableProtein1.model() is None or (self.tableProtein2.model() is None):
            pass
        else:
            self.AnalROCUI.pushButtonDatabase.clicked.connect(self.LoadProteinPair)
            self.AnalROCUI.pushButtonConfirm.clicked.connect(self.ShowAnalROC)
            self.AnalROCUI.pushButtonCancel.clicked.connect(self.AnalROCUI.close)
        
    
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
        
        F = MakeFigure(1.2, 1.2)
        F.axes.cla()
        F.ROCFigure(fpr, tpr, auroc)
        f = QtWidgets.QGraphicsScene()
        f.addWidget(F)
        self.AnalROCUI.graphicsView.setScene(f)
        self.AnalROCUI.pushButtonPval.clicked.connect(self.CalcProteinPairChange)


    def ResultDataROC(self, msg):
        self.resultDataROC.append(msg)
        

    def ProcessBarROC(self, msg):
        self.AnalROCUI.progressBar.setValue(int(msg))
        
        
    def CalcProteinPairChange(self):
        columns = self.columns
        pub_thres = self.AnalROCUI.spinBoxPub.value()
        proteinData1 = self.tableProtein1.model()._data
        proteinData2 = self.tableProtein2.model()._data    
        proteinPair = self.AnalROCUI.tableView.model()._data
        
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
        if 'Publications' in proteinPair.columns:
            proteinPair = proteinPair[proteinPair['Publications'] >= pub_thres]
        proteinPairDist = pd.DataFrame(self.resultProtPair)
        proteinPairDist.columns = ['Distance change', 'p-value', 'Distance Group1', 'Distance Group2']
        proteinPair = pd.concat([proteinPair, proteinPairDist], axis=1)
        proteinPair = proteinPair.sort_values(by = 'p-value')
        self.AnalROCUI.tableView.setModel(TableModel(proteinPair))
        

    def ResultProtPair(self, msg):
        self.resultProtPair.append(msg)

    
    def OpenPreprocessing(self):
        self.PreprocessUI.show()
        
    
    def OpenAnalTSA(self):
        self.AnalTSAUI.show()
        if self.tableProtein1.model() is None or (self.tableProtein2.model() is None):
            pass
        else:
            self.AnalTSAUI.ButtonConfirm.clicked.connect(self.ShowAnalTSA)
            self.AnalTSAUI.ButtonCancel.clicked.connect(self.AnalTSAUI.close)
    
    
    def ProcessBarTSA(self, msg):
        self.AnalTSAUI.progressBar.setValue(int(msg))
        
    
    def ResultDataTSA(self, msg):
        self.resultDataTSA.append(msg)
        # print(msg)
    
    
    def ShowAnalTSA(self):
        columns = self.columns
        self.ColumnSelectUI.close()
        self.resultDataTSA = []
        
        self.AnalTSAUI.tableWidgetProteinList.clear()
        self.AnalTSAUI.progressBar.setValue(0)
        
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
    
    
    def VisualizeTSA(self):
        proteinData1 = self.tableProtein1.model()._data
        proteinData2 = self.tableProtein2.model()._data
        
        prots = self.prots
        columns = self.columns
        
        res = pd.DataFrame(self.resultDataTSA)
        res.columns = ['Group1_R2', 'Group2_R2', 'Group1_Tm', 'Group2_Tm', 'delta_Tm']
    
        delta_Tm = res['delta_Tm']
        p_Val = []
        for i in range(len(res)):
            s = delta_Tm[i]
            pv = stats.t.sf((s - np.mean(delta_Tm)) / np.std(delta_Tm), len(delta_Tm)-1)
            p_Val.append(pv)
        score = -np.log10(np.array(p_Val)) * (res['Group1_R2'] * res['Group2_R2']) ** 2
    
        res['Accession'] = prots
        res['delta_Tm'] = delta_Tm
        res['p_Val (-log10)'] = -np.log10(p_Val)
        res['Score'] = score
        res = np.round(res, 3)
    
        res = res[['Accession', 'Score', 'p_Val (-log10)', 'delta_Tm', 'Group1_R2', 'Group2_R2', 'Group1_Tm', 'Group2_Tm']]
        TSA_table = res.sort_values(by = 'Score',axis = 0, ascending = False)
        
        self.TSA_table = TSA_table
        self.AnalTSAUI.tableWidgetProteinList.setRowCount(TSA_table.shape[0])
        self.AnalTSAUI.tableWidgetProteinList.setColumnCount(TSA_table.shape[1])
        self.AnalTSAUI.tableWidgetProteinList.setHorizontalHeaderLabels(TSA_table.columns)
        self.AnalTSAUI.tableWidgetProteinList.setVerticalHeaderLabels(TSA_table.index.astype(str))
        for i in range(TSA_table.shape[0]):
            for j in range(TSA_table.shape[1]):
                item = QtWidgets.QTableWidgetItem(str(TSA_table.iloc[i,j]))
                self.AnalTSAUI.tableWidgetProteinList.setItem(i, j, item)
        
        F = MakeFigure(1.2, 0.7)
        F.axes.cla()
        F.AverageTSAFigure(proteinData1, proteinData2, columns)
        f = QtWidgets.QGraphicsScene()
        f.addWidget(F)
        self.AnalTSAUI.graphicsViewAvgCurve.setScene(f)
        self.AnalTSAUI.ButtonShow.clicked.connect(self.ShowTSACurve)

    
    def ShowTSACurve(self):
        columns = self.columns
        proteinData1 = self.tableProtein1.model()._data
        proteinData2 = self.tableProtein2.model()._data
        
        header = [self.AnalTSAUI.tableWidgetProteinList.horizontalHeaderItem(i).text() for i in range(self.AnalTSAUI.tableWidgetProteinList.columnCount())]
        i = self.AnalTSAUI.tableWidgetProteinList.selectedItems()[0].row()
        j = header.index('Accession')
        ProteinAccession = self.AnalTSAUI.tableWidgetProteinList.item(i, j).text()

        F = MakeFigure(1.3, 1.3)
        F.axes.cla()
        F.SingleTSAFigure(proteinData1, proteinData2, columns, ProteinAccession)
        f = QtWidgets.QGraphicsScene()
        f.addWidget(F)
        self.AnalTSAUI.graphicsViewTSACurve.setScene(f)
        self.AnalTSAUI.pushButtonSave.clicked.connect(self.SaveTSAData)


    def SaveTSAData(self):
        options = QtWidgets.QFileDialog.Options()
        options |= QtWidgets.QFileDialog.DontUseNativeDialog
        fileName, _ = QtWidgets.QFileDialog.getSaveFileName(self,"QFileDialog.getOpenFileName()", "","All Files (*);;CSV Files (*.csv)", options=options)
        data = self.TSA_table
        data.to_csv(fileName)
    

if __name__ == '__main__':
    
    import sys
    app = QApplication(sys.argv)
    ui = TCPA_Main()
    ui.show()
    sys.exit(app.exec_())
    