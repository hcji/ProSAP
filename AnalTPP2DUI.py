# -*- coding: utf-8 -*-
"""
Created on Mon Nov 29 10:29:37 2021

@author: jihon
"""

import numpy as np
import pandas as pd

from PyQt5 import QtCore, QtGui, QtWidgets
from PyQt5.QtCore import Qt, QVariant
from PyQt5.QtWidgets import QApplication, QMainWindow, QGridLayout, QHBoxLayout

from AnalTPP2D import Ui_Form
from MakeFigure import MakeFigure
from matplotlib.backends.backend_qt5 import NavigationToolbar2QT as NavigationToolbar

from Utils import TableModel
from TPP2D import TPP2D

from RunningUI import Running_Win


'''
data = pd.read_csv('data/TPP2D/small_example.csv')

'''


class AnalTPP2DUI(QtWidgets.QWidget, Ui_Form):
    
    def __init__(self, parent=None): 
        super(AnalTPP2DUI, self).__init__(parent)
        self.setupUi(self)
        self.setWindowTitle("TPP2D Analysis")
        self.setWindowIcon(QtGui.QIcon("img/TPCA.ico"))
        self.Running_Win = Running_Win()
        
        self.comboBoxMethod.addItems(['Kurzawa'])
        
        self.pushButtonData.clicked.connect(self.LoadProteinFile)
        self.pushButtonOK.clicked.connect(self.DoPropress)
        self.pushButtonClose.clicked.connect(self.close)
        self.pushButtonSave.clicked.connect(self.SaveData)
        self.pushButtonShow.clicked.connect(self.PlotProteinHeatmap)
        
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
        
        self.data = None
        self.maxit = None
        self.B = None
        self.result = None
    
    
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
        self.tableViewData.setModel(TableModel(self.data))
        
        
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
    
    
    def DoPropress(self):
        self.WarnMsg('The processing may cost quite a lone time, please wait!')
        self.Running_Win.show()
        tpp2d = TPP2D(self.data)
        check = tpp2d.check()
        if not check:
            self.ErrorMsg('Necessary columns missing in your data, please check!')
            self.Running_Win.close()
            return None
        
        self.maxit = self.spinBoxMaxIt.value()
        self.B = self.spinBoxBoost.value()
        self.alpha = self.doubleSpinBoxAlpha.value()
        
        result = tpp2d.fit_data(self.maxit, self.B)
        if result is None:
            self.ErrorMsg('There is something wrong with R program, please check!')
            self.Running_Win.close()
            return
        elif len(result) == 0:
            self.ErrorMsg('There is something wrong when processing your data, please check!')
            self.Running_Win.close()
            return
        else:
            self.result = result
            self.hits = tpp2d.find_hits(self.alpha)
            self.hits['plot x'] = np.sign(self.result['slopeH1']) * np.sqrt(self.result['rssH0'] - self.result['rssH1'])
            self.hits['plot y'] = np.log2(self.result['F_statistic'] + 1)
            self.FillTable(self.hits)
            self.Running_Win.close()
            self.figureVolcano.TPP2D_Volcano(self.result, self.hits)
    
    
    def PlotProteinHeatmap(self):     
        try:
            data = self.tableViewData.model()._data
        except:
            return None
        header = [self.tableWidgetProteinList.horizontalHeaderItem(i).text() for i in range(self.tableWidgetProteinList.columnCount())]
        i = self.tableWidgetProteinList.selectedItems()[0].row()
        j = header.index('clustername')
        ProteinAccession = self.tableWidgetProteinList.item(i, j).text()
        self.figureHeatmap.TPP2D_protHeatmap(data, ProteinAccession)
    

    def TakeProteinList(self):
        ncol = self.tableWidgetProteinList.columnCount()
        nrow = self.tableWidgetProteinList.rowCount()
        header = [self.tableWidgetProteinList.horizontalHeaderItem(i).text() for i in range(ncol)]
        output = pd.DataFrame(np.zeros((nrow, ncol)))
        output.columns = header
        for i in range(nrow):
            for j in range(ncol):        
                v = self.tableWidgetProteinList.item(i, j).text()
                try:
                    v = float(v)
                except:
                    pass
                output.iloc[i,j] = v
        return output

    
    def SaveData(self):
        options = QtWidgets.QFileDialog.Options()
        options |= QtWidgets.QFileDialog.DontUseNativeDialog
        fileName, _ = QtWidgets.QFileDialog.getSaveFileName(self, "Save", ".csv","All Files (*);;CSV Files (*.csv)", options=options)
        if fileName:
            data = data = self.TakeProteinList()
            data.to_csv(fileName, index = False)
    
    

if __name__ == '__main__':
    import sys
    
    app = QApplication(sys.argv)
    ui = AnalTPP2D()
    ui.show()
    sys.exit(app.exec_())