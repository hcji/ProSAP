# -*- coding: utf-8 -*-
"""
Created on Wed Jun 23 14:00:41 2021

@author: hcji
"""


import numpy as np
import pandas as pd

from seaborn import boxplot

from PyQt5 import QtCore, QtGui, QtWidgets
from PyQt5.QtWidgets import QApplication, QMainWindow, QGridLayout
from matplotlib.backends.backend_qt5 import NavigationToolbar2QT as NavigationToolbar

from AnalQC import Ui_MainWindow
from Utils import TableModel
from MakeFigure import MakeFigure


def calcRSD(data, Missing=False):
    prot = data.iloc[:,0]
    vals = data.iloc[:,1:]
    for c in vals.columns:
        vals[c] = pd.to_numeric(vals[c], errors='coerce')
    
    missingNum = np.sum(vals.isnull().values, axis = 1)
    outputRSD = np.repeat(np.nan, len(missingNum))
    
    if Missing:
        keep = np.arange(len(missingNum))
    else:
        keep = np.where(missingNum == 0)[0]
        vals = vals.iloc[keep,:]
    
    sd = np.nanstd(vals, axis = 1)
    m = np.nanmean(vals, axis = 1)
    rsd = sd / m
    outputRSD[keep] = rsd
    
    return prot, missingNum, outputRSD


class AnalQCUI(QtWidgets.QMainWindow, Ui_MainWindow):
    
    def __init__(self, parent=None): 
        super(AnalQCUI, self).__init__(parent)
        self.setupUi(self)
        self.setWindowTitle("QC Analysis")
        self.setWindowIcon(QtGui.QIcon("img/TPCA.ico"))
          
        self.allData = []
        self.fileNames = None
        self.fileLabel = None
        
        self.figureRSD = MakeFigure(10, 10, dpi = 250)
        self.figureRSD_ntb = NavigationToolbar(self.figureRSD, self)
        self.gridlayoutRSD = QGridLayout(self.groupBox1)
        self.gridlayoutRSD.addWidget(self.figureRSD)
        self.gridlayoutRSD.addWidget(self.figureRSD_ntb)
        
        self.figureBox = MakeFigure(10, 10, dpi = 250)
        self.figureBox_ntb = NavigationToolbar(self.figureBox, self)
        self.gridlayoutBox = QGridLayout(self.groupBox2)
        self.gridlayoutBox.addWidget(self.figureBox)
        self.gridlayoutBox.addWidget(self.figureBox_ntb)
        
        self.ButtonOpen.clicked.connect(self.LoadProteinFile)
        self.ButtonConfirm.clicked.connect(self.ShowWithMissing)
        self.ButtonExclude.clicked.connect(self.ShowWithoutMissing)
        self.ButtonShow.clicked.connect(self.ShowProtBoxPlot)
        
        
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
            self.tableWidgetFile.setRowCount(len(fileNames))
            self.tableWidgetFile.setColumnCount(2)
            self.tableWidgetFile.setHorizontalHeaderLabels(['FilePath', 'Label'])
            self.tableWidgetFile.setVerticalHeaderLabels(np.arange(len(fileNames)).astype(str))
            for i in range(len(fileNames)):
                item = QtWidgets.QTableWidgetItem(str(fileNames[i]))
                self.tableWidgetFile.setItem(i, 0, item)
        self.LoadData()

    
    def LoadData(self):
        self.allData = []
        fileNames = [self.tableWidgetFile.item(i,0).text() for i in range(self.tableWidgetFile.rowCount())]
        # print(fileNames)
        for f in fileNames:
            if f.split('.')[1] == 'csv':
                data = pd.read_csv(f)
                # allData.append(data)
            elif f.text().split('.')[1] == 'xlsx':
                data = pd.read_excel(f)
            else:
                self.ErrorMsg('Cannot be load the selected file')
                return None
            
            if 'Accession' not in data.columns:
                self.ErrorMsg('Accession is not given in the data')
                return None
            else:
                self.allData.append(data)
    
    
    def CalcRSD(self, keepMissing):
        allData = self.allData
        if len(allData) == 0:
            return None
        '''
        fileNames = ['D:/project/Wuqiong/drugs/20210615_080148_K562_autoSISPROT_CETSA_36samples combined_noFAIMS_directDIA_20210613_Report.csv',
                     'D:/project/Wuqiong/drugs/20210615_081722_K562_autoSISPROT_CETSA_36samples combined_noFAIMS_DIA_20210611_Report.csv']
        fileLabel = ['Sample1', 'Sample2']
        '''
        fileNames = [self.tableWidgetFile.item(i,0) for i in range(self.tableWidgetFile.rowCount())]
        fileLabel = [self.tableWidgetFile.item(i,1) for i in range(self.tableWidgetFile.rowCount())]
        if None in fileLabel:
            fileLabel = ['file_{}'.format(i) for i in range(len(fileNames))]
        
        self.fileNames = fileNames
        self.fileLabel = fileLabel
        resTable, allRSD = None, []
        for i in range(len(fileLabel)):
            Accession, nMissing, rsd = calcRSD(allData[i], keepMissing)
            if i == 0:
                resTable = pd.DataFrame(zip(Accession, nMissing))
                resTable.columns = ['Accession', 'nMissing_{}'.format(fileLabel[i])]
            else:
                res = pd.DataFrame(zip(Accession, nMissing))
                res.columns = ['Accession', 'nMissing_{}'.format(fileLabel[i])]
                resTable = resTable.merge(res, how='outer', on = 'Accession')
            allRSD.append(rsd)
        return resTable, allRSD
    
    
    def ShowWithMissing(self):
        resTable, allRSD = self.CalcRSD(True)        
        fileLabel = self.fileLabel
        
        data_lab, data_rsd = [], []
        for i in range(len(allRSD)):
            data_lab += [fileLabel[i]] * len(allRSD[i])
            data_rsd += list(allRSD[i])
        data_plt = pd.DataFrame({'FileLabel': data_lab, 'RSD': data_rsd})
        
        boxplot(x="FileLabel", y="RSD", data=data_plt)
        print(data_plt)
        
        self.tableView.setModel(TableModel(resTable))
        self.figureRSD.PlotQCRSD(data_plt)
        

    def ShowWithoutMissing(self):
        resTable, allRSD = self.CalcRSD(False)
        fileLabel = self.fileLabel
        
        data_lab, data_rsd = [], []
        for i in range(len(allRSD)):
            data_lab += [fileLabel[i]] * len(allRSD[i])
            data_rsd += list(allRSD[i])
        data_plt = pd.DataFrame({'FileLabel': data_lab, 'RSD': data_rsd})
        
        print(data_plt)
        boxplot(x="FileLabel", y="RSD", data=data_plt)
        
        self.tableView.setModel(TableModel(resTable))
        self.figureRSD.PlotQCRSD(data_plt)
    
        
    def ShowProtBoxPlot(self):
        fileLabel = self.fileLabel
        header = [self.tableView.horizontalHeaderItem(i).text() for i in range(self.tableView.columnCount())]
        i = self.tableView.selectedIndexes()[0].row()
        j = list(header).index('Accession')
        p = self.tableView.item(i, j).text()
        
        allData = self.allData
        data_lab, data_val = [], []
        for i, data in enumerate(allData):
            prot = data.iloc[:,0]
            vals = data.iloc[:,1:]
            for c in vals.columns:
                vals[c] = pd.to_numeric(vals[c], errors='coerce')
            wh = np.where(prot == p)[0]
            data_val.append(vals.iloc[wh,:].values[0,:])
            data_lab.append([fileLabel[i]] * len(vals.iloc[wh,:].values[0,:]))
        data_plt = pd.DataFrame(zip(data_lab, data_val))
        data_plt.columns = ['FileLabel', 'Values']
        self.figureBox.PlotQCBox(data_plt)
    
        
if __name__ == '__main__':
    import sys
    
    app = QApplication(sys.argv)
    ui = AnalQCUI()
    ui.show()
    sys.exit(app.exec_())
    