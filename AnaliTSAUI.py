# -*- coding: utf-8 -*-
"""
Created on Fri May  7 08:36:21 2021

@author: jihon
"""

import numpy as np
import pandas as pd
from scipy.stats import ttest_ind

from PyQt5 import QtCore, QtGui, QtWidgets
from PyQt5.QtWidgets import QApplication, QMainWindow, QGridLayout, QHBoxLayout

from AnaliTSA import Ui_Form
from MakeFigure import MakeFigure
from matplotlib.backends.backend_qt5 import NavigationToolbar2QT as NavigationToolbar

from Utils import TableModel
from Thread import PreprocessThread
from MakeFigure import MakeFigure
from ColumnSelectUI import ColumnSelectUI


from rpy2 import robjects
from rpy2.robjects import numpy2ri, pandas2ri

numpy2ri.activate()
pandas2ri.activate()
robjects.r('''source('R/iTSA.R')''')

do_limma = robjects.globalenv['do_limma']
p_value_adjust = robjects.globalenv['p_value_adjust']


'''
proteinData = pd.read_csv('data/iTSA_TableS1_Proteomics.csv')
columns = ['V_log2.i._TMT_1_iTSA52', 'V_log2.i._TMT_3_iTSA52',
       'V_log2.i._TMT_5_iTSA52', 'V_log2.i._TMT_7_iTSA52',
       'V_log2.i._TMT_9_iTSA52', 'D_log2.i._TMT_2_iTSA52',
       'D_log2.i._TMT_4_iTSA52', 'D_log2.i._TMT_6_iTSA52',
       'D_log2.i._TMT_8_iTSA52', 'D_log2.i._TMT_10_iTSA52']
labels = [0,0,0,0,0,1,1,1,1,1]
'''

def do_ttest(X, y):
    i = np.where(y == 'Case')[0]
    j = np.where(y == 'Control')[0]
    pval = []
    for k in range(X.shape[0]):
        pval.append(ttest_ind(X.iloc[k, i], X.iloc[k, j]).pvalue)
    return np.expand_dims(np.array(pval), 1)


class iTSA:
    
    def __init__(self, method = 'Limma'):
        self.method = method
        
        
    def fit_data(self, X, y):
        self.X = X
        self.y = y
        
        if self.method == 'Limma':
            pass
        elif self.method == 'Ttest':
            self.res = do_ttest(X, y)
        else:
            raise IOError('{} is not a support method'.format(self.method))
        
        
    def fold_change(self):
        case_val = np.nanmean(self.X.loc[:, self.y == 'Case'], axis = 1)
        cont_val = np.nanmean(self.X.loc[:, self.y == 'Control'], axis = 1)
        return np.log2(case_val / cont_val)
        
        
    def p_values(self):
        return np.array(self.res)[:,0]
    
    
    def adjusted_p_values(self):
        return np.array(p_value_adjust(self.res))
    
    
    def summary(self):
        fc = self.fold_change()
        pval = self.p_values()
        apval = self.adjusted_p_values()
        return pd.DataFrame({'index': self.index, 'fold_change': fc, 'p_values': pval, 'adjusted_p_values': apval})



class AnaliTSAUI(QtWidgets.QWidget, Ui_Form):
    
    def __init__(self, parent=None): 
        super(AnaliTSAUI, self).__init__(parent)
        self.setupUi(self)
        self.setWindowTitle("iTSA Analysis")
        self.setWindowIcon(QtGui.QIcon("img/TPCA.ico"))

        self.figureTSA = MakeFigure(5, 5)
        self.figureTSA_ntb = NavigationToolbar(self.figureTSA, self)
        self.gridlayoutTSA = QGridLayout(self.groupBoxVolcano)
        self.gridlayoutTSA.addWidget(self.figureTSA)
        self.gridlayoutTSA.addWidget(self.figureTSA_ntb)
        self.ColumnSelectUI = ColumnSelectUI()
        self.comboBoxMerging.addItems(['comboBoxMethod'])
        
        self.pushButtonData.clicked.connect(self.LoadProteinFile)
        self.pushButtonOK.clicked.connect(self.DoPropress)
        
        self.data = None
        self.columns = None
        self.label = None
    
    
    def LoadProteinFile(self):
        options = QtWidgets.QFileDialog.Options()
        options |= QtWidgets.QFileDialog.DontUseNativeDialog
        fileName, _ = QtWidgets.QFileDialog.getOpenFileName(self,"QFileDialog.getOpenFileName()", "","All Files (*);;CSV Files (*.csv)", options=options)
        if fileName:
            if fileName.split('.')[1] == 'csv':
                self.data = pd.read_csv(fileName)
                self.tableViewData.setModel(TableModel(self.data))
                
                self.columns = self.data.columns
                for c in self.columns:
                    self.ColumnSelectUI.listWidget.addItem(c)
                self.ColumnSelectUI.show()
                self.ColumnSelectUI.ButtonColumnSelect.clicked.connect(self.SetLabel)
                self.ColumnSelectUI.ButtonColumnCancel.clicked.connect(self.ColumnSelectUI.close)
            else:
                self.ErrorMsg("Invalid format")
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
        if None in [self.tableWidgetTemp.item(i,1) for i in range(self.tableWidgetTemp.rowCount())]:
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
            lbs = np.unique(y)   
            i = np.where(y == lbs[0])[0]
            j = np.where(y == lbs[0])[1]
            self.pval = []
            for k in range(X.shape[0]):
                self.pval.append(ttest_ind(X.iloc[k, i], X.iloc[k, j]).pvalue)
        
        

if __name__ == '__main__':
    import sys
    
    app = QApplication(sys.argv)
    ui = AnaliTSAUI()
    ui.show()
    sys.exit(app.exec_())