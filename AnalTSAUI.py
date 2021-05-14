# -*- coding: utf-8 -*-
"""
Created on Thu Apr 15 15:57:14 2021

@author: hcji
"""

import numpy as np

from PyQt5 import QtCore, QtGui, QtWidgets
from PyQt5.QtCore import Qt, QVariant
from PyQt5.QtWidgets import QApplication, QMainWindow, QGridLayout, QHBoxLayout

from AnalTSA import Ui_Form
from MakeFigure import MakeFigure
from matplotlib.backends.backend_qt5 import NavigationToolbar2QT as NavigationToolbar

class AnalTSAUI(QtWidgets.QWidget, Ui_Form):
    
    def __init__(self, parent=None): 
        super(AnalTSAUI, self).__init__(parent)
        self.setupUi(self)
        self.setWindowTitle("TPP Analysis")
        self.setWindowIcon(QtGui.QIcon("img/TPCA.ico"))

        self.figureTSA = MakeFigure(10, 10, dpi = 250)
        self.figureTSA_ntb = NavigationToolbar(self.figureTSA, self)
        self.gridlayoutTSA = QGridLayout(self.groupBox)
        self.gridlayoutTSA.addWidget(self.figureTSA)
        self.gridlayoutTSA.addWidget(self.figureTSA_ntb)
        
        self.figureAvg = MakeFigure(10, 10, dpi = 150)
        self.gridlayoutAvg = QGridLayout(self.groupBoxAvg)
        self.gridlayoutAvg.addWidget(self.figureAvg)
        
        self.tableWidgetProteinList.setSortingEnabled(True)
        
        
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
                
        

if __name__ == '__main__':
    import sys
    
    app = QApplication(sys.argv)
    ui = AnalTSAUI()
    ui.show()
    sys.exit(app.exec_())