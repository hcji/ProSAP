# -*- coding: utf-8 -*-
"""
Created on Tue May 11 10:00:45 2021

@author: hcji
"""


from PyQt5 import QtCore, QtGui, QtWidgets
from PyQt5.QtWidgets import QApplication, QMainWindow, QGridLayout, QHBoxLayout

from NPTSA import Ui_Form
from MakeFigure import MakeFigure
from matplotlib.backends.backend_qt5 import NavigationToolbar2QT as NavigationToolbar

class NPTSAUI(QtWidgets.QWidget, Ui_Form):
    
    def __init__(self, parent=None): 
        super(NPTSAUI, self).__init__(parent)
        self.setupUi(self)
        self.setWindowTitle("NPTSA Analysis")
        self.setWindowIcon(QtGui.QIcon("img/TPCA.ico"))

        self.figureTSA = MakeFigure(10, 10, dpi=250)
        self.figureTSA_ntb = NavigationToolbar(self.figureTSA, self)
        self.gridlayoutTSA = QGridLayout(self.groupBox)
        self.gridlayoutTSA.addWidget(self.figureTSA)
        self.gridlayoutTSA.addWidget(self.figureTSA_ntb)
        
        self.figureAvg = MakeFigure(10, 10, dpi=150)
        self.gridlayoutAvg = QGridLayout(self.groupBoxAvg)
        self.gridlayoutAvg.addWidget(self.figureAvg)
        
        self.tableWidgetProteinList.setSortingEnabled(True)
        self.comboBox.addItems(['fitness', 'euclidean', 'squaredeuclidean', 'seuclidean', 'mahalanobis', 'cityblock', 'chebychev', 'cosine'])


if __name__ == '__main__':
    import sys
    
    app = QApplication(sys.argv)
    ui = NPTSAUI()
    ui.show()
    sys.exit(app.exec_())