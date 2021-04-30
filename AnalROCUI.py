# -*- coding: utf-8 -*-
"""
Created on Thu Apr 15 15:57:14 2021

@author: hcji
"""

from PyQt5 import QtCore, QtGui, QtWidgets
from PyQt5.QtWidgets import QApplication, QMainWindow, QGridLayout, QMessageBox, QMenu

from AnalROC import Ui_Form
from MakeFigure import MakeFigure
from matplotlib.backends.backend_qt5 import NavigationToolbar2QT as NavigationToolbar


class AnalROCUI(QtWidgets.QWidget, Ui_Form):
    
    def __init__(self, parent=None): 
        super(AnalROCUI, self).__init__(parent)
        self.setupUi(self)
        self.setWindowTitle("Protein pairs analysis")
        
        self.figureGroup1 = MakeFigure(3, 3)
        self.figureGroup2 = MakeFigure(3, 3)
        self.figureGroup1_ntb = NavigationToolbar(self.figureGroup1, self)
        self.figureGroup2_ntb = NavigationToolbar(self.figureGroup2, self)
        
        self.gridlayoutBox1 = QGridLayout(self.groupBox_1)
        self.gridlayoutBox2 = QGridLayout(self.groupBox_2)
        
        self.gridlayoutBox1.addWidget(self.figureGroup1)
        self.gridlayoutBox1.addWidget(self.figureGroup1_ntb)
        self.gridlayoutBox2.addWidget(self.figureGroup2)
        self.gridlayoutBox2.addWidget(self.figureGroup2_ntb)
        

if __name__ == '__main__':
    import sys
    
    app = QApplication(sys.argv)
    ui = AnalROCUI()
    ui.show()
    sys.exit(app.exec_())