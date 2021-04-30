# -*- coding: utf-8 -*-
"""
Created on Wed Apr 14 16:00:17 2021

@author: hcji
"""

from PyQt5 import QtCore, QtGui, QtWidgets
from PyQt5.QtWidgets import QApplication, QMainWindow

from ColumnSelect import Ui_Form

class ColumnSelectUI(QMainWindow, Ui_Form):
    
    def __init__(self, parent=None):
        super(ColumnSelectUI, self).__init__(parent)
        self.setupUi(self)
        self.clearList()
        self.setWindowIcon(QtGui.QIcon("img/TPCA.ico"))
        
    def clearList(self):
        self.listWidget.clear()
        

if __name__ == '__main__':
    import sys
    
    app = QApplication(sys.argv)
    ui = ColumnSelectUI()
    ui.show()
    sys.exit(app.exec_())