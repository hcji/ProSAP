# -*- coding: utf-8 -*-
"""
Created on Thu Apr 15 15:57:14 2021

@author: hcji
"""

from PyQt5 import QtCore, QtGui, QtWidgets
from PyQt5.QtWidgets import QApplication, QMainWindow

from AnalROC import Ui_Form

class AnalROCUI(QtWidgets.QWidget, Ui_Form):
    
    def __init__(self, parent=None): 
        super(AnalROCUI, self).__init__(parent)
        self.setupUi(self)


        

if __name__ == '__main__':
    import sys
    
    app = QApplication(sys.argv)
    ui = AnalROCUI()
    ui.show()
    sys.exit(app.exec_())