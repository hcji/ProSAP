# -*- coding: utf-8 -*-
"""
Created on Thu Apr 15 15:57:14 2021

@author: hcji
"""

from PyQt5 import QtCore, QtGui, QtWidgets
from PyQt5.QtWidgets import QApplication, QMainWindow

from AnalPvalComplex import Ui_Form

class AnalPvalComplexUI(QtWidgets.QWidget, Ui_Form):
    
    def __init__(self, parent=None): 
        super(AnalPvalComplexUI, self).__init__(parent)
        self.setupUi(self)


        

if __name__ == '__main__':
    import sys
    
    app = QApplication(sys.argv)
    ui = AnalPvalComplexUI()
    ui.show()
    sys.exit(app.exec_())