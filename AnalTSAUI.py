# -*- coding: utf-8 -*-
"""
Created on Thu Apr 15 15:57:14 2021

@author: hcji
"""

from PyQt5 import QtCore, QtGui, QtWidgets
from PyQt5.QtWidgets import QApplication, QMainWindow

from AnalTSA import Ui_Form

class AnalTSAUI(QtWidgets.QWidget, Ui_Form):
    
    def __init__(self, parent=None): 
        super(AnalTSAUI, self).__init__(parent)
        self.setupUi(self)
        self.setWindowTitle("TPP Analysis")


        

if __name__ == '__main__':
    import sys
    
    app = QApplication(sys.argv)
    ui = AnalTSAUI()
    ui.show()
    sys.exit(app.exec_())