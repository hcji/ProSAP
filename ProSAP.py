# -*- coding: utf-8 -*-
"""
Created on Mon Jun 21 07:57:08 2021

@author: jihon
"""


from PyQt5 import QtCore, QtGui, QtWidgets
from PyQt5.QtWidgets import QApplication, QMainWindow

from ProSAP_ import Ui_Form
from AnalTSAUI import AnalTSAUI
from AnaliTSAUI import AnaliTSAUI
from PreprocessUI import PreprocessUI
from AnalTPCAUI import AnalTPCAUI


class ProSAPUI(QMainWindow, Ui_Form):
    
    def __init__(self, parent=None):
        super(ProSAPUI, self).__init__(parent)
        self.setupUi(self)
        self.setWindowTitle("ProSAP")
        self.setWindowIcon(QtGui.QIcon("img/ProSAP.ico"))
        
        self.AnalTPCAUI = AnalTPCAUI()
        self.AnalTSAUI = AnalTSAUI()
        self.AnaliTSAUI = AnaliTSAUI()
        self.PreprocessUI = PreprocessUI()
        
        self.pushButtonPreprocess.clicked.connect(self.PreprocessUI.show)
        self.pushButtonITSA.clicked.connect(self.AnaliTSAUI.show)
        self.pushButtonTPP.clicked.connect(self.AnalTSAUI.show)
        self.pushButtonTPCA.clicked.connect(self.AnalTPCAUI.show)      
        
        

if __name__ == '__main__':
    import sys
    
    app = QApplication(sys.argv)
    ui = ProSAPUI()
    ui.show()
    sys.exit(app.exec_())