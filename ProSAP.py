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
from AnalTPP2DUI import AnalTPP2DUI


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
        self.AnalTPP2DUI = AnalTPP2DUI()
        
        self.pushButtonPreprocess.clicked.connect(self.PreprocessUI.show)
        self.pushButtonITSA.clicked.connect(self.AnaliTSAUI.show)
        self.pushButtonTPP.clicked.connect(self.AnalTSAUI.show)
        self.pushButtonTPCA.clicked.connect(self.AnalTPCAUI.show)
        self.pushButton2DTPP.clicked.connect(self.AnalTPP2DUI.show)   
        
        

if __name__ == '__main__':
    import sys
    
    app = QApplication(sys.argv)
    ui = ProSAPUI()
    ui.show()
    sys.exit(app.exec_())