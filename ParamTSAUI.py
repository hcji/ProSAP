# -*- coding: utf-8 -*-
"""
Created on Tue Jun 22 15:51:43 2021

@author: hcji
"""


from PyQt5 import QtCore, QtGui, QtWidgets
from PyQt5.QtWidgets import QApplication, QMainWindow
from ParamTSA import Ui_Form


class ParamsTSAUI(QtWidgets.QWidget, Ui_Form):
    
    def __init__(self, parent=None):
        super(ParamsTSAUI, self).__init__(parent)
        self.setupUi(self)



if __name__ == '__main__':
    import sys
    
    app = QApplication(sys.argv)
    ui = ParamsTSAUI()
    ui.show()
    sys.exit(app.exec_())