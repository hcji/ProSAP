# -*- coding: utf-8 -*-
"""
Created on Thu Jun 24 14:57:27 2021

@author: hcji
"""


import sys
from PyQt5 import QtCore, QtWidgets
from PyQt5.QtGui import QMovie
from PyQt5.QtWidgets import QApplication
from Running import Ui_Dialog
from Utils import get_pic
from memory_pic import loading_gif

class Running_Win(QtWidgets.QWidget, Ui_Dialog):
    
    def __init__(self, parent=None): 
        super(Running_Win, self).__init__(parent)
        self.setupUi(self)
        self.setWindowTitle("Running")
        # Ui_Dialog.setWindowFlags(QtCore.Qt.WindowStaysOnTopHint)
        
        get_pic(loading_gif, 'loading.gif')
        movie = QMovie('loading.gif')
        self.label.setMovie(movie)
        movie.start()  
 

if __name__ == "__main__":
    app = QApplication(sys.argv)
    win = Running_Win()
    win.show()
    sys.exit(app.exec_())