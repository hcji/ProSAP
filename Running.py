# -*- coding: utf-8 -*-
"""
Created on Mon Jun 21 15:53:38 2021

@author: jihon
"""


from PyQt5 import QtCore, QtGui, QtWidgets

class RunningDialog(object):
    def setupUi(self, Form):
        Form.setObjectName("Form")
        Form.resize(600, 350)
        
        self.label = QtWidgets.QLabel(Form)
        self.label.setGeometry(QtCore.QRect(0, 0, 600, 350))
        self.label.setObjectName("label")
        
        self.gif = QtGui.QMovie('img/loading.gif')
        self.label.setMovie(self.gif)
        self.gif.start()


if __name__ == "__main__":
    import sys
    app = QtWidgets.QApplication(sys.argv)
    Form = QtWidgets.QWidget()
    ui = RunningDialog()
    ui.setupUi(Form)
    Form.show()
    sys.exit(app.exec_())