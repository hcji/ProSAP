# -*- coding: utf-8 -*-

# Form implementation generated from reading ui file 'ColumnSelect.ui'
#
# Created by: PyQt5 UI code generator 5.9.2
#
# WARNING! All changes made in this file will be lost!

from PyQt5 import QtCore, QtGui, QtWidgets

class Ui_Form(object):
    def setupUi(self, Form):
        Form.setObjectName("Form")
        Form.resize(564, 470)
        self.listWidget = QtWidgets.QListWidget(Form)
        self.listWidget.setGeometry(QtCore.QRect(30, 60, 351, 381))
        self.listWidget.setObjectName("listWidget")
        self.labelColumnSelect = QtWidgets.QLabel(Form)
        self.labelColumnSelect.setGeometry(QtCore.QRect(30, 20, 381, 16))
        font = QtGui.QFont()
        font.setFamily("Arial")
        font.setPointSize(11)
        self.labelColumnSelect.setFont(font)
        self.labelColumnSelect.setObjectName("labelColumnSelect")
        self.ButtonColumnSelect = QtWidgets.QPushButton(Form)
        self.ButtonColumnSelect.setGeometry(QtCore.QRect(440, 360, 93, 28))
        self.ButtonColumnSelect.setObjectName("ButtonColumnSelect")
        self.ButtonColumnCancel = QtWidgets.QPushButton(Form)
        self.ButtonColumnCancel.setGeometry(QtCore.QRect(440, 400, 93, 28))
        self.ButtonColumnCancel.setObjectName("ButtonColumnCancel")

        self.retranslateUi(Form)
        QtCore.QMetaObject.connectSlotsByName(Form)

    def retranslateUi(self, Form):
        _translate = QtCore.QCoreApplication.translate
        Form.setWindowTitle(_translate("Form", "Form"))
        self.labelColumnSelect.setText(_translate("Form", "Select Abundance Columns for Analysis:"))
        self.ButtonColumnSelect.setText(_translate("Form", "Confirm"))
        self.ButtonColumnCancel.setText(_translate("Form", "Cancel"))


if __name__ == "__main__":
    import sys
    app = QtWidgets.QApplication(sys.argv)
    Form = QtWidgets.QWidget()
    ui = Ui_Form()
    ui.setupUi(Form)
    Form.show()
    sys.exit(app.exec_())

