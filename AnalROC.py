# -*- coding: utf-8 -*-

# Form implementation generated from reading ui file 'AnalROC.ui'
#
# Created by: PyQt5 UI code generator 5.15.4
#
# WARNING: Any manual changes made to this file will be lost when pyuic5 is
# run again.  Do not edit this file unless you know what you are doing.


from PyQt5 import QtCore, QtGui, QtWidgets


class Ui_Form(object):
    def setupUi(self, Form):
        Form.setObjectName("Form")
        Form.resize(981, 813)
        self.gridLayout = QtWidgets.QGridLayout(Form)
        self.gridLayout.setObjectName("gridLayout")
        self.horizontalLayout_4 = QtWidgets.QHBoxLayout()
        self.horizontalLayout_4.setObjectName("horizontalLayout_4")
        self.verticalLayout_8 = QtWidgets.QVBoxLayout()
        self.verticalLayout_8.setObjectName("verticalLayout_8")
        self.verticalLayout = QtWidgets.QVBoxLayout()
        self.verticalLayout.setObjectName("verticalLayout")
        self.verticalLayout_6 = QtWidgets.QVBoxLayout()
        self.verticalLayout_6.setObjectName("verticalLayout_6")
        self.label_4 = QtWidgets.QLabel(Form)
        font = QtGui.QFont()
        font.setFamily("Arial")
        font.setPointSize(13)
        self.label_4.setFont(font)
        self.label_4.setObjectName("label_4")
        self.verticalLayout_6.addWidget(self.label_4)
        self.horizontalLayout = QtWidgets.QHBoxLayout()
        self.horizontalLayout.setObjectName("horizontalLayout")
        self.verticalLayout_3 = QtWidgets.QVBoxLayout()
        self.verticalLayout_3.setObjectName("verticalLayout_3")
        self.label = QtWidgets.QLabel(Form)
        font = QtGui.QFont()
        font.setFamily("Times New Roman")
        font.setPointSize(11)
        self.label.setFont(font)
        self.label.setObjectName("label")
        self.verticalLayout_3.addWidget(self.label)
        self.label_7 = QtWidgets.QLabel(Form)
        font = QtGui.QFont()
        font.setFamily("Times New Roman")
        font.setPointSize(11)
        self.label_7.setFont(font)
        self.label_7.setObjectName("label_7")
        self.verticalLayout_3.addWidget(self.label_7)
        self.label_5 = QtWidgets.QLabel(Form)
        font = QtGui.QFont()
        font.setFamily("Times New Roman")
        font.setPointSize(11)
        self.label_5.setFont(font)
        self.label_5.setObjectName("label_5")
        self.verticalLayout_3.addWidget(self.label_5)
        self.label_6 = QtWidgets.QLabel(Form)
        font = QtGui.QFont()
        font.setFamily("Times New Roman")
        font.setPointSize(11)
        self.label_6.setFont(font)
        self.label_6.setObjectName("label_6")
        self.verticalLayout_3.addWidget(self.label_6)
        self.horizontalLayout.addLayout(self.verticalLayout_3)
        self.verticalLayout_2 = QtWidgets.QVBoxLayout()
        self.verticalLayout_2.setObjectName("verticalLayout_2")
        self.spinBoxPub = QtWidgets.QSpinBox(Form)
        self.spinBoxPub.setObjectName("spinBoxPub")
        self.verticalLayout_2.addWidget(self.spinBoxPub)
        self.comboBoxDataset = QtWidgets.QComboBox(Form)
        self.comboBoxDataset.setObjectName("comboBoxDataset")
        self.verticalLayout_2.addWidget(self.comboBoxDataset)
        self.comboBoxDistance = QtWidgets.QComboBox(Form)
        self.comboBoxDistance.setObjectName("comboBoxDistance")
        self.verticalLayout_2.addWidget(self.comboBoxDistance)
        self.spinBoxRandom = QtWidgets.QSpinBox(Form)
        self.spinBoxRandom.setMinimum(0)
        self.spinBoxRandom.setMaximum(999999999)
        self.spinBoxRandom.setProperty("value", 100000)
        self.spinBoxRandom.setObjectName("spinBoxRandom")
        self.verticalLayout_2.addWidget(self.spinBoxRandom)
        self.horizontalLayout.addLayout(self.verticalLayout_2)
        self.verticalLayout_6.addLayout(self.horizontalLayout)
        self.verticalLayout.addLayout(self.verticalLayout_6)
        self.horizontalLayout_2 = QtWidgets.QHBoxLayout()
        self.horizontalLayout_2.setObjectName("horizontalLayout_2")
        self.pushButtonDatabase = QtWidgets.QPushButton(Form)
        self.pushButtonDatabase.setObjectName("pushButtonDatabase")
        self.horizontalLayout_2.addWidget(self.pushButtonDatabase)
        self.pushButtonConfirm = QtWidgets.QPushButton(Form)
        self.pushButtonConfirm.setObjectName("pushButtonConfirm")
        self.horizontalLayout_2.addWidget(self.pushButtonConfirm)
        self.pushButtonCancel = QtWidgets.QPushButton(Form)
        self.pushButtonCancel.setObjectName("pushButtonCancel")
        self.horizontalLayout_2.addWidget(self.pushButtonCancel)
        self.verticalLayout.addLayout(self.horizontalLayout_2)
        self.verticalLayout_8.addLayout(self.verticalLayout)
        self.verticalLayout_4 = QtWidgets.QVBoxLayout()
        self.verticalLayout_4.setObjectName("verticalLayout_4")
        self.label_3 = QtWidgets.QLabel(Form)
        font = QtGui.QFont()
        font.setFamily("Arial")
        font.setPointSize(11)
        self.label_3.setFont(font)
        self.label_3.setObjectName("label_3")
        self.verticalLayout_4.addWidget(self.label_3)
        self.tableView = QtWidgets.QTableView(Form)
        self.tableView.setMinimumSize(QtCore.QSize(0, 200))
        self.tableView.setObjectName("tableView")
        self.verticalLayout_4.addWidget(self.tableView)
        self.horizontalLayout_3 = QtWidgets.QHBoxLayout()
        self.horizontalLayout_3.setObjectName("horizontalLayout_3")
        self.pushButtonPval = QtWidgets.QPushButton(Form)
        self.pushButtonPval.setObjectName("pushButtonPval")
        self.horizontalLayout_3.addWidget(self.pushButtonPval)
        self.pushButtonCurve = QtWidgets.QPushButton(Form)
        self.pushButtonCurve.setObjectName("pushButtonCurve")
        self.horizontalLayout_3.addWidget(self.pushButtonCurve)
        self.verticalLayout_4.addLayout(self.horizontalLayout_3)
        self.verticalLayout_8.addLayout(self.verticalLayout_4)
        self.verticalLayout_5 = QtWidgets.QVBoxLayout()
        self.verticalLayout_5.setObjectName("verticalLayout_5")
        self.groupBox = QtWidgets.QGroupBox(Form)
        self.groupBox.setMinimumSize(QtCore.QSize(0, 250))
        self.groupBox.setMaximumSize(QtCore.QSize(16777215, 300))
        font = QtGui.QFont()
        font.setFamily("Arial")
        font.setPointSize(11)
        self.groupBox.setFont(font)
        self.groupBox.setObjectName("groupBox")
        self.verticalLayout_5.addWidget(self.groupBox)
        self.verticalLayout_8.addLayout(self.verticalLayout_5)
        self.horizontalLayout_4.addLayout(self.verticalLayout_8)
        self.verticalLayout_11 = QtWidgets.QVBoxLayout()
        self.verticalLayout_11.setObjectName("verticalLayout_11")
        self.verticalLayout_10 = QtWidgets.QVBoxLayout()
        self.verticalLayout_10.setObjectName("verticalLayout_10")
        self.verticalLayout_7 = QtWidgets.QVBoxLayout()
        self.verticalLayout_7.setObjectName("verticalLayout_7")
        self.groupBox_1 = QtWidgets.QGroupBox(Form)
        font = QtGui.QFont()
        font.setFamily("Arial")
        font.setPointSize(11)
        self.groupBox_1.setFont(font)
        self.groupBox_1.setObjectName("groupBox_1")
        self.verticalLayout_7.addWidget(self.groupBox_1)
        self.verticalLayout_10.addLayout(self.verticalLayout_7)
        self.verticalLayout_9 = QtWidgets.QVBoxLayout()
        self.verticalLayout_9.setObjectName("verticalLayout_9")
        self.groupBox_2 = QtWidgets.QGroupBox(Form)
        font = QtGui.QFont()
        font.setFamily("Arial")
        font.setPointSize(11)
        self.groupBox_2.setFont(font)
        self.groupBox_2.setObjectName("groupBox_2")
        self.verticalLayout_9.addWidget(self.groupBox_2)
        self.verticalLayout_10.addLayout(self.verticalLayout_9)
        self.verticalLayout_11.addLayout(self.verticalLayout_10)
        self.progressBar = QtWidgets.QProgressBar(Form)
        self.progressBar.setProperty("value", 0)
        self.progressBar.setObjectName("progressBar")
        self.verticalLayout_11.addWidget(self.progressBar)
        self.horizontalLayout_4.addLayout(self.verticalLayout_11)
        self.gridLayout.addLayout(self.horizontalLayout_4, 0, 0, 1, 1)

        self.retranslateUi(Form)
        QtCore.QMetaObject.connectSlotsByName(Form)

    def retranslateUi(self, Form):
        _translate = QtCore.QCoreApplication.translate
        Form.setWindowTitle(_translate("Form", "Form"))
        self.label_4.setText(_translate("Form", "Parameters for protein pair analysis:"))
        self.label.setText(_translate("Form", "Publication cutoff"))
        self.label_7.setText(_translate("Form", "Select dataset for ROC"))
        self.label_5.setText(_translate("Form", "Distance calculation method:"))
        self.label_6.setText(_translate("Form", "Number of generated random paris:"))
        self.pushButtonDatabase.setText(_translate("Form", "Database"))
        self.pushButtonConfirm.setText(_translate("Form", "Confirm"))
        self.pushButtonCancel.setText(_translate("Form", "Cancel"))
        self.label_3.setText(_translate("Form", "Protein pair table"))
        self.pushButtonPval.setText(_translate("Form", "Calc Change"))
        self.pushButtonCurve.setText(_translate("Form", "Show Curve"))
        self.groupBox.setTitle(_translate("Form", "ROC Curve"))
        self.groupBox_1.setTitle(_translate("Form", "Group1 protein curve"))
        self.groupBox_2.setTitle(_translate("Form", "Group2 protein curve"))


if __name__ == "__main__":
    import sys
    app = QtWidgets.QApplication(sys.argv)
    Form = QtWidgets.QWidget()
    ui = Ui_Form()
    ui.setupUi(Form)
    Form.show()
    sys.exit(app.exec_())
