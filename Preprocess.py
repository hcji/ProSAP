# -*- coding: utf-8 -*-

# Form implementation generated from reading ui file 'Preprocess.ui'
#
# Created by: PyQt5 UI code generator 5.15.4
#
# WARNING: Any manual changes made to this file will be lost when pyuic5 is
# run again.  Do not edit this file unless you know what you are doing.


from PyQt5 import QtCore, QtGui, QtWidgets


class Ui_Form(object):
    def setupUi(self, Form):
        Form.setObjectName("Form")
        Form.resize(992, 775)
        self.gridLayout_2 = QtWidgets.QGridLayout(Form)
        self.gridLayout_2.setObjectName("gridLayout_2")
        self.verticalLayout_5 = QtWidgets.QVBoxLayout()
        self.verticalLayout_5.setObjectName("verticalLayout_5")
        self.horizontalLayout = QtWidgets.QHBoxLayout()
        self.horizontalLayout.setObjectName("horizontalLayout")
        self.gridLayout = QtWidgets.QGridLayout()
        self.gridLayout.setObjectName("gridLayout")
        self.label_3 = QtWidgets.QLabel(Form)
        font = QtGui.QFont()
        font.setFamily("Times New Roman")
        font.setPointSize(10)
        self.label_3.setFont(font)
        self.label_3.setObjectName("label_3")
        self.gridLayout.addWidget(self.label_3, 5, 2, 1, 2)
        self.comboBoxMV = QtWidgets.QComboBox(Form)
        self.comboBoxMV.setObjectName("comboBoxMV")
        self.gridLayout.addWidget(self.comboBoxMV, 11, 0, 1, 2)
        self.comboBoxNorm = QtWidgets.QComboBox(Form)
        self.comboBoxNorm.setObjectName("comboBoxNorm")
        self.gridLayout.addWidget(self.comboBoxNorm, 6, 2, 1, 2)
        self.label_11 = QtWidgets.QLabel(Form)
        font = QtGui.QFont()
        font.setFamily("Times New Roman")
        font.setPointSize(10)
        self.label_11.setFont(font)
        self.label_11.setObjectName("label_11")
        self.gridLayout.addWidget(self.label_11, 10, 2, 1, 2)
        self.pushButtonSave = QtWidgets.QPushButton(Form)
        self.pushButtonSave.setObjectName("pushButtonSave")
        self.gridLayout.addWidget(self.pushButtonSave, 16, 3, 1, 1)
        self.label_1 = QtWidgets.QLabel(Form)
        font = QtGui.QFont()
        font.setFamily("Arial")
        font.setPointSize(11)
        self.label_1.setFont(font)
        self.label_1.setObjectName("label_1")
        self.gridLayout.addWidget(self.label_1, 0, 0, 1, 4)
        self.pushButtonConfirm = QtWidgets.QPushButton(Form)
        self.pushButtonConfirm.setObjectName("pushButtonConfirm")
        self.gridLayout.addWidget(self.pushButtonConfirm, 16, 1, 1, 1)
        self.pushButtonOpen = QtWidgets.QPushButton(Form)
        self.pushButtonOpen.setObjectName("pushButtonOpen")
        self.gridLayout.addWidget(self.pushButtonOpen, 16, 0, 1, 1)
        self.pushButtonClear = QtWidgets.QPushButton(Form)
        self.pushButtonClear.setObjectName("pushButtonClear")
        self.gridLayout.addWidget(self.pushButtonClear, 16, 2, 1, 1)
        self.label_5 = QtWidgets.QLabel(Form)
        font = QtGui.QFont()
        font.setFamily("Times New Roman")
        font.setPointSize(10)
        self.label_5.setFont(font)
        self.label_5.setAlignment(QtCore.Qt.AlignLeading|QtCore.Qt.AlignLeft|QtCore.Qt.AlignVCenter)
        self.label_5.setObjectName("label_5")
        self.gridLayout.addWidget(self.label_5, 1, 0, 1, 2)
        self.comboBoxPSM = QtWidgets.QComboBox(Form)
        self.comboBoxPSM.setContextMenuPolicy(QtCore.Qt.DefaultContextMenu)
        self.comboBoxPSM.setCurrentText("")
        self.comboBoxPSM.setMinimumContentsLength(0)
        self.comboBoxPSM.setObjectName("comboBoxPSM")
        self.gridLayout.addWidget(self.comboBoxPSM, 2, 0, 1, 2)
        self.label_7 = QtWidgets.QLabel(Form)
        font = QtGui.QFont()
        font.setFamily("Times New Roman")
        font.setPointSize(10)
        self.label_7.setFont(font)
        self.label_7.setContextMenuPolicy(QtCore.Qt.NoContextMenu)
        self.label_7.setAlignment(QtCore.Qt.AlignLeading|QtCore.Qt.AlignLeft|QtCore.Qt.AlignVCenter)
        self.label_7.setObjectName("label_7")
        self.gridLayout.addWidget(self.label_7, 5, 0, 1, 2)
        self.label_6 = QtWidgets.QLabel(Form)
        font = QtGui.QFont()
        font.setFamily("Times New Roman")
        font.setPointSize(10)
        self.label_6.setFont(font)
        self.label_6.setAlignment(QtCore.Qt.AlignLeading|QtCore.Qt.AlignLeft|QtCore.Qt.AlignVCenter)
        self.label_6.setObjectName("label_6")
        self.gridLayout.addWidget(self.label_6, 1, 2, 1, 2)
        self.spinBoxPSMFilter = QtWidgets.QSpinBox(Form)
        self.spinBoxPSMFilter.setProperty("value", 2)
        self.spinBoxPSMFilter.setObjectName("spinBoxPSMFilter")
        self.gridLayout.addWidget(self.spinBoxPSMFilter, 2, 2, 1, 2)
        self.comboBoxMerging = QtWidgets.QComboBox(Form)
        self.comboBoxMerging.setObjectName("comboBoxMerging")
        self.gridLayout.addWidget(self.comboBoxMerging, 6, 0, 1, 2)
        self.label = QtWidgets.QLabel(Form)
        font = QtGui.QFont()
        font.setFamily("Times New Roman")
        font.setPointSize(10)
        self.label.setFont(font)
        self.label.setObjectName("label")
        self.gridLayout.addWidget(self.label, 10, 0, 1, 2)
        self.doubleSpinBoxMVFilter = QtWidgets.QDoubleSpinBox(Form)
        self.doubleSpinBoxMVFilter.setSingleStep(0.02)
        self.doubleSpinBoxMVFilter.setProperty("value", 0.2)
        self.doubleSpinBoxMVFilter.setObjectName("doubleSpinBoxMVFilter")
        self.gridLayout.addWidget(self.doubleSpinBoxMVFilter, 11, 2, 1, 2)
        self.doubleSpinBoxRSDFilter = QtWidgets.QDoubleSpinBox(Form)
        self.doubleSpinBoxRSDFilter.setSingleStep(0.02)
        self.doubleSpinBoxRSDFilter.setProperty("value", 0.2)
        self.doubleSpinBoxRSDFilter.setObjectName("doubleSpinBoxRSDFilter")
        self.gridLayout.addWidget(self.doubleSpinBoxRSDFilter, 13, 0, 1, 2)
        self.label_8 = QtWidgets.QLabel(Form)
        font = QtGui.QFont()
        font.setFamily("Times New Roman")
        font.setPointSize(10)
        self.label_8.setFont(font)
        self.label_8.setAlignment(QtCore.Qt.AlignLeading|QtCore.Qt.AlignLeft|QtCore.Qt.AlignVCenter)
        self.label_8.setObjectName("label_8")
        self.gridLayout.addWidget(self.label_8, 12, 0, 1, 2)
        self.label_10 = QtWidgets.QLabel(Form)
        font = QtGui.QFont()
        font.setFamily("Times New Roman")
        font.setPointSize(10)
        self.label_10.setFont(font)
        self.label_10.setContextMenuPolicy(QtCore.Qt.NoContextMenu)
        self.label_10.setAlignment(QtCore.Qt.AlignLeading|QtCore.Qt.AlignLeft|QtCore.Qt.AlignVCenter)
        self.label_10.setObjectName("label_10")
        self.gridLayout.addWidget(self.label_10, 12, 2, 1, 2)
        self.comboBoxReference = QtWidgets.QComboBox(Form)
        self.comboBoxReference.setObjectName("comboBoxReference")
        self.gridLayout.addWidget(self.comboBoxReference, 13, 2, 1, 2)
        self.horizontalLayout.addLayout(self.gridLayout)
        self.splitter = QtWidgets.QSplitter(Form)
        self.splitter.setOrientation(QtCore.Qt.Horizontal)
        self.splitter.setObjectName("splitter")
        self.layoutWidget = QtWidgets.QWidget(self.splitter)
        self.layoutWidget.setObjectName("layoutWidget")
        self.verticalLayout = QtWidgets.QVBoxLayout(self.layoutWidget)
        self.verticalLayout.setContentsMargins(0, 0, 0, 0)
        self.verticalLayout.setObjectName("verticalLayout")
        self.label_2 = QtWidgets.QLabel(self.layoutWidget)
        font = QtGui.QFont()
        font.setFamily("Arial")
        font.setPointSize(11)
        self.label_2.setFont(font)
        self.label_2.setObjectName("label_2")
        self.verticalLayout.addWidget(self.label_2)
        self.ListFile = QtWidgets.QListWidget(self.layoutWidget)
        self.ListFile.setObjectName("ListFile")
        self.verticalLayout.addWidget(self.ListFile)
        self.layoutWidget1 = QtWidgets.QWidget(self.splitter)
        self.layoutWidget1.setObjectName("layoutWidget1")
        self.verticalLayout_2 = QtWidgets.QVBoxLayout(self.layoutWidget1)
        self.verticalLayout_2.setContentsMargins(0, 0, 0, 0)
        self.verticalLayout_2.setObjectName("verticalLayout_2")
        self.label_9 = QtWidgets.QLabel(self.layoutWidget1)
        font = QtGui.QFont()
        font.setFamily("Arial")
        font.setPointSize(11)
        self.label_9.setFont(font)
        self.label_9.setObjectName("label_9")
        self.verticalLayout_2.addWidget(self.label_9)
        self.tableWidgetTemp = QtWidgets.QTableWidget(self.layoutWidget1)
        self.tableWidgetTemp.setObjectName("tableWidgetTemp")
        self.tableWidgetTemp.setColumnCount(0)
        self.tableWidgetTemp.setRowCount(0)
        self.verticalLayout_2.addWidget(self.tableWidgetTemp)
        self.horizontalLayout.addWidget(self.splitter)
        self.verticalLayout_5.addLayout(self.horizontalLayout)
        self.splitter_2 = QtWidgets.QSplitter(Form)
        self.splitter_2.setOrientation(QtCore.Qt.Horizontal)
        self.splitter_2.setObjectName("splitter_2")
        self.layoutWidget2 = QtWidgets.QWidget(self.splitter_2)
        self.layoutWidget2.setObjectName("layoutWidget2")
        self.verticalLayout_4 = QtWidgets.QVBoxLayout(self.layoutWidget2)
        self.verticalLayout_4.setContentsMargins(0, 0, 0, 0)
        self.verticalLayout_4.setObjectName("verticalLayout_4")
        self.groupBox = QtWidgets.QGroupBox(self.layoutWidget2)
        font = QtGui.QFont()
        font.setFamily("Arial")
        font.setPointSize(11)
        self.groupBox.setFont(font)
        self.groupBox.setObjectName("groupBox")
        self.verticalLayout_4.addWidget(self.groupBox)
        self.layoutWidget3 = QtWidgets.QWidget(self.splitter_2)
        self.layoutWidget3.setObjectName("layoutWidget3")
        self.verticalLayout_3 = QtWidgets.QVBoxLayout(self.layoutWidget3)
        self.verticalLayout_3.setContentsMargins(0, 0, 0, 0)
        self.verticalLayout_3.setObjectName("verticalLayout_3")
        self.label_4 = QtWidgets.QLabel(self.layoutWidget3)
        font = QtGui.QFont()
        font.setFamily("Arial")
        font.setPointSize(11)
        self.label_4.setFont(font)
        self.label_4.setObjectName("label_4")
        self.verticalLayout_3.addWidget(self.label_4)
        self.tableView = QtWidgets.QTableView(self.layoutWidget3)
        self.tableView.setMinimumSize(QtCore.QSize(500, 0))
        self.tableView.setObjectName("tableView")
        self.verticalLayout_3.addWidget(self.tableView)
        self.progressBar = QtWidgets.QProgressBar(self.layoutWidget3)
        self.progressBar.setProperty("value", 0)
        self.progressBar.setObjectName("progressBar")
        self.verticalLayout_3.addWidget(self.progressBar)
        self.verticalLayout_5.addWidget(self.splitter_2)
        self.gridLayout_2.addLayout(self.verticalLayout_5, 0, 0, 1, 1)

        self.retranslateUi(Form)
        QtCore.QMetaObject.connectSlotsByName(Form)

    def retranslateUi(self, Form):
        _translate = QtCore.QCoreApplication.translate
        Form.setWindowTitle(_translate("Form", "Form"))
        self.label_3.setText(_translate("Form", "Normlization"))
        self.label_11.setText(_translate("Form", "Missing value filter ratio"))
        self.pushButtonSave.setText(_translate("Form", "Save as..."))
        self.label_1.setText(_translate("Form", "Parameters"))
        self.pushButtonConfirm.setText(_translate("Form", "Confirm"))
        self.pushButtonOpen.setText(_translate("Form", "Open"))
        self.pushButtonClear.setText(_translate("Form", "Clear"))
        self.label_5.setText(_translate("Form", "PSM column name"))
        self.label_7.setText(_translate("Form", "Merging metrics"))
        self.label_6.setText(_translate("Form", "PSM filter threshold"))
        self.label.setText(_translate("Form", "Missing value imputation"))
        self.label_8.setText(_translate("Form", "RSD filter threshold"))
        self.label_10.setText(_translate("Form", "Reference columns"))
        self.label_2.setText(_translate("Form", "File List"))
        self.label_9.setText(_translate("Form", "Temperature Seting"))
        self.groupBox.setTitle(_translate("Form", "RSD Distribution"))
        self.label_4.setText(_translate("Form", "Data Viewer"))


if __name__ == "__main__":
    import sys
    app = QtWidgets.QApplication(sys.argv)
    Form = QtWidgets.QWidget()
    ui = Ui_Form()
    ui.setupUi(Form)
    Form.show()
    sys.exit(app.exec_())
