# -*- coding: utf-8 -*-

# Form implementation generated from reading ui file 'AnalTPP2D.ui'
#
# Created by: PyQt5 UI code generator 5.15.4
#
# WARNING: Any manual changes made to this file will be lost when pyuic5 is
# run again.  Do not edit this file unless you know what you are doing.


from PyQt5 import QtCore, QtGui, QtWidgets


class Ui_Form(object):
    def setupUi(self, Form):
        Form.setObjectName("Form")
        Form.resize(1101, 782)
        self.gridLayout = QtWidgets.QGridLayout(Form)
        self.gridLayout.setObjectName("gridLayout")
        self.splitter = QtWidgets.QSplitter(Form)
        self.splitter.setOrientation(QtCore.Qt.Horizontal)
        self.splitter.setObjectName("splitter")
        self.widget = QtWidgets.QWidget(self.splitter)
        self.widget.setObjectName("widget")
        self.verticalLayout = QtWidgets.QVBoxLayout(self.widget)
        self.verticalLayout.setContentsMargins(0, 0, 0, 0)
        self.verticalLayout.setObjectName("verticalLayout")
        self.label_1 = QtWidgets.QLabel(self.widget)
        font = QtGui.QFont()
        font.setFamily("Times New Roman")
        font.setPointSize(12)
        font.setBold(False)
        font.setWeight(50)
        self.label_1.setFont(font)
        self.label_1.setObjectName("label_1")
        self.verticalLayout.addWidget(self.label_1)
        self.horizontalLayout = QtWidgets.QHBoxLayout()
        self.horizontalLayout.setObjectName("horizontalLayout")
        self.label_2 = QtWidgets.QLabel(self.widget)
        font = QtGui.QFont()
        font.setFamily("Times New Roman")
        font.setPointSize(11)
        self.label_2.setFont(font)
        self.label_2.setObjectName("label_2")
        self.horizontalLayout.addWidget(self.label_2)
        self.comboBoxMethod = QtWidgets.QComboBox(self.widget)
        self.comboBoxMethod.setObjectName("comboBoxMethod")
        self.horizontalLayout.addWidget(self.comboBoxMethod)
        self.verticalLayout.addLayout(self.horizontalLayout)
        self.horizontalLayout_2 = QtWidgets.QHBoxLayout()
        self.horizontalLayout_2.setObjectName("horizontalLayout_2")
        self.label_9 = QtWidgets.QLabel(self.widget)
        font = QtGui.QFont()
        font.setFamily("Times New Roman")
        font.setPointSize(11)
        self.label_9.setFont(font)
        self.label_9.setObjectName("label_9")
        self.horizontalLayout_2.addWidget(self.label_9)
        self.spinBoxMaxIt = QtWidgets.QSpinBox(self.widget)
        self.spinBoxMaxIt.setMinimum(2)
        self.spinBoxMaxIt.setMaximum(1000)
        self.spinBoxMaxIt.setProperty("value", 200)
        self.spinBoxMaxIt.setObjectName("spinBoxMaxIt")
        self.horizontalLayout_2.addWidget(self.spinBoxMaxIt)
        self.verticalLayout.addLayout(self.horizontalLayout_2)
        self.horizontalLayout_3 = QtWidgets.QHBoxLayout()
        self.horizontalLayout_3.setObjectName("horizontalLayout_3")
        self.label_12 = QtWidgets.QLabel(self.widget)
        font = QtGui.QFont()
        font.setFamily("Times New Roman")
        font.setPointSize(11)
        self.label_12.setFont(font)
        self.label_12.setObjectName("label_12")
        self.horizontalLayout_3.addWidget(self.label_12)
        self.spinBoxBoost = QtWidgets.QSpinBox(self.widget)
        self.spinBoxBoost.setMinimum(2)
        self.spinBoxBoost.setMaximum(1000)
        self.spinBoxBoost.setProperty("value", 2)
        self.spinBoxBoost.setObjectName("spinBoxBoost")
        self.horizontalLayout_3.addWidget(self.spinBoxBoost)
        self.verticalLayout.addLayout(self.horizontalLayout_3)
        self.horizontalLayout_6 = QtWidgets.QHBoxLayout()
        self.horizontalLayout_6.setObjectName("horizontalLayout_6")
        self.label_10 = QtWidgets.QLabel(self.widget)
        font = QtGui.QFont()
        font.setFamily("Times New Roman")
        font.setPointSize(11)
        self.label_10.setFont(font)
        self.label_10.setObjectName("label_10")
        self.horizontalLayout_6.addWidget(self.label_10)
        self.doubleSpinBoxAlpha = QtWidgets.QDoubleSpinBox(self.widget)
        self.doubleSpinBoxAlpha.setSingleStep(0.01)
        self.doubleSpinBoxAlpha.setProperty("value", 0.1)
        self.doubleSpinBoxAlpha.setObjectName("doubleSpinBoxAlpha")
        self.horizontalLayout_6.addWidget(self.doubleSpinBoxAlpha)
        self.verticalLayout.addLayout(self.horizontalLayout_6)
        self.horizontalLayout_5 = QtWidgets.QHBoxLayout()
        self.horizontalLayout_5.setObjectName("horizontalLayout_5")
        self.pushButtonData = QtWidgets.QPushButton(self.widget)
        font = QtGui.QFont()
        font.setFamily("Times New Roman")
        font.setPointSize(11)
        self.pushButtonData.setFont(font)
        self.pushButtonData.setObjectName("pushButtonData")
        self.horizontalLayout_5.addWidget(self.pushButtonData)
        self.pushButtonOK = QtWidgets.QPushButton(self.widget)
        font = QtGui.QFont()
        font.setFamily("Times New Roman")
        font.setPointSize(11)
        self.pushButtonOK.setFont(font)
        self.pushButtonOK.setObjectName("pushButtonOK")
        self.horizontalLayout_5.addWidget(self.pushButtonOK)
        self.pushButtonClose = QtWidgets.QPushButton(self.widget)
        font = QtGui.QFont()
        font.setFamily("Times New Roman")
        font.setPointSize(11)
        self.pushButtonClose.setFont(font)
        self.pushButtonClose.setObjectName("pushButtonClose")
        self.horizontalLayout_5.addWidget(self.pushButtonClose)
        self.verticalLayout.addLayout(self.horizontalLayout_5)
        self.labelDatabase = QtWidgets.QLabel(self.widget)
        font = QtGui.QFont()
        font.setFamily("Times New Roman")
        font.setPointSize(12)
        self.labelDatabase.setFont(font)
        self.labelDatabase.setObjectName("labelDatabase")
        self.verticalLayout.addWidget(self.labelDatabase)
        self.tableWidgetProteinList = QtWidgets.QTableWidget(self.widget)
        self.tableWidgetProteinList.setMinimumSize(QtCore.QSize(0, 500))
        self.tableWidgetProteinList.setObjectName("tableWidgetProteinList")
        self.tableWidgetProteinList.setColumnCount(0)
        self.tableWidgetProteinList.setRowCount(0)
        self.verticalLayout.addWidget(self.tableWidgetProteinList)
        self.horizontalLayout_4 = QtWidgets.QHBoxLayout()
        self.horizontalLayout_4.setObjectName("horizontalLayout_4")
        spacerItem = QtWidgets.QSpacerItem(40, 20, QtWidgets.QSizePolicy.Expanding, QtWidgets.QSizePolicy.Minimum)
        self.horizontalLayout_4.addItem(spacerItem)
        self.pushButtonSave = QtWidgets.QPushButton(self.widget)
        font = QtGui.QFont()
        font.setFamily("Times New Roman")
        font.setPointSize(11)
        self.pushButtonSave.setFont(font)
        self.pushButtonSave.setObjectName("pushButtonSave")
        self.horizontalLayout_4.addWidget(self.pushButtonSave)
        self.pushButtonShow = QtWidgets.QPushButton(self.widget)
        font = QtGui.QFont()
        font.setFamily("Times New Roman")
        font.setPointSize(11)
        self.pushButtonShow.setFont(font)
        self.pushButtonShow.setObjectName("pushButtonShow")
        self.horizontalLayout_4.addWidget(self.pushButtonShow)
        self.verticalLayout.addLayout(self.horizontalLayout_4)
        self.widget1 = QtWidgets.QWidget(self.splitter)
        self.widget1.setObjectName("widget1")
        self.verticalLayout_5 = QtWidgets.QVBoxLayout(self.widget1)
        self.verticalLayout_5.setContentsMargins(0, 0, 0, 0)
        self.verticalLayout_5.setObjectName("verticalLayout_5")
        self.verticalLayout_2 = QtWidgets.QVBoxLayout()
        self.verticalLayout_2.setObjectName("verticalLayout_2")
        self.label_11 = QtWidgets.QLabel(self.widget1)
        font = QtGui.QFont()
        font.setFamily("Times New Roman")
        font.setPointSize(12)
        self.label_11.setFont(font)
        self.label_11.setObjectName("label_11")
        self.verticalLayout_2.addWidget(self.label_11)
        self.tableViewData = QtWidgets.QTableView(self.widget1)
        self.tableViewData.setMinimumSize(QtCore.QSize(0, 250))
        self.tableViewData.setObjectName("tableViewData")
        self.verticalLayout_2.addWidget(self.tableViewData)
        self.verticalLayout_5.addLayout(self.verticalLayout_2)
        self.tabWidget = QtWidgets.QTabWidget(self.widget1)
        self.tabWidget.setMinimumSize(QtCore.QSize(0, 300))
        font = QtGui.QFont()
        font.setFamily("Times New Roman")
        font.setPointSize(12)
        self.tabWidget.setFont(font)
        self.tabWidget.setObjectName("tabWidget")
        self.vocano = QtWidgets.QWidget()
        font = QtGui.QFont()
        font.setPointSize(12)
        self.vocano.setFont(font)
        self.vocano.setObjectName("vocano")
        self.gridLayout_5 = QtWidgets.QGridLayout(self.vocano)
        self.gridLayout_5.setObjectName("gridLayout_5")
        self.groupBoxVolcano = QtWidgets.QGroupBox(self.vocano)
        self.groupBoxVolcano.setMinimumSize(QtCore.QSize(0, 0))
        font = QtGui.QFont()
        font.setFamily("Arial")
        font.setPointSize(13)
        self.groupBoxVolcano.setFont(font)
        self.groupBoxVolcano.setTitle("")
        self.groupBoxVolcano.setObjectName("groupBoxVolcano")
        self.gridLayout_5.addWidget(self.groupBoxVolcano, 0, 0, 1, 1)
        self.tabWidget.addTab(self.vocano, "")
        self.heatmap = QtWidgets.QWidget()
        self.heatmap.setObjectName("heatmap")
        self.gridLayout_2 = QtWidgets.QGridLayout(self.heatmap)
        self.gridLayout_2.setObjectName("gridLayout_2")
        self.groupBoxHeatmap = QtWidgets.QGroupBox(self.heatmap)
        self.groupBoxHeatmap.setMinimumSize(QtCore.QSize(0, 0))
        font = QtGui.QFont()
        font.setFamily("Arial")
        font.setPointSize(13)
        self.groupBoxHeatmap.setFont(font)
        self.groupBoxHeatmap.setTitle("")
        self.groupBoxHeatmap.setObjectName("groupBoxHeatmap")
        self.gridLayout_2.addWidget(self.groupBoxHeatmap, 0, 0, 1, 1)
        self.tabWidget.addTab(self.heatmap, "")
        self.verticalLayout_5.addWidget(self.tabWidget)
        self.gridLayout.addWidget(self.splitter, 0, 0, 1, 1)

        self.retranslateUi(Form)
        self.tabWidget.setCurrentIndex(0)
        QtCore.QMetaObject.connectSlotsByName(Form)

    def retranslateUi(self, Form):
        _translate = QtCore.QCoreApplication.translate
        Form.setWindowTitle(_translate("Form", "Form"))
        self.label_1.setText(_translate("Form", "Parameters for TPP2D analysis:"))
        self.label_2.setText(_translate("Form", "Select methods"))
        self.label_9.setText(_translate("Form", "Maximum Iter"))
        self.label_12.setText(_translate("Form", "Bootstrap Round"))
        self.label_10.setText(_translate("Form", "Significance Threshold"))
        self.pushButtonData.setText(_translate("Form", "Data"))
        self.pushButtonOK.setText(_translate("Form", "Confirm"))
        self.pushButtonClose.setText(_translate("Form", "Cancel"))
        self.labelDatabase.setText(_translate("Form", "Hit Targets Table"))
        self.pushButtonSave.setText(_translate("Form", "Save"))
        self.pushButtonShow.setText(_translate("Form", "Show Result"))
        self.label_11.setText(_translate("Form", "Data viewer"))
        self.tabWidget.setTabText(self.tabWidget.indexOf(self.vocano), _translate("Form", "Volcano"))
        self.tabWidget.setTabText(self.tabWidget.indexOf(self.heatmap), _translate("Form", "Heatmap"))


if __name__ == "__main__":
    import sys
    app = QtWidgets.QApplication(sys.argv)
    Form = QtWidgets.QWidget()
    ui = Ui_Form()
    ui.setupUi(Form)
    Form.show()
    sys.exit(app.exec_())
