# -*- coding: utf-8 -*-
"""
Created on Fri Apr 16 16:23:29 2021

@author: hcji
"""


import matplotlib
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler
from scipy.optimize import curve_fit
from adjustText import adjust_text

matplotlib.use("Qt5Agg")
from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg as FigureCanvas
from matplotlib.figure import Figure
from PyQt5 import QtCore, QtGui, QtWidgets
from seaborn import violinplot, boxplot, scatterplot, color_palette, heatmap
from Utils import meltCurve


class MakeFigure(FigureCanvas):
    def __init__(self,width=5, height=5, dpi=300):
        self.fig = Figure(figsize=(width, height), dpi=dpi)
        self.fig.subplots_adjust(top=0.95,bottom=0.2,left=0.15,right=0.85)
        super(MakeFigure,self).__init__(self.fig) 
        self.axes = self.fig.add_subplot(111)
        self.axes.spines['bottom'].set_linewidth(0.5)
        self.axes.spines['left'].set_linewidth(0.5)
        self.axes.spines['right'].set_linewidth(0.5)
        self.axes.spines['top'].set_linewidth(0.5)
        self.axes.tick_params(labelsize=5)
        FigureCanvas.setSizePolicy(self, QtWidgets.QSizePolicy.Expanding, QtWidgets.QSizePolicy.Expanding)
        FigureCanvas.updateGeometry(self)

    
    def ProteinComplexFigure(self, proteinSubunit, proteinData, colNames):
        self.fig.subplots_adjust(top=0.95,bottom=0.2,left=0.15,right=0.8)
        prots = proteinSubunit.split(',')
        prots = [p.replace('(', '') for p in prots]
        prots = [p.replace(')', '') for p in prots]

        try:
            temps = np.array([float(t.replace('T', '')) for t in colNames])
            temps_ = np.arange(temps[0], temps[-1], 0.1)
        except:
            raise ValueError('invalid column names')

        pltData = dict()
        for p in prots:
            if p in list(proteinData['Accession']):
                vals = proteinData.loc[proteinData.loc[:, 'Accession'] == p, colNames]
                pltData[p] = vals.values[0,:]
        
        self.axes.cla()
        for p, vec in pltData.items():
            paras = curve_fit(meltCurve, temps, vec, bounds=(0, [float('inf'), float('inf'), 0.3]))[0]
            self.axes.scatter(temps, vec, marker='.', label = p, s = 15)
            self.axes.plot(temps_, meltCurve(temps_, paras[0], paras[1], paras[2]), lw = 1)
            
        self.axes.set_xlabel('Temperature (℃)', fontsize=6)
        self.axes.set_ylabel('Abundances', fontsize=6)
        self.axes.tick_params(labelsize=5)
        self.axes.legend(fontsize=4.5, bbox_to_anchor=(1,1), loc="upper left")
        self.draw()
    
    
    def AverageTSAFigure(self, proteinData1, proteinData2, colNames):
        try:
            temps = np.array([float(t.replace('T', '')) for t in colNames])
        except:
            raise ValueError('invalid column names')
        
        vec_1 = np.mean(proteinData1.loc[:, colNames], axis = 0)
        vec_2 = np.mean(proteinData2.loc[:, colNames], axis = 0)
        
        self.axes.cla()
        self.axes.plot(temps, vec_1, label = 'Group 1')
        self.axes.plot(temps, vec_2, label = 'Group 2')
        
        self.axes.tick_params(labelsize=4)
        self.axes.set_xlabel('Temperature (℃)', fontsize=5)
        self.axes.set_ylabel('Abundances', fontsize=5)
        self.axes.legend(fontsize=4)    
        self.draw()
        
    
    def SingleTSAFigure(self, proteinData1, proteinData2, colNames, ProteinAccession, proteinData3=None, proteinData4=None):
        try:
            temps = np.array([float(t.replace('T', '')) for t in colNames])
            temps_ = np.arange(temps[0], temps[-1], 0.1)
        except:
            raise ValueError('invalid column names')
            
        vec_1 = proteinData1.loc[proteinData1.loc[:, 'Accession'] == ProteinAccession, colNames].values[0,:]
        vec_2 = proteinData2.loc[proteinData2.loc[:, 'Accession'] == ProteinAccession, colNames].values[0,:]
        
        paras1 = curve_fit(meltCurve, temps, vec_1, bounds=(0, [float('inf'), float('inf'), 0.3]))[0]
        paras2 = curve_fit(meltCurve, temps, vec_2, bounds=(0, [float('inf'), float('inf'), 0.3]))[0]
        
        self.axes.cla()
        self.axes.scatter(temps, vec_1, marker='.', label = 'Group 1_{}'.format(ProteinAccession), color='b', s = 10)
        self.axes.scatter(temps, vec_2, marker='.', label = 'Group 2_{}'.format(ProteinAccession), color='r', s = 10)
        self.axes.plot(temps_, meltCurve(temps_, paras1[0], paras1[1], paras1[2]), color='b', lw=1)
        self.axes.plot(temps_, meltCurve(temps_, paras2[0], paras2[1], paras2[2]), color='r', lw=1)
        
        if (proteinData3 is not None) and (proteinData4 is not None):
            vec_3 = proteinData3.loc[proteinData3.loc[:, 'Accession'] == ProteinAccession, colNames].values[0,:]
            vec_4 = proteinData4.loc[proteinData4.loc[:, 'Accession'] == ProteinAccession, colNames].values[0,:]
            paras3 = curve_fit(meltCurve, temps, vec_3, bounds=(0, [float('inf'), float('inf'), 0.3]))[0]
            paras4 = curve_fit(meltCurve, temps, vec_4, bounds=(0, [float('inf'), float('inf'), 0.3]))[0]
            self.axes.scatter(temps, vec_3, marker='.', label = 'Group 1_r2_{}'.format(ProteinAccession), color='b', s = 10)
            self.axes.scatter(temps, vec_4, marker='.', label = 'Group 2_r2_{}'.format(ProteinAccession), color='r', s = 10)
            self.axes.plot(temps_, meltCurve(temps_, paras3[0], paras3[1], paras3[2]), color='b', linestyle='--', lw=1)
            self.axes.plot(temps_, meltCurve(temps_, paras4[0], paras4[1], paras4[2]), color='r', linestyle='--', lw=1)    
        
        self.axes.tick_params(labelsize=4)
        self.axes.set_xlabel('Temperature (℃)', fontsize=5)
        self.axes.set_ylabel('Abundances', fontsize=5)
        self.axes.legend(fontsize=3)
        self.draw()
        
        
    def RankTSAResults(self, resultTable):
        self.axes.cla()
        self.axes.scatter(1 + np.arange(len(resultTable)), resultTable.loc[:,'Score'], s = 3)
        for i in range(min(len(resultTable.index), 10)):
            j = resultTable.index[i]
            x, y, s = i, resultTable.loc[j, 'Score'], resultTable.loc[j, 'Accession'].split(';')[0]
            self.axes.text(x, y, s, fontsize = 4, color='r')
        self.draw()
    
    
    def RSDHistFigure(self, rsdList):
        rsdList = [i for i in rsdList if not np.isnan(i)]
        self.axes.cla()
        self.axes.hist(rsdList, 100)
        self.axes.tick_params(labelsize=4)
        self.axes.set_xlabel('RSD', fontsize=4)
        self.axes.set_ylabel('Number', fontsize=4)
        self.draw()
    
    
    def ROCFigure(self, fpr, tpr, auroc):
        self.axes.cla()
        self.axes.plot(fpr, tpr, label='AUC = {}'.format(auroc), color = 'red', lw=0.7)
        self.axes.plot([0, 1], [0, 1], color='black', linestyle='--', lw=0.7)
        self.axes.set_xlabel('False Positive Rate', fontsize = 5)
        self.axes.set_ylabel('True Positive Rate', fontsize = 5)
        self.axes.tick_params(labelsize=4)
        self.axes.legend(fontsize=3)
        self.draw()
        
        
    def ProteinPairFigure(self, p1, p2, proteinData, colNames):
        prots = [p1, p2]
        try:
            temps = np.array([float(t.replace('T', '')) for t in colNames])
            temps_ = np.arange(temps[0], temps[-1], 0.1)
        except:
            raise ValueError('invalid column names')

        pltData = dict()
        for p in prots:
            if p in list(proteinData['Accession']):
                vals = proteinData.loc[proteinData.loc[:, 'Accession'] == p, colNames]
                pltData[p] = vals.values[0,:]
                
        self.axes.cla()
        for p, vec in pltData.items():
            paras = curve_fit(meltCurve, temps, vec, bounds=(0, [float('inf'), float('inf'), 0.3]))[0]
            self.axes.scatter(temps, vec, marker='.', label = p, s = 5)
            self.axes.plot(temps_, meltCurve(temps_, paras[0], paras[1], paras[2]), lw=1)
        
        self.axes.set_xlabel('Temperature (℃)', fontsize=6)
        self.axes.set_ylabel('Abundances', fontsize=6)
        self.axes.legend(fontsize=4)
        self.draw()
        
        
    def iTSA_Volcano(self, iTSA_result, fc_thres, pv_thres, show_marker=False):        
        fc = iTSA_result['logFC']
        pv = iTSA_result['-logAdjPval']
        lb = iTSA_result['Accession']
        
        group = []
        for i in range(len(pv)):
            if (abs(fc[i]) < np.log2(fc_thres)) and (pv[i] < -np.log10(pv_thres)):
                group.append('Not sig')
            elif (abs(fc[i]) >= np.log2(fc_thres)) and (pv[i] < -np.log10(pv_thres)):
                group.append('Fold change')
            elif (abs(fc[i]) < np.log2(fc_thres)) and (pv[i] >= -np.log10(pv_thres)):
                group.append('Score')
            else:
                group.append('Both sig')
        pltdata = pd.DataFrame({'LB':lb, 'FC': fc, 'PV': pv, 'G': group})
        # sig = np.where(np.logical_and(np.abs(fc) >= np.log2(fc_thres), pv >= -np.log10(pv_thres)))[0]
        
        self.axes.cla()
        scatterplot(data=pltdata, x="FC", y="PV", hue="G", palette='tab10', legend=False, marker='.', alpha=0.7, edgecolor='none', ax=self.axes)
        
        if show_marker:
            markers = pltdata[pltdata['G'] == 'Both sig']
            markers = markers.iloc[:min(len(markers), 10),:]
            texts = []
            for i in markers.index:
                x, y, s = markers.loc[i, 'FC'], markers.loc[i, 'PV'], markers.loc[i, 'LB'].split(';')[0]
                texts.append(self.axes.text(x, y, s, fontsize=3.5))
        
        '''
        adjust_text(texts, force_points=0.2, force_text=0.2,
                    expand_points=(1, 1), expand_text=(1, 1),
                    arrowprops=dict(arrowstyle="-", color='black', lw=0.5), ax=self.axes)
        '''
        self.axes.axvline(x = np.log2(fc_thres),ls = '--', color = 'black', lw=0.5)
        self.axes.axvline(x = -np.log2(fc_thres),ls = '--', color = 'black', lw=0.5)
        self.axes.axhline(y = -np.log10(pv_thres), ls = '--', color = 'black', lw=0.5)
        self.axes.set_xlabel('Log FC', fontsize = 5)
        self.axes.set_ylabel('-Log Adj P', fontsize = 5)
        self.axes.tick_params(labelsize=4)
        self.draw()
    
    
    def PCAPlot(self, X, y):
        pca = PCA(n_components=2)
        X_s = StandardScaler().fit_transform(X.T)
        X_r = pca.fit(X_s).transform(X_s)
        label = np.unique(y)
        target_names = ['group_{}'.format(i) for i in label]
        
        self.axes.cla()
        for i in range(len(label)):
            self.axes.scatter(X_r[y == label[i], 0], X_r[y == label[i], 1], alpha=.8, lw=1, label=target_names[i], s=10)
        self.axes.set_xlabel('PC 1', fontsize = 5)
        self.axes.set_ylabel('PC 2', fontsize = 5)
        self.axes.tick_params(labelsize=5)
        self.draw()        
    
    
    def BarChart(self, X, y):
        self.axes.cla()
        cm = plt.cm.get_cmap('rainbow')
        flierprops = dict(markersize = 2)
        bplot = self.axes.boxplot(np.log2(X), patch_artist=True, flierprops=flierprops)
        self.axes.set_xticklabels(list(X.columns), rotation = 90) 
        self.axes.set_xlabel('Sample', fontsize = 6)
        self.axes.set_ylabel('Log2 Intensity', fontsize = 6)
        colors = [cm(val / len(X.columns)) for val in range(len(X.columns))]
        for patch, color in zip(bplot['boxes'], colors):
            patch.set_facecolor(color)
        self.draw()
        
    
    def CorrHeatMap(self, X):
        self.axes.cla()
        X_copy = X.copy()
        for i in range(X.shape[0]):
            X_copy.iloc[i,:] /= np.nanmax(X_copy.iloc[i,:])
        
        corr = np.round(np.corrcoef(X_copy.T), 2)
        self.axes.imshow(corr, cmap="YlOrBr")
        self.axes.set_xticks(np.arange(corr.shape[0]))
        self.axes.set_yticks(np.arange(corr.shape[0]))
        self.axes.set_xticklabels(list(X.columns), fontsize = 6, rotation = 90)
        self.axes.set_yticklabels(list(X.columns), fontsize = 6)
        for i in range(corr.shape[0]):
            for j in range(corr.shape[0]):
                self.axes.text(i, j, corr[i, j], ha="center", va="center", color="black", fontsize=3)
        self.draw()
        
        
    def PlotQCRSD(self, dataRSD):
        self.axes.cla()
        violinplot(x="Method", y="RSD", data=dataRSD, ax=self.axes)
        self.draw()
    
    
    def PlotQCBox(self, databox):
        self.axes.cla()
        boxplot(x="Method", y="Values", data=databox, ax=self.axes)
        self.draw()
        
        
    def TPP2D_Volcano(self, fdr_df, hits):
        x = np.sign(fdr_df['slopeH1']) * np.sqrt(fdr_df['rssH0'] - fdr_df['rssH1'])
        y = np.log2(fdr_df['F_statistic'] + 1)
        l = fdr_df['clustername'].values
        
        group = []
        for ll in l:
            if ll in hits['clustername'].values:
                group.append('Hits')
            else:
                group.append('Others')
        pltdata = pd.DataFrame({'x':x, 'y': y, 'l':l, 'G': group})
        # sig = np.where(np.logical_and(np.abs(fc) >= np.log2(fc_thres), pv >= -np.log10(pv_thres)))[0]
        
        self.axes.cla()
        scatterplot(data=pltdata, x="x", y="y", hue="G", palette='tab10', legend=False, alpha=0.7, edgecolor='none', marker='.', ax=self.axes)   
        '''
        markers = pltdata[pltdata['G'] == 'Hits']
        markers = markers.iloc[:min(len(markers), 10),:]
        texts = []
        for i in markers.index:
            x, y, s = markers.loc[i, 'x'], markers.loc[i, 'y'], markers.loc[i, 'l'].split(';')[0]
            texts.append(self.axes.text(x, y, s, fontsize=3))
        
        p = adjust_text(texts, force_points=0.2, force_text=0.2,
                    expand_points=(1, 1), expand_text=(1, 1),
                    arrowprops=dict(arrowstyle="-", color='black', lw=0.5), ax=self.axes)
        '''
        self.axes.set_xlabel('sign(k) sqrt(RSS0-RSS1)', fontsize = 5)
        self.axes.set_ylabel('np.log2 (F_statistic + 1)', fontsize = 5)
        self.axes.tick_params(labelsize=4)
        self.draw()
        
        
    def TPP2D_protHeatmap(self, data, ProteinAccession):
        pltdata = data[data['clustername'] == ProteinAccession]
        conc = np.unique(pltdata['conc'])
        temp = np.unique(pltdata['temperature'])
        img = np.zeros((len(conc), len(temp)))
        for i in pltdata.index:
            a = np.where(conc == pltdata.loc[i,'conc'])[0][0]
            b = np.where(temp == pltdata.loc[i,'temperature'])[0][0]
            img[a, b] = pltdata.loc[i,'rel_value']
        img = pd.DataFrame(img)
        img.index = conc
        img.columns = temp
        heatmap(img, ax=self.axes, cbar=False)
        self.axes.tick_params(labelsize = 6)
        self.axes.set_xlabel('temperture', fontsize = 6)
        self.axes.set_ylabel('drug concentration', fontsize = 6)
        self.draw()
        