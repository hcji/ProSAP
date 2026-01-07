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

    def iTSA_Volcano(self, iTSA_result, fc_thres, pv_thres, show_marker=True, top_n=15):
        """修复版本：解决第一次作图挤在一起和标签重叠问题"""

        fc = iTSA_result['logFC']
        pv = iTSA_result['-logAdjPval']
        lb = iTSA_result['Accession']

        # 清除之前的图形
        self.axes.cla()

        # 确保图形尺寸正确（修复第一次作图挤在一起的问题）
        self.fig.set_size_inches(20, 20, forward=True)
        self.fig.set_dpi(300)
        self.fig.tight_layout()

        # 1. 绘图基础 - 使用简化的三色方案
        colors = np.where((pv >= -np.log10(pv_thres)) & (fc >= np.log2(fc_thres)), '#E41A1C',
                          np.where((pv >= -np.log10(pv_thres)) & (fc <= -np.log2(fc_thres)), '#377EB8',
                                   '#B0B0B0'))

        # 绘制散点图
        scatter = self.axes.scatter(fc, pv,
                                    c=colors,
                                    alpha=0.7,
                                    edgecolors='black',
                                    linewidths=0.3,
                                    s=25,
                                    zorder=5)

        # 2. 添加阈值线
        self.axes.axvline(x=np.log2(fc_thres), ls='--', color='#666666', lw=1, alpha=0.8, zorder=1)
        self.axes.axvline(x=-np.log2(fc_thres), ls='--', color='#666666', lw=1, alpha=0.8, zorder=1)
        self.axes.axhline(y=-np.log10(pv_thres), ls='--', color='#666666', lw=1, alpha=0.8, zorder=1)

        # 3. 设置合理的坐标轴范围（解决挤在一起的问题）
        if len(fc) > 0:
            # 计算合适的范围
            x_min, x_max = np.min(fc), np.max(fc)
            y_min, y_max = np.min(pv), np.max(pv)

            # 确保有足够的空间
            x_range = max(2.0, x_max - x_min)  # 最小x范围2.0
            y_range = max(2.0, y_max - y_min)  # 最小y范围2.0

            x_center = (x_min + x_max) / 2
            y_center = (y_min + y_max) / 2

            self.axes.set_xlim(x_center - x_range / 2 * 1.2, x_center + x_range / 2 * 1.2)
            self.axes.set_ylim(max(0, y_center - y_range / 2 * 1.2), y_center + y_range / 2 * 1.2)

        # 4. 改进的智能标签放置算法
        if show_marker and top_n > 0:
            # 只选择显著变化的点进行标记
            sig_mask = (pv >= -np.log10(pv_thres)) & ((fc >= np.log2(fc_thres)) | (fc <= -np.log2(fc_thres)))
            sig_indices = np.where(sig_mask)[0]

            if len(sig_indices) > 0:
                # 按显著性排序
                sig_scores = pv[sig_mask] * np.abs(fc[sig_mask])
                sorted_indices = sig_indices[np.argsort(-sig_scores)]

                # 限制最大标签数量
                max_labels = min(top_n, len(sorted_indices))

                # 获取当前坐标轴范围
                xlim = self.axes.get_xlim()
                ylim = self.axes.get_ylim()
                x_range = xlim[1] - xlim[0]
                y_range = ylim[1] - ylim[0]

                # 估计标签大小（像素转数据坐标）
                # 使用更保守的估计
                label_width_data = x_range * 0.08  # 8% of x range
                label_height_data = y_range * 0.04  # 4% of y range

                # 使用四叉树空间分区避免重叠（简化版）
                grid_size_x = label_width_data * 2  # 网格大小是标签宽度的2倍
                grid_size_y = label_height_data * 2

                # 创建网格占用表
                x_bins = int(x_range / grid_size_x) + 1
                y_bins = int(y_range / grid_size_y) + 1
                grid_occupied = np.zeros((x_bins, y_bins), dtype=bool)

                texts = []
                used_positions = []

                def get_grid_index(x, y):
                    """将数据坐标转换为网格索引"""
                    x_idx = int((x - xlim[0]) / grid_size_x)
                    y_idx = int((y - ylim[0]) / grid_size_y)
                    return np.clip(x_idx, 0, x_bins - 1), np.clip(y_idx, 0, y_bins - 1)

                def is_position_available(x, y):
                    """检查位置是否可用"""
                    x_idx, y_idx = get_grid_index(x, y)

                    # 检查中心网格
                    if grid_occupied[x_idx, y_idx]:
                        return False

                    # 检查相邻网格
                    for dx in [-1, 0, 1]:
                        for dy in [-1, 0, 1]:
                            nx, ny = x_idx + dx, y_idx + dy
                            if 0 <= nx < x_bins and 0 <= ny < y_bins:
                                if grid_occupied[nx, ny]:
                                    # 检查实际距离
                                    occupied_x = xlim[0] + nx * grid_size_x + grid_size_x / 2
                                    occupied_y = ylim[0] + ny * grid_size_y + grid_size_y / 2
                                    if (abs(x - occupied_x) < label_width_data * 1.2 and
                                            abs(y - occupied_y) < label_height_data * 1.2):
                                        return False
                    return True

                def mark_position_used(x, y):
                    """标记位置为已使用"""
                    x_idx, y_idx = get_grid_index(x, y)
                    grid_occupied[x_idx, y_idx] = True

                # 处理每个要标记的点
                labels_added = 0
                for idx in sorted_indices:
                    if labels_added >= max_labels:
                        break

                    x_data, y_data = fc[idx], pv[idx]
                    protein_name = lb[idx].split(';')[0] if ';' in lb[idx] else lb[idx]

                    # 尝试多个标签位置
                    positions_to_try = []

                    # 生成候选位置：圆形排列
                    n_candidates = 12
                    for i in range(n_candidates):
                        angle = 2 * np.pi * i / n_candidates
                        distance = label_width_data * 2.5  # 起始距离

                        # 逐步增加距离
                        for dist_multiplier in [1.0, 1.5, 2.0, 3.0]:
                            x_candidate = x_data + distance * dist_multiplier * np.cos(angle)
                            y_candidate = y_data + distance * dist_multiplier * np.sin(angle)

                            # 确保在坐标轴范围内
                            if (xlim[0] <= x_candidate <= xlim[1] and
                                    ylim[0] <= y_candidate <= ylim[1]):
                                positions_to_try.append((x_candidate, y_candidate, angle))

                    # 按距离排序（先尝试近的）
                    positions_to_try.sort(key=lambda pos: ((pos[0] - x_data) ** 2 + (pos[1] - y_data) ** 2))

                    # 尝试每个候选位置
                    position_found = False
                    for x_try, y_try, angle in positions_to_try:
                        if is_position_available(x_try, y_try):
                            # 确定对齐方式
                            if -np.pi / 4 <= angle < np.pi / 4:
                                ha = 'left'
                            elif np.pi / 4 <= angle < 3 * np.pi / 4:
                                ha = 'center'
                            else:
                                ha = 'right'

                            if -np.pi / 2 <= angle < np.pi / 2:
                                va = 'bottom'
                            else:
                                va = 'top'

                            # 创建标签
                            text = self.axes.text(x_try, y_try, protein_name,
                                                  fontsize=4,
                                                  verticalalignment=va,
                                                  horizontalalignment=ha,
                                                  bbox=dict(boxstyle="round,pad=0.2",
                                                            facecolor="yellow",
                                                            alpha=0.9,
                                                            edgecolor='black',
                                                            linewidth=0.5),
                                                  zorder=10)

                            # 添加连接线
                            self.axes.plot([x_data, x_try], [y_data, y_try],
                                           color='gray',
                                           linewidth=0.5,
                                           alpha=0.6,
                                           zorder=6,
                                           linestyle='--')

                            texts.append(text)
                            used_positions.append((x_try, y_try))
                            mark_position_used(x_try, y_try)
                            position_found = True
                            labels_added += 1
                            break

                    if not position_found:
                        # 如果找不到位置，跳过这个标签
                        continue

        # 5. 设置标签和标题
        self.axes.set_xlabel('Log2 Fold Change\n← Destabilized | Stabilized →',
                             fontsize=8, labelpad=8)
        self.axes.set_ylabel('-Log10 Adjusted P-value\n(Statistical Significance)',
                             fontsize=8, labelpad=8)
        self.axes.tick_params(labelsize=7)

        # 6. 创建图例
        from matplotlib.patches import Patch

        # 清除现有图例
        if self.axes.get_legend():
            self.axes.get_legend().remove()

        # 创建图例句柄
        legend_elements = [
            Patch(facecolor='#E41A1C', alpha=0.7, edgecolor='black', linewidth=0.3,
                  label='Thermally stabilized'),
            Patch(facecolor='#377EB8', alpha=0.7, edgecolor='black', linewidth=0.3,
                  label='Thermally destabilized'),
            Patch(facecolor='#B0B0B0', alpha=0.7, edgecolor='black', linewidth=0.3,
                  label='Not significant')
        ]

        # 添加图例（放在图外右侧）
        legend = self.axes.legend(handles=legend_elements,
                                  fontsize=6,
                                  loc='upper left',
                                  bbox_to_anchor=(1.02, 1.0),
                                  borderaxespad=0.0,
                                  framealpha=0.2,  # 更低的透明度
                                  fancybox=True,
                                  handlelength=1.0,
                                  handleheight=1.0,
                                  labelspacing=0.2)

        # 7. 添加统计信息框
        total_proteins = len(iTSA_result)
        stabilized = np.sum((pv >= -np.log10(pv_thres)) & (fc >= np.log2(fc_thres)))
        destabilized = np.sum((pv >= -np.log10(pv_thres)) & (fc <= -np.log2(fc_thres)))

        stats_text = (f"Total: {total_proteins}\n"
                      f"Stabilized: {stabilized}\n"
                      f"Destabilized: {destabilized}")

        self.axes.text(0.02, 0.98, stats_text,
                       transform=self.axes.transAxes,
                       fontsize=5,
                       verticalalignment='top',
                       horizontalalignment='left',
                       bbox=dict(boxstyle="round,pad=0.3",
                                 facecolor="white",
                                 alpha=0.7,
                                 edgecolor='gray',
                                 linewidth=0.5))

        # 8. 设置网格
        self.axes.grid(True, alpha=0.15, linestyle=':', linewidth=0.5, zorder=0)

        # 9. 设置坐标轴背景色
        self.axes.set_facecolor('white')

        # 10. 强制重绘并调整布局
        self.fig.canvas.draw_idle()
        self.fig.tight_layout(rect=[0, 0, 0.82, 1])  # 右侧留出18%空间

        # 11. 立即重绘
        self.draw()

        return len(texts) if 'texts' in locals() else 0
    
    
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
        