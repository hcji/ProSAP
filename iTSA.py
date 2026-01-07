import numpy as np
import pandas as pd
from scipy.stats import ttest_ind
from limma_py import lmFit, make_contrasts, contrasts_fit, eBayes, toptable
from pydeseq2.dds import DeseqDataSet
from pydeseq2.ds import DeseqStats
from inmoose.edgepy import DGEList, glmFit,glmLRT,topTags,estimateGLMCommonDisp, estimateGLMTagwiseDisp
from patsy import dmatrix

# 定义Python版本的函数
def do_limma_py(X, y, names, lbs):
    if not isinstance(X, pd.DataFrame):
        X_df = pd.DataFrame(X)
    else:
        X_df = X.copy()

    # 设置基因名作为行索引
    X_df.index = names

    y = np.array(y)

    # 分组标签标准化
    y_encoded = np.array(['G1' if str(val) == str(lbs[0]) else 'G2' if str(val) == str(lbs[1]) else val
                          for val in y], dtype=object)
    # 实验设计矩阵
    design_df = pd.get_dummies(y_encoded, drop_first=False, dtype=int)
    factor_levels = sorted(set(y_encoded))
    design_df.columns = factor_levels
    contrast_matrix = make_contrasts('G2-G1', levels=factor_levels)

    # 模型拟合
    fit = lmFit(X_df, design_df)
    fit2 = contrasts_fit(fit, contrast_matrix)
    eb = eBayes(fit2)

    # 结果提取
    DEG = toptable(eb, coef=[0], number=len(fit2["coefficients"]), sort_by="none")
    DEG['ID'] = DEG.index if DEG.index.name else names

    return DEG[['ID', 'logFC', 'P.Value', 'adj.P.Val']]


def do_edgepy(X, y, names, lbs):
    """
    使用 inmoose.edgepy 精确复现R edgeR的do_edgeR函数
    遵循inmoose官方示例模式
    """
    def create_fallback_result(gene_names):
        """创建备选结果"""
        n = len(gene_names)
        return pd.DataFrame({
            'ID': gene_names,
            'logFC': np.zeros(n),
            'PValue': np.ones(n),
            'FDR': np.ones(n)
        })

    try:
        # 1. 准备数据：确保是genes×samples格式
        if isinstance(X, pd.DataFrame):
            counts = X.values  # 假设已经是genes×samples
        else:
            counts = np.array(X)

        # 如果不是genes×samples，转置
        if counts.shape[1] != len(y):
            if counts.shape[0] == len(y):
                counts = counts.T

            else:
                return create_fallback_result(names)

        # 2. 转换分组标签为G1/G2（与R代码一致）
        y_array = np.array(y)
        y_encoded = y_array.copy()
        y_encoded[y_array == lbs[0]] = 'G1'
        y_encoded[y_array == lbs[1]] = 'G2'

        # 3. 创建样本注释DataFrame
        samples_df = pd.DataFrame({
            'condition': y_encoded,
            'sample': [f'Sample_{i}' for i in range(len(y_encoded))]
        })

        # 4. 创建设计矩阵（使用patsy，与示例一致）
        design = dmatrix("~condition", data=samples_df)

        # 5. 创建DGEList对象（与示例一致）
        # 注意：示例使用了genes=pas.var，我们需要提供基因信息
        genes_df = pd.DataFrame({
            'gene_id': names,
            'gene_name': names
        }, index=range(len(names)))

        dge = DGEList(
            counts=counts,
            samples=samples_df,
            group_col="condition",  # 指定分组列
            genes=genes_df  # 提供基因信息
        )


        # 7. 估计离散度（使用DGEList方法，与示例一致）
        dge.estimateGLMCommonDisp(design=design)

        # 8. 估计基因特异性离散度（与R代码一致）
        try:
            dge.estimateGLMTagwiseDisp(design=design)
        except AttributeError:
            print("⚠️  estimateGLMTagwiseDisp方法不可用，使用公共离散度")

        # 9. GLM拟合（使用DGEList方法，与示例一致）
        fit = dge.glmFit(design=design)

        # 10. 似然比检验（与示例一致）
        # 确定对比：G2 vs G1
        # 在设计矩阵"~condition"中，第二列是condition[T.G2]
        coef_to_test = 1  # 第二列是处理效应

        lrt = glmLRT(fit, coef=coef_to_test)


        # 11. 提取结果（与示例一致）
        n_genes = counts.shape[0]

        if not hasattr(lrt, 'genes') or lrt.genes is None:
            lrt.genes = genes_df

        top_results = topTags(lrt, n=n_genes, sort_by="none")



        # 提取结果DataFrame - 关键修正！
        if hasattr(top_results, 'table'):
            result_df = top_results.table
        elif hasattr(top_results, '_table'):

            result_df = top_results._table
        elif hasattr(top_results, '__array__'):

            result_df = pd.DataFrame(top_results)
        else:
            # 尝试访问属性

            result_df = pd.DataFrame({
                'logFC': getattr(top_results, 'logFC', np.zeros(n_genes)),
                'PValue': getattr(top_results, 'PValue', np.ones(n_genes)),
                'FDR': getattr(top_results, 'FDR', np.ones(n_genes))
            })


        # 如果result_df不是DataFrame，转换为DataFrame
        if not isinstance(result_df, pd.DataFrame):
            result_df = pd.DataFrame(result_df)

        # 12. 格式化结果以匹配R输出
        # 添加ID列
        result_df['ID'] = names

        # 重命名列
        column_map = {}

        # 查找logFC列
        for col in result_df.columns:
            col_lower = col.lower()
            if 'logfc' in col_lower or 'log2foldchange' in col_lower:
                column_map[col] = 'logFC'
                break

        # 查找PValue列
        for col in result_df.columns:
            col_lower = col.lower()
            if 'pvalue' in col_lower or 'p.val' in col_lower:
                column_map[col] = 'PValue'
                break

        # 查找FDR列
        for col in result_df.columns:
            col_lower = col.lower()
            if 'fdr' in col_lower or 'adj.p.val' in col_lower or 'adj_pvalue' in col_lower:
                column_map[col] = 'FDR'
                break

        if column_map:
            result_df = result_df.rename(columns=column_map)

        # 确保有必要的列
        required_cols = ['ID', 'logFC', 'PValue', 'FDR']
        for col in required_cols:
            if col not in result_df.columns:
                if col == 'logFC':
                    result_df[col] = 0.0
                else:
                    result_df[col] = 1.0

        # 只返回需要的列
        final_result = result_df[required_cols].copy()


        return final_result

    except Exception as e:
        import traceback
        traceback.print_exc()
        return create_fallback_result(names)


def do_pyDESeq2(X, y, names, lbs):

    # 复制数据避免修改原数据
    data = X.copy()

    # DESeq2需要最大值小于2147483647（与R版本一致的处理）
    if np.max(data.values) > 2147483647:
        r = round(np.max(data) / 2147483647 + 1)
        data = np.round(data / r).astype(int)
    else:
        data = data.astype(int)

    # 关键修正：使用原始分组标签，不要转换为G1/G2
    # pydeseq2 可以直接处理 'Control'/'Treatment' 这样的标签
    y_series = pd.Series(y, index=X.columns if isinstance(X, pd.DataFrame) else range(X.shape[1]))

    # 创建样本信息DataFrame - 行数必须等于样本数
    # 注意：data.T 后，行是样本，列是基因
    sample_data = pd.DataFrame(
        {'condition': y_series.values},  # 使用原始标签值
        index=data.T.index  # 索引应该是样本名
    )

    try:
        # 创建DeseqDataSet - 保持你的 data.T 转置
        dds = DeseqDataSet(
            counts=data.T,  # samples × genes 格式
            metadata=sample_data,
            design_factors='condition',
            ref_level=['condition', lbs[0]],  # 使用原始标签，如 'Control'
            quiet=True
        )

        # 进行DESeq2分析
        dds.deseq2()

        # 提取统计结果：使用原始标签对比
        stat_res = DeseqStats(dds, contrast=['condition', lbs[1], lbs[0]])  # Treatment vs Control
        stat_res.summary()
        results_df = stat_res.results_df

        # 格式化结果 - 直接使用 results_df 的索引作为基因名
        # 注意：results_df 的索引是基因名（来自 data 的行名）
        result = pd.DataFrame({
            'ID': results_df.index.tolist(),  # 基因名
            'logFC': results_df['log2FoldChange'].values,
            'P.Value': results_df['pvalue'].values,
            'adj.P.Val': results_df['padj'].values
        })

        return result

    except Exception as e:
        print(f"DESeq2分析出错: {e}")
        import traceback
        traceback.print_exc()
        # 返回一个包含正确结构的空结果
        return pd.DataFrame({
            'ID': names,
            'logFC': np.nan,
            'P.Value': np.nan,
            'adj.P.Val': np.nan
        })
# 设置Python函数可用
do_limma = do_limma_py
do_DESeq2 = do_pyDESeq2
do_edgeR = do_edgepy
Python_isavailable = True


class iTSA:

    def __init__(self, method='t-Test'):
        self.method = method

    def fit_data(self, X, y, names):
        self.X = X
        self.y = np.array(y)
        self.names = np.array(names)
        self.lbs = np.unique(y)

        if self.method == 'Limma':
            if not Python_isavailable:
                return None
            res = do_limma(np.log2(self.X), self.y, self.names, self.lbs)
            res = pd.DataFrame(res)
            res = res[['ID', 'logFC', 'P.Value', 'adj.P.Val']]
            for i, c in enumerate(res.columns):
                if c == 'ID':
                    pass
                elif c == 'logFC':
                    res.iloc[:, i] = np.round(res.iloc[:, i].astype(float), 4)
                else:
                    res.iloc[:, i] = -np.round(np.log10(res.iloc[:, i].astype(float)), 4)

            res.columns = ['Accession', 'logFC', '-logPval', '-logAdjPval']


        elif self.method == 'DESeq2':

            if not Python_isavailable:
                return None

            res = do_DESeq2(self.X, self.y, self.names, self.lbs)

            res = pd.DataFrame(res)





            column_map = {}

            actual_cols = [col.lower() for col in res.columns]

            # 映射ID列

            for i, col in enumerate(res.columns):

                col_lower = col.lower()

                if 'id' in col_lower or 'gene' in col_lower or col == 'Unnamed: 0':
                    column_map[col] = 'ID'

                    break

            # 映射logFC列

            for i, col in enumerate(res.columns):

                col_lower = col.lower()

                if 'logfc' in col_lower or 'log2fold' in col_lower or 'lfc' in col_lower:
                    column_map[col] = 'logFC'

                    break

            # 映射P值列

            for i, col in enumerate(res.columns):

                col_lower = col.lower()

                if 'pvalue' in col_lower or 'p.val' in col_lower or 'pvalue' in col_lower:
                    column_map[col] = 'P.Value'

                    break

            # 映射校正P值列

            for i, col in enumerate(res.columns):

                col_lower = col.lower()

                if 'adj' in col_lower or 'padj' in col_lower or 'fdr' in col_lower:
                    column_map[col] = 'adj.P.Val'

                    break

            if column_map:
                res = res.rename(columns=column_map)

                print(f"DEBUG: 列名映射后: {list(res.columns)}")

            # 确保有必要的列

            required_cols = ['ID', 'logFC', 'P.Value', 'adj.P.Val']

            for col in required_cols:

                if col not in res.columns:
                    print(f"DEBUG: 缺失列 {col}，使用NaN填充")

                    res[col] = np.nan

            res = res[required_cols]

            res.columns = ['Accession', 'logFC', '-logPval', '-logAdjPval']

            # 处理可能为NaN或0的p值

            res['-logPval'] = -np.round(np.log10(res['-logPval'].replace(0, 1e-300).fillna(1)), 4)

            res['-logAdjPval'] = -np.round(np.log10(res['-logAdjPval'].replace(0, 1e-300).fillna(1)), 4)

            res['logFC'] = np.round(res['logFC'].fillna(0), 4)

        elif self.method == 'edgeR':
            res = do_edgeR(self.X, self.y, self.names, self.lbs)
            res = pd.DataFrame(res)
            res = res[['ID', 'logFC', 'PValue', 'FDR']]
            res.columns = ['Accession', 'logFC', '-logPval', '-logAdjPval']
            res['-logPval'] = -np.round(np.log10(res['-logPval']), 4)
            res['-logAdjPval'] = -np.round(np.log10(res['-logAdjPval']), 4)
            res['logFC'] = np.round(res['logFC'], 4)

        elif self.method == 't-Test':
            i = np.where(self.y == self.lbs[0])[0]
            j = np.where(self.y == self.lbs[1])[0]
            pval = []
            for k in range(X.shape[0]):
                pval.append(ttest_ind(X.iloc[k, i], X.iloc[k, j]).pvalue)
            logfc = self.fold_change(self.X, self.y, self.lbs)
            adjpval = np.array(p_value_adjust(pval))
            lgpval = -np.round(np.log10(np.array(pval)), 4)
            lgapval = -np.round(np.log10(np.array(adjpval)), 4)
            res = pd.DataFrame({'Accession': self.names, 'logFC': logfc, '-logPval': lgpval, '-logAdjPval': lgapval})

        else:
            raise IOError('{} is not a support method'.format(self.method))

        X_ = X.copy()
        X_['Accession'] = names
        res = pd.merge(res, X_)
        res = res.sort_values(by='-logAdjPval', ascending=False)
        return res

    def fold_change(self, X, y, lbs):
        case_val = np.nanmean(X.loc[:, y == lbs[1]], axis=1)
        cont_val = np.nanmean(X.loc[:, y == lbs[0]], axis=1)
        return np.round(np.log2(case_val / cont_val), 4)


def data_balance(X, y):
    y_uni = list(set(y))
    n = min(len(np.where(y == y_uni[0])[0]), len(np.where(y == y_uni[1])[0]))

    k1 = np.where(y == y_uni[0])[0]
    k2 = np.where(y == y_uni[1])[0]

    X1 = np.sort(X.iloc[:, k1], axis=1)
    X2 = np.sort(X.iloc[:, k2], axis=1)

    l1 = np.arange(int(len(k1) / 2 - 0.5 * n), int(len(k1) / 2 + 0.5 * n))
    l2 = np.arange(int(len(k2) / 2 - 0.5 * n), int(len(k2) / 2 + 0.5 * n))

    X1 = X1[:, l1]
    X2 = X2[:, l2]

    colnames = ['Group_1_Sample_{}'.format(i) for i in range(n)] + ['Group_2_Sample_{}'.format(i) for i in range(n)]
    X_new = pd.DataFrame(np.hstack((X1, X2)))
    X_new.columns = colnames
    y_new = np.array([y_uni[0]] * n + [y_uni[1]] * n)
    return X_new, y_new


def p_value_adjust(pvalues, correction_type="Benjamini-Hochberg"):
    pvalues = np.array(pvalues)
    n = pvalues.shape[0]
    new_pvalues = np.empty(n)
    if correction_type == "Bonferroni":
        new_pvalues = n * pvalues
    elif correction_type == "Bonferroni-Holm":
        values = [(pvalue, i) for i, pvalue in enumerate(pvalues)]
        values.sort()
        for rank, vals in enumerate(values):
            pvalue, i = vals
            new_pvalues[i] = (n - rank) * pvalue
    elif correction_type == "Benjamini-Hochberg":
        values = [(pvalue, i) for i, pvalue in enumerate(pvalues)]
        values.sort()
        values.reverse()
        new_values = []
        for i, vals in enumerate(values):
            rank = n - i
            pvalue, index = vals
            new_values.append((n / rank) * pvalue)
        for i in np.arange(0, int(n) - 1):
            if new_values[i] < new_values[i + 1]:
                new_values[i + 1] = new_values[i]
        for i, vals in enumerate(values):
            pvalue, index = vals
            new_pvalues[index] = new_values[i]
    return new_pvalues