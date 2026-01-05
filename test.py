# debug_limma.py
import pandas as pd
import numpy as np
import sys

sys.path.append('.')

print("🔍 调试 limma_py 维度问题")
print("=" * 60)

# 创建简单测试数据
np.random.seed(123)
X = np.random.rand(10, 6) * 1000  # 10基因×6样本
names = [f"Gene_{i}" for i in range(10)]
y = ['Control', 'Control', 'Control', 'Treatment', 'Treatment', 'Treatment']
lbs = ['Control', 'Treatment']

X_df = pd.DataFrame(X, index=names)

print(f"数据形状: {X_df.shape}")
print(f"分组: {y}")

# 测试不同的设计矩阵
from patsy import dmatrix

# 1. 无截距设计矩阵（与R代码一致）
y_encoded = np.array(y)
y_encoded[y_encoded == 'Control'] = 'G1'
y_encoded[y_encoded == 'Treatment'] = 'G2'

design_df = pd.DataFrame({'group': y_encoded})

print("\n1. 无截距设计矩阵 (~0+group):")
design_no_intercept = dmatrix('~0+group', data=design_df)
print(f"   形状: {design_no_intercept.shape}")
print(f"   列名: {design_no_intercept.design_info.column_names}")

print("\n2. 有截距设计矩阵 (~group):")
design_with_intercept = dmatrix('~group', data=design_df)
print(f"   形状: {design_with_intercept.shape}")
print(f"   列名: {design_with_intercept.design_info.column_names}")

# 测试lmFit
try:
    from limma_py import lmFit

    print("\n3. 测试lmFit:")

    # 测试无截距
    print("   无截距设计矩阵:")
    fit1 = lmFit(X_df, design_no_intercept)
    print(f"     系数形状: {fit1.coefficients.shape}")

    # 测试有截距
    print("   有截距设计矩阵:")
    fit2 = lmFit(X_df, design_with_intercept)
    print(f"     系数形状: {fit2.coefficients.shape}")

except Exception as e:
    print(f"   错误: {e}")