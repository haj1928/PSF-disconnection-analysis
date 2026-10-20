# -*- coding: utf-8 -*-
"""
构建"断连本身"对照向量（区域级, 100 维）——疲劳特异性对照的核心输入
=====================================================================
背景（2026-09-02）：用户方案——用"每患者断连总量"替换 salience 做对照，
      证明分子特征（OXPHOS/小胶质/α4β2 等）是"疲劳维度"特有而非"断连本身"伴随物。
      因 salience 是区域级（100 区），对照向量也必须是区域级（用户方案维度校正）。

方法（Billot et al. 2022 Cortex 的 VDSM vs VLSM 对照思想；区别于疲劳模式的无监督模式）：
  对照向量 1「断连 PC1 载荷」：124×100 断连矩阵无监督 PCA 第一主成分载荷（100 维）。
    代表"断连变异的主空间模式"，完全不含 FSS/疲劳信息（推荐，最干净的对照）。
  对照向量 2「区域平均断连」：跨 124 患者每区断连均值（100 维）。
    代表"断连在哪里发生"的总体分布（简单直观的第二个对照）。

关键检查：corr(对照向量, salience) —— 若 r>0.9 说明断连主模式与疲劳模式几乎重合，
  对照无区分度；预期中等相关（疲劳模式 ≠ 断连平均模式）→ 对照有效。

输入: F:/dcmtonii/PSF_SDC_parcel_discon_N124_cortex100_newmasks_20260915.csv（124×100, 列=区名, 与 salience 对齐）
      F:/dcmtonii/run62_crosssectional/PSF_pls_behavior/PSF_pls_discon_salience_residY_age_sex_cortex100_newmasks_20260915.csv
输出: F:/dcmtonii/run62_crosssectional/PSF_ahba_100/PSF_discon_control_vectors_newmasks_20260915.csv
"""
import numpy as np
import pandas as pd

BASE = 'F:/dcmtonii'
SAL_CSV = BASE + '/run62_crosssectional/PSF_pls_behavior/PSF_pls_discon_salience_residY_age_sex_cortex100_newmasks_20260915.csv'
DC_CSV = BASE + '/PSF_SDC_parcel_discon_N124_cortex100_newmasks_20260915.csv'
OUT_CSV = BASE + '/run62_crosssectional/PSF_ahba_100/PSF_discon_control_vectors_newmasks_20260915.csv'

# ---- 1. 加载 ----
sal = pd.read_csv(SAL_CSV)
dc = pd.read_csv(DC_CSV)
assert set(sal['region']) == set(dc.columns[1:]), 'salience 区与断连矩阵列不一致'
cols = list(dc.columns[1:])            # 100 区, 顺序 = id 1-100
X = dc[cols].values.astype(float)      # 124 × 100
print(f'[1] 断连矩阵 {X.shape} | salience {len(sal)} 区, 顺序一致: {list(sal["region"]) == cols}', flush=True)

# ---- 2. 对照向量 1：PCA 第一主成分载荷 ----
Xc = X - X.mean(0)                     # 每区中心化（PCA 惯例）
U, S, Vt = np.linalg.svd(Xc, full_matrices=False)
pc1_load = Vt[0]                       # 100 维载荷（第一主成分方向）
var_exp = S[0] ** 2 / (S ** 2).sum()
print(f'[2] PCA: PC1 解释方差 {var_exp:.1%}', flush=True)

# ---- 3. 对照向量 2：区域平均断连 ----
mean_discon = X.mean(0)
print(f'[3] 区域平均断连 完成 | 范围 [{mean_discon.min():.4f}, {mean_discon.max():.4f}]', flush=True)

# ---- 4. 与 salience 的相关（区分度检查）----
sal_vals = sal['salience'].values.astype(float)
r_pc1 = np.corrcoef(pc1_load, sal_vals)[0, 1]
r_mean = np.corrcoef(mean_discon, sal_vals)[0, 1]
r_pc1_mean = np.corrcoef(pc1_load, mean_discon)[0, 1]
print(f'[4] corr(PC1载荷, salience) = {r_pc1:+.3f}', flush=True)
print(f'    corr(平均断连, salience) = {r_mean:+.3f}', flush=True)
print(f'    corr(PC1载荷, 平均断连) = {r_pc1_mean:+.3f} (两个对照向量彼此独立度)', flush=True)
if abs(r_pc1) > 0.9 or abs(r_mean) > 0.9:
    print('  ⚠️ 对照与 salience 高度相关（>0.9），对照区分度存疑！', flush=True)
else:
    print('  ✓ 对照区分度 OK（与 salience 中等相关，非共线）', flush=True)

# ---- 5. 保存（region 顺序 = salience 顺序）----
out = pd.DataFrame({'region': sal['region'], 'salience': sal_vals,
                    'discon_pc1': pc1_load, 'discon_mean': mean_discon})
out.to_csv(OUT_CSV, index=False, encoding='utf-8-sig')
print(f'[5] 已保存: {OUT_CSV}', flush=True)
print('\n对照向量前 5 区:')
print(out.head(5).to_string(index=False))
