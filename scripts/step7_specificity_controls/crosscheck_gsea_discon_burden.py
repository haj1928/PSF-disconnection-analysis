# -*- coding: utf-8 -*-
"""
disconnection burden 作排序向量的 GSEA 官方口径对照（gseapy prerank）
=====================================================================
背景（2026-09-03, 用户拍板 GSEA 主口径 = gseapy 官方）：Step7 特异性对照
  表 B 原用自实现（Tian mean-ES）跑 burden 版 0 存活；现统一为官方口径
  （Subramanian 2005 加权 KS, gseapy.prerank, 1000 置换, seed=0, 与
  `run_gsea_gseapy_w1_cortex100_residY.py` 完全同参数），仅把排序向量
  从 salience 换成 disconnection burden（区域级断连总量, 不含 FSS）。
  度量（2026-09-04 更新）：与主口径一致用 Pearson（w1∝Pearson, 单响应 PLS
  载荷 ∝ Xᵀy），保证"主 vs 对照"同源可比。

输出: PSF_ahba_100/gseapy_prerank_discon_burden_newmasks_20260915_{GO,KEGG}.csv
"""
import numpy as np
import pandas as pd
import gseapy as gp
from scipy.stats import rankdata

BASE = 'F:/dcmtonii/run62_crosssectional'
GMT = BASE + '/PSF_references/gene_sets'
OUT = BASE + '/PSF_ahba_100'
N_PERM = 1000
SEED = 0

# ---- 1. burden 向量（= discon_mean 列, 官方对照表型）----
ctrl = pd.read_csv(OUT + '/PSF_discon_control_vectors_newmasks_20260915.csv')
expr = pd.read_csv(OUT + '/PSF_ahba_100_expression.csv', index_col=0)
burden_map = dict(zip(ctrl['region'], ctrl['discon_mean']))
expr_valid = ~expr.isna().all(axis=1)
common = [r for r in expr.index if r in burden_map and expr_valid[r]]
common_expr = expr.loc[common].astype(float)
common_burden = np.array([burden_map[r] for r in common])
print(f'有效区 {len(common)} | burden 范围 [{common_burden.min():.4f}, {common_burden.max():.4f}]')

# ---- 2. r_gene = Pearson(burden, 表达)（2026-09-04 统一输入: 与主口径 w1∝Pearson 同源）----
ce = common_expr.values.astype(float)                        # (R, G)
ce_z = (ce - ce.mean(0)) / (ce.std(0, ddof=1) + 1e-12)       # 每基因跨区 Z-score
bz = (common_burden - common_burden.mean()) / (common_burden.std(ddof=1) + 1e-12)
r_gene = (ce_z.T @ bz) / (len(common) - 1)                   # (G,)
genes = common_expr.columns.astype(str)
rnk = pd.DataFrame({'gene': genes, 'score': r_gene}).dropna()
print(f'r_gene (Pearson) 就绪, n={len(rnk)} | 并列值 {int(rnk["score"].duplicated().sum())}')

# ---- 3. gseapy prerank（与 salience 主分析完全同参数）----
for name, fn in [('GO', 'GO_Biological_Process_2021.gmt'), ('KEGG', 'KEGG_2021_Human.gmt')]:
    pre = gp.prerank(rnk=rnk, gene_sets=f'{GMT}/{fn}', outdir=None,
                     min_size=5, max_size=300,
                     permutation_num=N_PERM, seed=SEED, threads=4,
                     permutation_type='gene_set', no_plot=True,
                     verbose=False)
    res = pre.res2d.reset_index(drop=True)
    sig = (res['FDR q-val'] < 0.05).sum()
    print(f'[{name}] 通路 {len(res)} | FDR<0.05: {sig}')
    res.to_csv(f'{OUT}/gseapy_prerank_discon_burden_newmasks_20260915_{name}.csv', index=False,
               encoding='utf-8-sig')

print('完成。burden 官方口径对照输出见', OUT)
