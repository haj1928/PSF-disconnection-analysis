# -*- coding: utf-8 -*-
"""
GSEA 官方口径（gseapy 加权 KS）显著通路的 salience 表型置换（spin 1000×）
=========================================================================
背景（2026-09-03, 用户拍板）：GSEA 主口径统一为 gseapy 官方（Subramanian 2005
  加权 KS, `gseapy_prerank_cortex100_residY_{GO,KEGG}.csv`）；自实现 Tian mean-ES
  版不再作为报告口径。本脚本对官方 FDR<0.05 的全部显著通路做 salience 空间置换
  （spin），检验其 ES 是否显著超出"保留空间结构的随机 salience"（表型层 null,
  Alexander-Bloch 2018 / Seidlitz 2020 spin 框架）。

度量说明（2026-09-04 统一输入更新）：GSEA 主口径输入已由 Spearman r_gene 切换为
  Step3 的 PLS 载荷 w1（= 每基因 Pearson(salience, 表达)，单响应 PLS 载荷 ∝ Xᵀy,
  实测 corr(w1, Pearson)=1.000000）。本脚本同步改为 Pearson 度量以与 gseapy 官方
  ES 同源核对（ES 对 score 全局缩放不变，w1 与 Pearson 排序一致 → ES 完全复现）。

实现：
  1. r_gene = Pearson(salience, AHBA 表达)（表达每基因跨 99 区 Z-score）与主口径同源
  2. spin null：1000 个球面旋转 salience（Schaefer100 fsa5, seed=0,
     文件 `PSF_sal_perm_spin.csv` 已由 build_salience_permutations.py 生成）
     → 每个重算 r_gene → 对每条显著通路算加权 KS ES
  3. 加权 KS ES 复现 Subramanian 2005（p=1）：排序基因列表按 r_gene 降序，
     hit 增量 |r|/Σ|r|_S、miss 增量 1/(N-NH)，ES=游走最大偏差（带符号）
  4. 验证：观测 ES 与 gseapy 官方输出 ES 列一致（容差 1e-4）
  5. p_spin = (1 + Σ|null|≥|obs|) / (1 + 1000) 双尾（同受体层 spin 判读）

输出: PSF_ahba_100/permutation_newmasks_20260915/PSF_gsea_pathway_spin_official_result.csv
      + report（存活数、代谢核心通路明细）
"""
import numpy as np
import pandas as pd
import os, time, re
from scipy.stats import rankdata

BASE = 'F:/dcmtonii/run62_crosssectional'
OUT = BASE + '/PSF_ahba_100/permutation_newmasks_20260915'
os.makedirs(OUT, exist_ok=True)
GMT = BASE + '/PSF_references/gene_sets'

# ---------------- 1. 输入 ----------------
sal_df = pd.read_csv(BASE + '/PSF_pls_behavior/PSF_pls_discon_salience_residY_age_sex_cortex100_newmasks_20260915.csv')
parcels = sal_df['region'].tolist()
sal_100 = sal_df['salience'].values.astype(float)

expr = pd.read_csv(BASE + '/PSF_ahba_100/PSF_ahba_100_expression.csv', index_col=0)
keep = [i for i, r in enumerate(parcels) if not expr.loc[r].isna().all()]
print(f'[1] 有效皮层区 {len(keep)} (剔除 RH_SomMot_2)', flush=True)
expr_v = expr.iloc[keep].astype(float)          # 99 × 15633
genes = np.array([str(g) for g in expr_v.columns])
sal99 = sal_100[keep]                            # (99,)

# spin 置换 salience（1000×100, 第0行=表头）
spin = pd.read_csv(OUT + '/PSF_sal_perm_spin.csv')
spin_sal = spin.values.astype(float)             # (1000, 100)
spin99 = spin_sal[:, keep]                       # (1000, 99)
sal_all = np.vstack([sal99, spin99])             # (1001, 99)  row0=观测
print(f'[2] spin 置换: {spin99.shape[0]} 个 | 观测+null = {sal_all.shape[0]} 行', flush=True)

# ---------------- 2. r_gene（Pearson 99 区, 与主口径 w1∝Pearson 同源）向量化 ----------------
def zscore_axis0(M):
    """按列(axis0, 每基因跨区) Z-score: M shape (R, G)"""
    mu = M.mean(0, keepdims=True); sd = M.std(0, ddof=1, keepdims=True) + 1e-12
    return (M - mu) / sd

def zscore_rows(M):
    """逐行(每旋转的 99 区) Z-score: M shape (K, R)"""
    return (M - M.mean(1, keepdims=True)) / (M.std(1, ddof=1, keepdims=True) + 1e-12)

# 表达侧 Z-score（固定, 99 区 per gene）
z_expr = zscore_axis0(expr_v.values.astype(float))          # (99, G)
z_sal = zscore_rows(sal_all)                                # (1001, 99)
r_gene_all = (z_sal @ z_expr) / (len(keep) - 1)             # (1001, G) row0=观测
print(f'[3] r_gene (Pearson): {r_gene_all.shape}', flush=True)

# ---------------- 3. 加载 gseapy 官方显著通路 ----------------
sig_paths = []
for gmt, gname in [('GO', 'GO_Biological_Process_2021'), ('KEGG', 'KEGG_2021_Human')]:
    res = pd.read_csv(BASE + f'/PSF_ahba_100/gseapy_prerank_cortex100_residY_w1_newmasks_20260915_{gmt}.csv')
    sig = res[res['FDR q-val'] < 0.05]
    sets = {}
    with open(f'{GMT}/{gname}.gmt', encoding='utf-8', errors='ignore') as f:
        for line in f:
            parts = line.rstrip('\n').split('\t')
            if len(parts) >= 3:
                sets[parts[0]] = parts[2:]
    for _, r in sig.iterrows():
        term = r['Term']
        sgenes = sets.get(term)
        if sgenes is None:
            # 名称匹配（Term 可能带括号描述, GMT 主名在 Name 列）
            alt = sets.get(r['Name'])
            sgenes = alt if alt is not None else None
        if sgenes is None:
            print(f'  !! 未匹配 GMT: {term[:60]}')
            continue
        sig_paths.append({'gmt': gmt, 'term': term, 'es_obs_gseapy': r['ES'],
                          'nes': r['NES'], 'fdr': r['FDR q-val'], 'genes': sgenes})
print(f'[4] gseapy 官方显著通路: {len(sig_paths)} 条', flush=True)

# ---------------- 4. 加权 KS ES（Subramanian 2005, p=1）向量化 ----------------
K, G = r_gene_all.shape
gidx = np.arange(G)[None, :]
abs_r = np.abs(r_gene_all)
order = np.argsort(-r_gene_all, axis=1)                    # (1001, G) 降序
sorted_gidx = np.take_along_axis(np.broadcast_to(gidx, (K, G)), order, axis=1)
sorted_abs = np.take_along_axis(abs_r, order, axis=1)

def weighted_ks_es(hit_mask, nr_denom):
    """hit_mask: (K,G) bool（排序空间命中）; 返回 (K,) 带符号 ES"""
    nh = int(hit_mask[0].sum())
    hit_add = np.where(hit_mask, sorted_abs / np.maximum(nr_denom[:, None], 1e-12), 0.0)
    miss_add = np.where(hit_mask, 0.0, 1.0 / (G - nh))
    diff = np.cumsum(hit_add - miss_add, axis=1)
    posmax = diff.max(1)
    negmin = diff.min(1)
    return np.where(np.abs(posmax) >= np.abs(negmin), posmax, negmin)

rows = []
t0 = time.time()
n_miss = 0
for i, p in enumerate(sig_paths):
    idx = np.array([j for j, g in enumerate(genes) if g in p['genes']], dtype=int)
    if not (5 <= len(idx) <= 300) or len(idx) == 0:
        n_miss += 1
        continue
    in_set = np.isin(sorted_gidx, idx)                     # (K,G)
    nr = (sorted_abs * in_set).sum(1)
    es_all = weighted_ks_es(in_set, nr)                    # (K,)  [0]=观测
    es_obs = es_all[0]
    es_null = es_all[1:]
    p_spin = (1 + np.sum(np.abs(es_null) >= abs(es_obs))) / (K)
    # 2026-10-09：新增同向单侧（论文定稿口径；报告 44 §0）
    p_spin_one = ((1 + np.sum(es_null >= es_obs)) / K) if es_obs >= 0 else ((1 + np.sum(es_null <= es_obs)) / K)
    rows.append({'gmt': p['gmt'], 'pathway': p['term'], 'ES_obs': es_obs,
                 'p_spin_one': p_spin_one,
                 'ES_obs_gseapy': p['es_obs_gseapy'], 'NES': p['nes'],
                 'fdr': p['fdr'], 'ES_null_median': float(np.median(es_null)),
                 'p_spin': p_spin})
    if (i + 1) % 30 == 0:
        print(f'  [{i+1}/{len(sig_paths)}] {time.time()-t0:.0f}s', flush=True)

res = pd.DataFrame(rows)
# 验证观测 ES 与 gseapy 一致性
res['es_diff'] = np.abs(res['ES_obs'] - res['ES_obs_gseapy'])
print(f'\n[5] 完成 {len(res)} 条 | ES 与 gseapy max|diff| = {res["es_diff"].max():.6f} '
      f'(中位 {res["es_diff"].median():.2e})', flush=True)

# ---------------- 5. 输出 ----------------
res = res.sort_values('p_spin').reset_index(drop=True)
res.to_csv(OUT + '/PSF_gsea_pathway_spin_official_result.csv', index=False,
           encoding='utf-8-sig')

surv = res[res['p_spin'] < 0.05]
n_one = int((res['p_spin_one'] < 0.05).sum())
n_go = int((res['gmt'] == 'GO').sum()); n_kegg = int((res['gmt'] == 'KEGG').sum())
with open(OUT + '/PSF_gsea_pathway_spin_official_report.txt', 'w', encoding='utf-8') as f:
    f.write(f'GSEA 官方口径通路层表型置换（spin 1000×, seed=0, Pearson/w1 度量 2026-09-04；行对齐质心+同向单侧 2026-10-09）\n')
    f.write(f'显著通路: GO {n_go} + KEGG {n_kegg} = {len(res)} | 有效检验: {len(res)}\n')
    f.write(f'p_spin_one<0.05 存活（同向单侧·稿件口径）: {n_one}/{len(res)}\n')
    f.write(f'p_spin<0.05 存活（双侧·对照）: {len(surv)}/{len(res)}\n\n')
    f.write('=== 未存活 (p_spin_one>=0.05, 同向单侧) ===\n')
    for _, r in res[res['p_spin_one'] >= 0.05].iterrows():
        f.write(f"  {r['pathway'][:70]} | p_spin_one={r['p_spin_one']:.4f} FDR={r['fdr']:.4f}\n")
    f.write('\n=== 代谢核心通路 (top) ===\n')
    met_kw = 'oxidat|glycoly|electron transport|respir|ATP synth|mitochondr|aerobic|glucose|proteasome|ribosome|translat|aminoacyl'
    met = res[res['pathway'].str.contains(met_kw, case=False, na=False)]
    for _, r in met.head(15).iterrows():
        f.write(f"  {r['pathway'][:70]} | ES={r['ES_obs']:+.4f} p_spin={r['p_spin']:.4f} FDR={r['fdr']:.4f}\n")
    f.write(f'\n代谢类通路 {len(met)} 条 | p_spin<0.05: {(met["p_spin"]<0.05).sum()}\n')
print('完成。输出:', OUT)
