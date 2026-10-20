# -*- coding: utf-8 -*-
"""
Step 7'（皮层 100 区·主分析版）：断连 salience 的 GSEA（新口径复核）
=====================================================================
背景（2026-09-01）：主分析定稿为"仅残差 Y(FSS9) 对[年龄,性别]"皮层 100 区，
      产生新 salience（PSF_pls_discon_salience_residY_age_sex_cortex100_newmasks_20260915.csv）。
      原 Step7 GSEA（报告 §7.6）用的是全脑 125 区旧口径，免疫信号
      （IFN-γ/抗原呈递/小胶质激活/突触修剪）需在新口径复核是否保留，
      以决定"免疫+代谢双轴"叙事是否成立。

方法（与旧版 run_gsea_discon.py 完全一致，保证可比）：
  - 每基因 r_gene = Spearman(断连 salience, 基因表达) 跨 99 皮层区 (15633,)
  - ES = 基因集内 r_gene 的均值（Tian et al. 2005 competitive test）
  - null = 基因置换 10,000 次, seed=0, 双尾 |ES|>=|obs|
  - BH-FDR（GO 与 KEGG 分别校正）；基因集大小过滤 5-300
基因集: Enrichr GO_Biological_Process_2021 (6035) / KEGG_2021_Human (320)
输入: PSF_pls_discon_salience_residY_age_sex_cortex100_newmasks_20260915.csv (100 区)
      PSF_ahba_100/PSF_ahba_100_expression.csv (100×15633, RH_SomMot_2 全 NaN 剔除)
输出: PSF_ahba_100/PSF_discon_GSEA_cortex100_residY_newmasks_20260915_{GO,KEGG}.csv + report
"""
import numpy as np, pandas as pd, os, time
from scipy.stats import rankdata

BASE = 'F:/dcmtonii/run62_crosssectional'
SAL_CSV = BASE + '/PSF_pls_behavior/PSF_pls_discon_salience_residY_age_sex_cortex100_newmasks_20260915.csv'
EXPR_CSV = BASE + '/PSF_ahba_100/PSF_ahba_100_expression.csv'
GMT = BASE + '/PSF_references/gene_sets'
OUT = BASE + '/PSF_ahba_100'
N_PERM = 10000
SEED = 0
MIN_N, MAX_N = 5, 300

def bh_fdr(pvals):
    p = np.asarray(pvals, float); n = len(p)
    order = np.argsort(p); ranked = p[order]
    fdr = ranked * n / (np.arange(1, n + 1))
    fdr = np.minimum.accumulate(fdr[::-1])[::-1]
    out = np.empty(n); out[order] = np.clip(fdr, 0, 1)
    return out

def load_gmt(fn):
    sets = {}
    with open(fn, encoding='utf-8', errors='ignore') as f:
        for line in f:
            parts = line.rstrip('\n').split('\t')
            if len(parts) >= 3:
                sets[parts[0]] = parts[2:]
    return sets

# ---- 1. 对齐：新 salience(100区) ∩ 表达(100区, 剔全NaN) → 99 区 ----
sal = pd.read_csv(SAL_CSV)
expr = pd.read_csv(EXPR_CSV, index_col=0)
sal_map = dict(zip(sal['region'], sal['salience']))
expr_valid = ~expr.isna().all(axis=1)
common = [r for r in expr.index if r in sal_map and expr_valid[r]]
common_expr = expr.loc[common].astype(float)
common_sal = np.array([sal_map[r] for r in common])
genes = np.array([str(g) for g in expr.columns])
print(f'[1] 对齐区: {len(common)}（100 - RH_SomMot_2 全NaN）| salience {common_sal.shape} expr {common_expr.shape}', flush=True)

# ---- 2. r_gene = Spearman(秩相关) ----
rs = rankdata(common_sal); rs_z = (rs - rs.mean()) / (rs.std(ddof=1) + 1e-12)
re = np.apply_along_axis(rankdata, 0, common_expr.values)
re_z = (re - re.mean(0)) / (re.std(0, ddof=1) + 1e-12)
r_gene = (rs_z @ re_z) / (len(common) - 1)
print(f'[2] r_gene 完成, n={len(r_gene)}, 范围 [{r_gene.min():.3f}, {r_gene.max():.3f}]', flush=True)

rng = np.random.default_rng(SEED)
report = []
for gmt_name in ['GO', 'KEGG']:
    fn = f'{GMT}/{"GO_Biological_Process_2021" if gmt_name=="GO" else "KEGG_2021_Human"}.gmt'
    sets = load_gmt(fn)
    rows = []
    t0 = time.time()
    for k, (sname, sgenes) in enumerate(sets.items()):
        idx = np.array([i for i, g in enumerate(genes) if g in sgenes])
        if not (MIN_N <= len(idx) <= MAX_N):
            continue
        obs = r_gene[idx].mean()
        sel = rng.integers(0, len(genes), size=(N_PERM, len(idx)))
        null = r_gene[sel].mean(1)
        p = (1 + np.sum(np.abs(null) >= abs(obs))) / (1 + N_PERM)
        rows.append(dict(pathway=sname, n=len(idx), ES=float(obs), perm_p=float(p)))
    df = pd.DataFrame(rows)
    df['fdr'] = bh_fdr(df['perm_p'].values)
    df = df.sort_values('perm_p')
    df.to_csv(f'{OUT}/PSF_discon_GSEA_cortex100_residY_newmasks_20260915_{gmt_name}.csv', index=False, encoding='utf-8-sig')
    n_sig = int((df['fdr'] < 0.05).sum())
    n_nom = int((df['perm_p'] < 0.05).sum())
    print(f'[{gmt_name}] {len(df)} 通路(过滤后) | 名义 p<0.05: {n_nom} | FDR<0.05: {n_sig} | {time.time()-t0:.0f}s', flush=True)
    print(f'  top FDR 通路:', flush=True)
    for _, r in df[df['fdr'] < 0.05].head(20).iterrows():
        print(f'    {r["pathway"][:70]:<72s} ES={r["ES"]:+.4f} p={r["perm_p"]:.4f} FDR={r["fdr"]:.4f}', flush=True)
    report.append((gmt_name, len(df), n_nom, n_sig))

with open(f'{OUT}/PSF_discon_GSEA_cortex100_residY_newmasks_20260915_report.txt', 'w', encoding='utf-8') as f:
    f.write('Step7\' 断连 salience GSEA 新口径复核（皮层100区, 仅残差Y salience）\n')
    f.write('方法: 每基因 Spearman(salience, 表达) r_gene; ES=基因集内均值; 基因置换10k seed=0 双尾; BH-FDR\n')
    f.write(f'输入: salience(residY, {len(common)}皮层区) × AHBA 15633 基因\n')
    f.write(f'基因集: Enrichr GO_Biological_Process_2021 / KEGG_2021_Human（命中 5-300 过滤）\n\n')
    for name, n, n_nom, n_sig in report:
        f.write(f'{name}: {n} 通路 | 名义 p<0.05: {n_nom} | FDR<0.05: {n_sig}\n')
    f.write('\nES 为正 = 通路基因与断连 salience 正相关（断连-疲劳相关区域富集该通路）\n')
    f.write('ES 为负 = 负相关\n')
print('完成。', flush=True)
