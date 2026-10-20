# -*- coding: utf-8 -*-
"""
细胞类型富集 · HPA 脑表达过滤敏感性（Dataset5 官方面板版, 2026-09-03 固化）
=====================================================================
背景：审稿 M2（脑表达过滤）敏感性在 Dataset5 面板下此前为命令行执行，
      本脚本将其固化为正式产物（逻辑与 `run_celltype_enrich_hpafilter.py`
      完全一致，仅面板换成 Seidlitz Dataset5 官方 7 类）。
方法：
  - 排序：PLS 第一载荷 w1（复用 Step4b 产物, 15633 基因全排序）
  - 过滤：HPA consensus normal tissue, cerebral cortex nTPM>=1 = detected
    （HPA 官方阈值; 掩码文件 hpa_cortex_detected_genes_consensus.csv, 12,753）
  - 宇宙：原 15633 ∩ HPA detected → 按原 w1 降序重排 rank 1..N
  - 富集：Dataset5 7 类基因集中位秩 + PRAND 10,000 次（seed=0）双侧 + 组内 BH-FDR
输入: PSF_pls_gene_ranking_cortex100_residY_newmasks_20260915.csv
      hpa_cortex_detected_genes_consensus.csv（资源_面板与掩码/）
      seidlitz2020_Dataset5_7class_ahba_aligned.json（资源_面板与掩码/）
输出: PSF_ahba_100/PSF_celltype_enrich_seidlitz_Dataset5_hpafilter_newmasks_20260915.csv + report
"""
import numpy as np
import pandas as pd
import json

BASE = 'F:/dcmtonii/run62_crosssectional'
OUT = BASE + '/PSF_ahba_100'
RK_CSV = OUT + '/PSF_pls_gene_ranking_cortex100_residY_newmasks_20260915.csv'
FLT_CSV = OUT + '/hpa_cortex_detected_genes_consensus.csv'
PANEL_JSON = BASE + '/PSF_references/Seidlitz2020/seidlitz2020_Dataset5_7class_ahba_aligned.json'
N_PERM = 10000
SEED = 0

# ---- 1. 排序 + HPA 过滤宇宙 ----
rk = pd.read_csv(RK_CSV)
flt = set(pd.read_csv(FLT_CSV)['gene'].astype(str))
keep = rk['gene'].astype(str).isin(flt)
rk_f = rk[keep].sort_values('loading', ascending=False).reset_index(drop=True)
rk_f['rank'] = np.arange(1, len(rk_f) + 1)
genes = rk_f['gene'].astype(str).tolist()
rank_of = rk_f['rank'].values.astype(int)
g2i = {g: i for i, g in enumerate(genes)}
print(f'HPA 过滤后基因宇宙: {len(genes)}（原 15633, 掩码 12,753）', flush=True)

# ---- 2. Dataset5 官方面板 ----
PANEL = json.load(open(PANEL_JSON, encoding='utf-8'))
print(f'面板 7 类: { {k: len(v) for k, v in PANEL.items()} }', flush=True)

# ---- 3. PRAND（与主分析版逐字一致）----
rng = np.random.default_rng(SEED)
center = np.median(rank_of)
rows = []
for ct, gl in PANEL.items():
    idx = [g2i[g] for g in gl if g in g2i]
    k = len(idx)
    if k == 0:
        continue
    obs = np.median(rank_of[idx])
    null_idx = rng.integers(0, len(genes), size=(N_PERM, k))
    null_rank = np.median(rank_of[null_idx], axis=1)
    dist_obs = abs(obs - center)
    dist_null = np.abs(null_rank - center)
    p_twotail = (1 + np.sum(dist_null >= dist_obs)) / (1 + N_PERM)
    direction = '前端(正共变)' if obs < center else '后端(负共变)'
    rows.append({'cell_type': ct, 'n_marker': k, 'median_rank': int(obs),
                 'null_median': float(np.median(null_rank)), 'p_prand_twotail': p_twotail,
                 'direction': direction})
    print(f'  {ct:<22} n={k:>4} 中位秩={obs:>7.0f} p={p_twotail:.4f} {direction}', flush=True)

df = pd.DataFrame(rows)
# BH-FDR（7 类组内）
pvals = df['p_prand_twotail'].values
n = len(pvals)
order_p = np.argsort(pvals)
fdr = pvals[order_p] * n / np.arange(1, n + 1)
fdr = np.minimum.accumulate(fdr[::-1])[::-1]
fdr_out = np.empty(n)
fdr_out[order_p] = np.clip(fdr, 0, 1)
df['fdr_bh'] = fdr_out
df = df.sort_values('p_prand_twotail').reset_index(drop=True)

# ---- 4. 输出 ----
df.to_csv(OUT + '/PSF_celltype_enrich_seidlitz_Dataset5_hpafilter_newmasks_20260915.csv', index=False,
          encoding='utf-8-sig')
with open(OUT + '/PSF_celltype_enrich_seidlitz_Dataset5_hpafilter_report_newmasks_20260915.txt', 'w',
          encoding='utf-8') as f:
    f.write('Dataset5 官方面板 · HPA 脑表达过滤敏感性（PRAND 10k, seed=0, 双侧, 组内 BH-FDR）\n')
    f.write(f'HPA 过滤宇宙: {len(genes)} 基因\n\n')
    for _, r in df.iterrows():
        mark = ' ***' if r['fdr_bh'] < 0.05 else ''
        f.write(f"  {r['cell_type']:<22} n={int(r['n_marker']):>5} "
                f"中位秩={int(r['median_rank']):>7} p={r['p_prand_twotail']:.4f} "
                f"FDR={r['fdr_bh']:.4f}{mark}  {r['direction']}\n")
print('完成。输出:', OUT + '/PSF_celltype_enrich_seidlitz_Dataset5_hpafilter_newmasks_20260915.csv')
