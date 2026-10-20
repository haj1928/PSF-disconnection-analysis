# -*- coding: utf-8 -*-
"""
细胞类型富集（Seidlitz Dataset5 官方面板版·新主口径）
=====================================================================
背景（2026-09-02）：用户拍板——细胞富集面板改用 Seidlitz Supplementary Dataset 5
      （58 个细胞特异基因集，5 个单细胞研究编译：Lake 2018 / Habib 2018 /
       Li 2018 / Darmanis 2015 / Zhang 2016），合并 7 canonical 类。
      替代原 63 个自选经典 marker 面板（原面板无文献出处注释）。
方法（与主分析 run_celltype_enrich_seidlitz_prand_cortex100_residY.py 完全一致）：
  - 排序：PLS 第一载荷 w1（复用 Step4b 产物，15633 基因全排序）
  - 面板：seidlitz2020_Dataset5_7class_ahba_aligned.json（官方 MOESM8 解析，
    Class ∈ {Astro,Endo,Micro,Neuro-Ex,Neuro-In,Oligo,OPC}，AHBA 命中基因）
  - 富集：7 类基因集中位秩 + PRAND 10,000 次（seed=0）双侧 + 组内 BH-FDR
输入: PSF_pls_gene_ranking_cortex100_residY_newmasks_20260915.csv + 面板 JSON
输出: PSF_ahba_100/PSF_celltype_enrich_seidlitz_Dataset5_newmasks_20260915.csv + report
"""
import numpy as np
import pandas as pd
import json

BASE = 'F:/dcmtonii/run62_crosssectional'
RK_CSV = BASE + '/PSF_ahba_100/PSF_pls_gene_ranking_cortex100_residY_newmasks_20260915.csv'
PANEL_JSON = BASE + '/PSF_references/Seidlitz2020/seidlitz2020_Dataset5_7class_ahba_aligned.json'
OUT_CSV = BASE + '/PSF_ahba_100/PSF_celltype_enrich_seidlitz_Dataset5_newmasks_20260915.csv'
OUT_REPORT = BASE + '/PSF_ahba_100/PSF_celltype_enrich_seidlitz_Dataset5_newmasks_20260915_report.txt'
N_PERM = 10000
SEED = 0

# ---- 1. 读排序 + 官方面板 ----
rk = pd.read_csv(RK_CSV)
genes = rk['gene'].astype(str).tolist()
rank_of = rk['rank'].values.astype(int)
print(f'[1] 排序加载: {len(genes)} 基因', flush=True)

with open(PANEL_JSON, encoding='utf-8') as f:
    MARKERS = json.load(f)
for ct, gl in MARKERS.items():
    print(f'  面板 {ct:<22} {len(gl)} 基因(AHBA命中)', flush=True)

# ---- 2. 中位秩 + PRAND ----
rng = np.random.default_rng(SEED)
g2i = {g: i for i, g in enumerate(genes)}
rows = []
center = np.median(rank_of)
print(f'[2] 全基因组中位秩 = {center:.0f}', flush=True)

for ct, gl in MARKERS.items():
    idx = [g2i[g] for g in gl if g in g2i]
    if not idx:
        print(f'  [跳过] {ct}: marker 全不在基因库'); continue
    k = len(idx)
    obs_rank = np.median(rank_of[idx])
    null_idx = rng.integers(0, len(genes), size=(N_PERM, k))
    null_rank = np.median(rank_of[null_idx], axis=1)
    dist_obs = abs(obs_rank - center)
    dist_null = np.abs(null_rank - center)
    p_twotail = (1 + np.sum(dist_null >= dist_obs)) / (1 + N_PERM)
    direction = '前端(正共变)' if obs_rank < center else '后端(负共变)'
    rows.append({'cell_type': ct, 'n_marker': k, 'median_rank': int(obs_rank),
                 'null_median': float(np.median(null_rank)),
                 'null_95CI': f"[{np.percentile(null_rank,2.5):.0f},{np.percentile(null_rank,97.5):.0f}]",
                 'p_prand_twotail': p_twotail, 'direction': direction})
    print(f'  {ct:<22} n={k:>4} 中位秩={obs_rank:>7.0f} null={np.median(null_rank):>6.0f} p={p_twotail:.4f} {direction}', flush=True)

# ---- 3. 组内 BH-FDR ----
df = pd.DataFrame(rows)
pvals = df['p_prand_twotail'].values
n = len(pvals)
order_p = np.argsort(pvals)
fdr = pvals[order_p] * n / np.arange(1, n + 1)
fdr = np.minimum.accumulate(fdr[::-1])[::-1]
fdr_out = np.empty(n); fdr_out[order_p] = np.clip(fdr, 0, 1)
df['fdr_bh'] = fdr_out
df = df.sort_values('p_prand_twotail').reset_index(drop=True)
df.to_csv(OUT_CSV, index=False, encoding='utf-8-sig')

rep = [
    '=' * 78,
    '细胞类型富集（Seidlitz Dataset5 官方面板·新主口径）：中位秩 + PRAND',
    '表型 Y = salience_residY | X = AHBA 100 区表达排序（15633 基因, Step4b 产物）',
    '面板: Seidlitz 2020 Supplementary Dataset 5（58 基因集, 5 单细胞研究 Lake/Habib/Li/Darmanis/Zhang）→ 7 canonical 类',
    f'PRAND: {N_PERM} 次基因置换 seed={SEED} | 双侧 | BH-FDR: 7 类组内',
    '-' * 78,
    f"{'细胞类型':<22}{'n':>5}{'中位秩':>9}{'null中位':>9}{'p(PRAND)':>10}  方向",
]
for _, r in df.iterrows():
    rep.append(f"{r['cell_type']:<22}{int(r['n_marker']):>5}{r['median_rank']:>9.0f}{r['null_median']:>9.0f}{r['p_prand_twotail']:>10.4f}  {r['direction']}")
rep.append('-' * 78)
rep.append('组内 BH-FDR 校正后：')
for _, r in df.iterrows():
    mark = ' ***' if r['fdr_bh'] < 0.05 else (' *' if r['p_prand_twotail'] < 0.05 else '')
    rep.append(f"  {r['cell_type']:<22} p={r['p_prand_twotail']:.4f} FDR={r['fdr_bh']:.4f}{mark}")
rep.append('')
rep.append('对照: 63 自选 marker 面板版（旧主口径）Micro p=0.0003 FDR=0.0021')
rep.append('数据源: Seidlitz et al. 2020 Nat Commun 11:3358 Supplementary Dataset 5 (MOESM8)')
rep.append('        作者 GitHub 同源文件 celltypes_PSP.csv (Jaccard vs MOESM8 7类合并 = 0.995-0.998)')
rep.append('        参考链: Chu et al. 2025 Biol Psychiatry 97:290-301 (ref 42 = Seidlitz 2020)')
with open(OUT_REPORT, 'w', encoding='utf-8') as f:
    f.write('\n'.join(rep))
print('\n' + '\n'.join(rep[:12]))
print(f'\n[3] 已保存:\n  {OUT_CSV}\n  {OUT_REPORT}')
