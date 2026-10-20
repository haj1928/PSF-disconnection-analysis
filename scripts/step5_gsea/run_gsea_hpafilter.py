# -*- coding: utf-8 -*-
"""
GSEA 通路富集（HPA 脑表达过滤敏感性版）：审稿 M2 回应
=====================================================================
背景（2026-09-02）：Seidlitz 2020 在富集前用 HPA 剔除 cortex not-detected 基因
      （15043→7971）。本脚本复刻该过滤（HPA consensus, cortex nTPM>=1 detected），
      在过滤后基因宇宙（12753）内重跑 GSEA 官方链路，检验代谢/免疫通路是否
      仍显著（排除"非脑表达/低可靠基因稀释 null"的偏倚）。
方法：与主分析 run_gsea_discon_cortex100_residY.py 完全同参数，仅基因宇宙不同：
  - r_gene = Spearman(salience, expr) 对过滤后基因
  - ES = 基因集内 r_gene 均值；null = 基因置换 10,000 次 seed=0 双尾
  - BH-FDR（GO 与 KEGG 分别）；基因集大小过滤 5-300
输入: salience_residY(99区) × expr(15633) × hpa_cortex_detected_genes_consensus.csv
输出: PSF_ahba_100/PSF_discon_GSEA_hpafilter_newmasks_20260915_{GO,KEGG}.csv + report
"""
import numpy as np, pandas as pd, os, time
from scipy.stats import rankdata

BASE = 'F:/dcmtonii/run62_crosssectional'
SAL_CSV = BASE + '/PSF_pls_behavior/PSF_pls_discon_salience_residY_age_sex_cortex100_newmasks_20260915.csv'
EXPR_CSV = BASE + '/PSF_ahba_100/PSF_ahba_100_expression.csv'
FILTER_CSV = BASE + '/PSF_ahba_100/hpa_cortex_detected_genes_consensus.csv'
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

# ---- 1. 对齐 99 区 + HPA 过滤基因子集 ----
sal = pd.read_csv(SAL_CSV)
expr = pd.read_csv(EXPR_CSV, index_col=0)
flt = set(pd.read_csv(FILTER_CSV)['gene'].astype(str))
sal_map = dict(zip(sal['region'], sal['salience']))
expr_valid = ~expr.isna().all(axis=1)
common = [r for r in expr.index if r in sal_map and expr_valid[r]]
common_expr = expr.loc[common].astype(float)
common_sal = np.array([sal_map[r] for r in common])
# 基因宇宙：expr 列 ∩ HPA cortex detected（保持 expr 列序）
genes_all = np.array([str(g) for g in expr.columns])
keep_mask = np.array([g in flt for g in genes_all])
genes = genes_all[keep_mask]
common_expr = common_expr.loc[:, genes]
print(f'[1] 对齐区 {len(common)} | 基因宇宙: {len(genes_all)} → HPA 过滤后 {len(genes)}（-{len(genes_all)-len(genes)}）', flush=True)

# ---- 2. r_gene ----
rs = rankdata(common_sal); rs_z = (rs - rs.mean()) / (rs.std(ddof=1) + 1e-12)
re = np.apply_along_axis(rankdata, 0, common_expr.values)
re_z = (re - re.mean(0)) / (re.std(0, ddof=1) + 1e-12)
r_gene = (rs_z @ re_z) / (len(common) - 1)
print(f'[2] r_gene 完成 n={len(r_gene)}', flush=True)

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
    df.to_csv(f'{OUT}/PSF_discon_GSEA_hpafilter_newmasks_20260915_{gmt_name}.csv', index=False, encoding='utf-8-sig')
    n_sig = int((df['fdr'] < 0.05).sum())
    n_nom = int((df['perm_p'] < 0.05).sum())
    print(f'[{gmt_name}] {len(df)} 通路 | 名义 p<0.05: {n_nom} | FDR<0.05: {n_sig} | {time.time()-t0:.0f}s', flush=True)
    for _, r in df[df['fdr'] < 0.05].head(15).iterrows():
        print(f'    {r["pathway"][:70]:<72s} ES={r["ES"]:+.4f} p={r["perm_p"]:.4f} FDR={r["fdr"]:.4f}', flush=True)
    report.append((gmt_name, len(df), n_nom, n_sig))

with open(f'{OUT}/PSF_discon_GSEA_hpafilter_newmasks_20260915_report.txt', 'w', encoding='utf-8') as f:
    f.write('GSEA（HPA 脑表达过滤敏感性版）\n')
    f.write('方法: Spearman r_gene; ES=基因集内均值; 基因置换10k seed=0 双尾; BH-FDR（与主分析同参）\n')
    f.write(f'HPA 过滤: consensus normal tissue cerebral cortex nTPM>=1 detected（官方阈值 nTPM<1=Not detected）\n')
    f.write(f'宇宙: {len(genes_all)} → {len(genes)} 基因 | 区: {len(common)}\n\n')
    for name, n, n_nom, n_sig in report:
        f.write(f'{name}: {n} 通路 | 名义 p<0.05: {n_nom} | FDR<0.05: {n_sig}\n')
    f.write('\n对照：未过滤版（主分析）GO FDR<0.05=183 / KEGG=40\n')
    f.write('参考文献: Seidlitz 2020 Nat Commun 11:3358; Uhlen 2015 Science PMID 25613900\n')
print('完成。', flush=True)
