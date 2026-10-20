# -*- coding: utf-8 -*-
"""
Step7c 对照版（2026-09-02）：断连"对照向量"的 GSEA——疲劳特异性对照
=====================================================================
背景：用断连本身（PC1/mean）替换疲劳 salience 重跑 GSEA，检验通路层
      （OXPHOS/糖酵解等）是否疲劳维度特有。
方法（与主分析版 run_gsea_discon_cortex100_residY.py 逐字一致）：
  - 每基因 r_gene = Spearman(对照向量, 基因表达) 跨 99 皮层区 (15633,)
  - ES = 基因集内 r_gene 均值（Tian 2005 competitive test）
  - null = 基因置换 10,000 次, seed=0, 双尾 |ES|>=|obs|
  - BH-FDR（GO 与 KEGG 分别校正）；基因集大小过滤 5-300
输入: PSF_ahba_100/PSF_discon_control_vectors_newmasks_20260915.csv × PSF_ahba_100/PSF_ahba_100_expression.csv
输出: PSF_ahba_100/PSF_discon_GSEA_control_newmasks_20260915_{pc1,mean}_{GO,KEGG}.csv + report
"""
import numpy as np, pandas as pd, os, time
from scipy.stats import rankdata

BASE = 'F:/dcmtonii/run62_crosssectional'
CTRL = BASE + '/PSF_ahba_100/PSF_discon_control_vectors_newmasks_20260915.csv'
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

ctrl = pd.read_csv(CTRL)
expr = pd.read_csv(EXPR_CSV, index_col=0)
expr_valid = ~expr.isna().all(axis=1)
common = [r for r in expr.index if r in set(ctrl['region']) and expr_valid[r]]
common_expr = expr.loc[common].astype(float)
genes = np.array([str(g) for g in expr.columns])
print(f'[1] 对齐区: {len(common)}', flush=True)

for vec_name in ['discon_pc1', 'discon_mean']:
    vec = ctrl.set_index('region').loc[common, vec_name].values.astype(float)
    rs = rankdata(vec); rs_z = (rs - rs.mean()) / (rs.std(ddof=1) + 1e-12)
    re = np.apply_along_axis(rankdata, 0, common_expr.values)
    re_z = (re - re.mean(0)) / (re.std(0, ddof=1) + 1e-12)
    r_gene = (rs_z @ re_z) / (len(common) - 1)
    print(f'[2] {vec_name} r_gene 完成, 范围 [{r_gene.min():.3f}, {r_gene.max():.3f}]', flush=True)

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
        df.to_csv(f'{OUT}/PSF_discon_GSEA_control_newmasks_20260915_{vec_name}_{gmt_name}.csv',
                  index=False, encoding='utf-8-sig')
        n_sig = int((df['fdr'] < 0.05).sum())
        print(f'[{vec_name}|{gmt_name}] {len(df)} 通路 | FDR<0.05: {n_sig} | {time.time()-t0:.0f}s', flush=True)
        for _, r in df[df['fdr'] < 0.05].head(15).iterrows():
            print(f'    {r["pathway"][:66]:<68s} ES={r["ES"]:+.4f} FDR={r["fdr"]:.4f}', flush=True)
        report.append((gmt_name, len(df), n_sig))

    with open(f'{OUT}/PSF_discon_GSEA_control_newmasks_20260915_{vec_name}_report.txt', 'w', encoding='utf-8') as f:
        f.write(f'Step7c 对照版：断连对照向量 {vec_name} GSEA（皮层100区）\n')
        f.write('方法: Spearman(对照向量, 表达) r_gene; ES=基因集内均值; 基因置换10k seed=0 双尾; BH-FDR\n')
        for name, n, n_sig in report:
            f.write(f'{name}: {n} 通路 | FDR<0.05: {n_sig}\n')
print('全部完成。', flush=True)
