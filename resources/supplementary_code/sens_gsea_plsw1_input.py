# -*- coding: utf-8 -*-
"""敏感性：GSEA 改用 Step3 的 PLS w1 作 prerank score（严格"统一输入"），对比主口径 Spearman r_gene 版结论
同参：gseapy.prerank, 1000 perm, seed=0, gene_set 置换, min 5/max 300, ascending=False
对比：FDR<0.05 显著通路集一致性 + NES 相关
"""
import numpy as np, pandas as pd
import gseapy as gp

BASE = 'F:/dcmtonii/run62_crosssectional'
GMT = BASE + '/PSF_references/gene_sets'
OUT = BASE + '/PSF_ahba_100'
N_PERM = 1000
SEED = 0

# ---- 1. PLS w1 排序（Step3 产物，细胞富集同源输入）----
rk = pd.read_csv(BASE + '/PSF_ahba_100/PSF_pls_gene_ranking_cortex100_residY_newmasks_20260915.csv')
rnk_pls = pd.DataFrame({'gene': rk['gene'].astype(str), 'score': rk['loading'].astype(float)})
print(f'[1] PLS w1 排序就绪: {len(rnk_pls)} 基因', flush=True)

summary = {}
for name in ['GO', 'KEGG']:
    fn = f'{GMT}/{"GO_Biological_Process_2021" if name=="GO" else "KEGG_2021_Human"}.gmt'
    pre = gp.prerank(rnk=rnk_pls, gene_sets=fn, outdir=None,
                     min_size=5, max_size=300,
                     permutation_num=N_PERM, seed=SEED, threads=4,
                     permutation_type='gene_set', no_plot=True,
                     ascending=False)
    res = pre.res2d.copy()
    res.to_csv(f'{OUT}/gseapy_prerank_cortex100_residY_newmasks_20260915_{{name}}_plsw1_input.csv', index=False)

    # 主口径对照（Spearman r_gene 版，已有产物）
    ref = pd.read_csv(f'{OUT}/gseapy_prerank_cortex100_residY_{name}.csv')
    # 按 Term 对齐
    m = pd.merge(ref[['Term', 'NES', 'FDR q-val']].rename(columns={'NES': 'NES_ref', 'FDR q-val': 'FDR_ref'}),
                 res[['Term', 'NES', 'FDR q-val']].rename(columns={'NES': 'NES_pls', 'FDR q-val': 'FDR_pls'}),
                 on='Term', how='inner')
    sig_ref = set(m.loc[m['FDR_ref'] < 0.05, 'Term'])
    sig_pls = set(m.loc[m['FDR_pls'] < 0.05, 'Term'])
    jac = len(sig_ref & sig_pls) / len(sig_ref | sig_pls) if (sig_ref | sig_pls) else 1.0
    n_missing_ref = len(sig_ref) - len(sig_ref & sig_pls)
    n_new_pls = len(sig_pls) - len(sig_ref & sig_pls)
    corr_nes = np.corrcoef(m['NES_ref'], m['NES_pls'])[0, 1]
    # 代谢核心抽查
    core = m[m['Term'].str.contains('oxidative phosphorylation|ribosom|glycolysis|electron transport|respiratory electron|ATP synth', case=False)]
    core_rows = core[['Term', 'NES_ref', 'FDR_ref', 'NES_pls', 'FDR_pls']].sort_values('NES_pls', ascending=False).head(6)
    print(f'[{name}] ref显著 {len(sig_ref)} | PLS输入显著 {len(sig_pls)} | Jaccard={jac:.3f} | ref丢 {n_missing_ref} | 新增 {n_new_pls} | NES corr={corr_nes:.4f}', flush=True)
    print(f'  代谢核心对照:'); print(core_rows.to_string(index=False), flush=True)
    summary[name] = dict(n_ref=len(sig_ref), n_pls=len(sig_pls), jac=jac, n_missing=n_missing_ref, n_new=n_new_pls, nes_corr=corr_nes)

print('\n=== 汇总 ===')
for k, v in summary.items():
    print(f'{k}: {v}')
print('完成。')
