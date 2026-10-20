# -*- coding: utf-8 -*-
"""
GSEA 统一输入版（w1 = Step3 PLS 载荷, 皮层 100 区·主分析口径）
=====================================================================
背景（2026-09-04, 统一口径决策候选）：
  原 GSEA 输入为脚本内重算的 r_gene = Spearman(salience, 表达)（crosscheck 版），
  与 Step3/Step4 消费的 PLS w1 载荷排序文件（PSF_pls_gene_ranking_cortex100_residY_newmasks_20260915.csv）
  非同一文件（corr(w1, Spearman)≈0.986, 排序近等价但含 4.23% 并列值→gseapy ties 警告）。
  现改为直接读 Step3 排序 CSV 的 w1 载荷作为 prerank score——
  (1) 实现"Step3 PLS 排序 = 细胞富集 + GSEA 统一输入"（文件级一致）；
  (2) w1 为 SVD 连续浮点载荷, 15,633 值全唯一, 从源头消除 ties 报警。
  数学依据：单响应 PLS(X 每列 Z-score) 载荷 w1 ∝ Xᵀy = 每基因 Pearson(salience, 表达)。

方法: gseapy.prerank（Subramanian 2005 官方加权 KS 游走 ES + 基因置换 1000×, seed=0,
       min_size=5, max_size=300, permutation_type='gene_set'）——与官方主口径同参。
两个宇宙：
  [A] 全宇宙 15,633（主口径）
  [B] HPA cortex detected 过滤 12,753（Seidlitz 2020 同法, 官方阈值 nTPM≥1, 敏感性）

输出: PSF_ahba_100/gseapy_prerank_cortex100_residY_w1_newmasks_20260915_{GO,KEGG}.csv        （全宇宙）
      PSF_ahba_100/gseapy_prerank_cortex100_residY_w1_newmasks_20260915_hpafilter_{GO,KEGG}.csv （HPA 宇宙）
"""
import time
import pandas as pd
import gseapy as gp

BASE = 'F:/dcmtonii/run62_crosssectional'
GMT = BASE + '/PSF_references/gene_sets'
OUT = BASE + '/PSF_ahba_100'
RK_CSV = OUT + '/PSF_pls_gene_ranking_cortex100_residY_newmasks_20260915.csv'
HPA_CSV = OUT + '/hpa_cortex_detected_genes_consensus.csv'
N_PERM = 1000
SEED = 0

# ---- 1. 统一输入：读 Step3 PLS w1 载荷排序（rank 1 = w1 最正端）----
rk = pd.read_csv(RK_CSV)
rnk_full = rk[['gene', 'loading']].rename(columns={'loading': 'score'}).astype({'gene': str})
rnk_full = rnk_full[rnk_full['score'].notna()]
print(f'[A] 全宇宙 w1: n={len(rnk_full)} | score 唯一值 {rnk_full["score"].nunique()} '
      f'| 并列 {int(rnk_full["score"].duplicated().sum())}', flush=True)

flt = set(pd.read_csv(HPA_CSV)['gene'].astype(str))
rnk_hpa = rnk_full[rnk_full['gene'].isin(flt)].reset_index(drop=True)
print(f'[B] HPA 宇宙 w1: n={len(rnk_hpa)} | score 唯一值 {rnk_hpa["score"].nunique()} '
      f'| 并列 {int(rnk_hpa["score"].duplicated().sum())}', flush=True)

# ---- 2. gseapy prerank（与官方同参）----
def run(label, rnk, suffix):
    for name, fn in [('GO', 'GO_Biological_Process_2021.gmt'), ('KEGG', 'KEGG_2021_Human.gmt')]:
        t0 = time.time()
        pre = gp.prerank(rnk=rnk, gene_sets=f'{GMT}/{fn}', outdir=None,
                         min_size=5, max_size=300,
                         permutation_num=N_PERM, seed=SEED, threads=4,
                         permutation_type='gene_set', no_plot=True,
                         verbose=False)
        res = pre.res2d.reset_index(drop=True)
        res.to_csv(f'{OUT}/gseapy_prerank_cortex100_residY_w1_newmasks_20260915{suffix}_{name}.csv',
                   index=False, encoding='utf-8-sig')
        n_sig = int((res['FDR q-val'] < 0.05).sum())
        print(f'[{label}][{name}] 通路 {len(res)} | FDR<0.05: {n_sig} | {time.time()-t0:.0f}s', flush=True)
        for kw in ['oxidative phosphorylation', 'electron transport', 'respiratory',
                   'mitochondr', 'ATP', 'ribosom', 'glycol', 'microglial',
                   'interferon', 'antigen', 'myelin', 'oligodendrocyte', 'immune']:
            hit = res[res['Term'].str.contains(kw, case=False)]
            if len(hit):
                row = hit.sort_values('FDR q-val').iloc[0]
                print(f'    [{kw}] {row["Term"][:52]:<54s} NES={row["NES"]:+.3f} '
                      f'FDR={row["FDR q-val"]:.4f}', flush=True)

run('[A] full', rnk_full, '')
run('[B] hpa ', rnk_hpa, '_hpafilter')
print('完成。', flush=True)
