# -*- coding: utf-8 -*-
"""
LODO 留一供体 PLS 分析（2026-09-02）——AHBA 供体稳健性检验
=====================================================================
背景：用户要求补做 Seidlitz 2020 的 leave-one-donor-out（LODO）分析，
      验证 PLS 载荷（下游所有分子分析的基石）不受单一 AHBA 供体支配。

【文献背书】（全部联网核实）
1. Seidlitz J, et al. Transcriptomic and cellular decoding of regional brain
   vulnerability to neurogenetic disorders. Nat Commun 2020;11:3358.
   DOI 10.1038/s41467-020-17051-5. PMID 32620757.
   原文（Results, LODO 段）：
   "To ensure robustness across the six donors in the AHBA, we performed a
   leave-one-donor-out PLS analysis for each syndrome's brain map, which
   showed highly consistent PLS loadings across all permutations relative
   to the PLS loadings derived from the whole AHBA dataset
   (+X mean r = 0.98; +Y mean r = 0.84; -X mean r = 0.97; +21 mean r = 0.98;
    -22q11 mean r = 0.96; -11p13 mean r = 0.98; all N = 15,043 genes)"
   → 金标准：LODO 载荷 vs 全供体载荷 Pearson r ∈ [0.84, 0.98]。
2. Markello RD, et al. Standardizing workflows in imaging transcriptomics
   with the abagen toolbox. eLife 2021;10:e72129.（abagen 官方）
   原文推荐："analyses are repeated six times, omitting one donor each
   time, to ensure that none of the donors are unduly influencing analytic
   estimates (Vogel et al., 2020; Arnatkeviciute et al., 2019)"
3. 作者官方代码仓库 github.com/jms290/PolySyn_MSNs 仅含绘图脚本与合并后
   数据（无 donor 级数据、无 LODO 代码）→ 无法直接复用作者代码；
   本脚本按原文方法自实现（abagen 官方工具生成 donor 级表达）。

【方法】（逐字对照 Seidlitz 2020 + abagen 官方流程）
  ① 6 个 AHBA 供体（9861/10021/12876/14380/15496/15697）
  ② 每次排除 1 个供体，对剩余 5 供体跑完整 abagen pipeline
     （全默认参数，与 gen_ahba_100.py 一致）→ LODO 表达矩阵（皮层 100 区）
  ③ 全供体表达 = 现有 PSF_ahba_100_expression.csv（6 供体 abagen 产物）
  ④ 每矩阵：每基因跨区 Z-score（Seidlitz 原文）→ R pls SIMPLS ncomp=1
     （Y = 断连-疲劳 salience residY 版）→ 载荷 w1（15633 基因）
  ⑤ corr(LODO w1, 全供体 w1) → 6 个 r，报告范围/均值
  ⑥ 对照金标准 r ∈ [0.84, 0.98]；并验证全供体载荷与主分析
     PSF_pls_gene_ranking_cortex100_residY_newmasks_20260915.csv 一致（自洽检查）

输入: F:/dcmtonii/abagen-data（6 供体）/ extdata/Schaefer_Yeo_Plus_Subcort/
      PSF_ahba_100/PSF_ahba_100_expression.csv（全供体）
      PSF_pls_behavior/PSF_pls_discon_salience_residY_age_sex_cortex100_newmasks_20260915.csv
输出: PSF_ahba_100/lodo/PSF_ahba_LODO_{donor}.csv（6 个 LODO 表达矩阵）
      PSF_ahba_100/lodo/PSF_lodo_pls_loadings.csv（7 组载荷）
      PSF_ahba_100/lodo/PSF_lodo_report.txt（结果报告）
=====================================================================
"""
import os
import numpy as np
import pandas as pd
import subprocess
import tempfile

BASE = 'F:/dcmtonii/run62_crosssectional'
ATLAS = BASE + '/extdata/Schaefer_Yeo_Plus_Subcort/100Parcels7Networks.nii.gz'
INFO = BASE + '/extdata/Schaefer_Yeo_Plus_Subcort/100Parcels7Networks_abagen_info.csv'
DATA = 'F:/dcmtonii/abagen-data'
OUT = BASE + '/PSF_ahba_100/lodo_newmasks_20260915'
RSCRIPT = 'I:/R/bin/Rscript.exe'
SAL_CSV = BASE + '/PSF_pls_behavior/PSF_pls_discon_salience_residY_age_sex_cortex100_newmasks_20260915.csv'
FULL_EXPR = BASE + '/PSF_ahba_100/PSF_ahba_100_expression.csv'
MAIN_RANK = BASE + '/PSF_ahba_100/PSF_pls_gene_ranking_cortex100_residY_newmasks_20260915.csv'

DONORS = ['9861', '10021', '12876', '14380', '15496', '15697']
os.makedirs(OUT, exist_ok=True)

import abagen


def run_abagen(donor_subset, tag):
    """对给定供体子集跑完整 abagen pipeline，返回皮层 100 区表达矩阵（index=label）"""
    info_full = pd.read_csv(INFO)
    cortex_ids = info_full[info_full['structure'] == 'cortex']['id'].tolist()
    print(f'[{tag}] abagen 处理供体 {donor_subset} ...', flush=True)
    expr = abagen.get_expression_data(
        ATLAS, atlas_info=info_full, data_dir=DATA,
        donors=donor_subset, verbose=1, n_proc=1,
    )
    expr_cortex = expr.loc[[i for i in cortex_ids if i in expr.index]].copy()
    id2label = info_full.set_index('id')['label'].to_dict()
    expr_cortex.index = [id2label.get(int(i), i) for i in expr_cortex.index]
    expr_cortex.index.name = 'region'
    return expr_cortex


def pls_loadings(expr_df, sal, tmp_dir):
    """每基因 Z-score → R pls SIMPLS（Y=salience）→ 载荷 w1（与 step4b 完全同法）"""
    # 对齐：表达有效区 ∩ salience 非 NaN 区
    expr_valid = ~expr_df.isna().all(axis=1)
    common = [r for r in expr_df.index if r in sal.index and expr_valid[r]
              and not np.isnan(sal[r])]
    common_expr = expr_df.loc[common].astype(float)
    common_sal = sal.loc[common].values.astype(float)
    gmean = common_expr.mean(axis=0); gstd = common_expr.std(axis=0, ddof=1)
    expr_z = (common_expr - gmean) / (gstd + 1e-12)

    expr_csv = os.path.join(tmp_dir, 'e.csv').replace('\\', '/')
    y_csv = os.path.join(tmp_dir, 'y.csv').replace('\\', '/')
    load_csv = os.path.join(tmp_dir, 'l.csv').replace('\\', '/')
    expr_z.T.to_csv(expr_csv, index=True)
    pd.DataFrame({'y': common_sal}).to_csv(y_csv, index=False)
    r_script = r'''
    suppressMessages(library(pls))
    expr <- as.matrix(read.csv("@EXPR@", row.names=1))
    y <- read.csv("@YFILE@")[,1]
    fit <- plsr(y ~ t(expr), ncomp=1, method="simpls", center=TRUE, scale=FALSE)
    w1 <- as.vector(coef(fit))
    write.csv(data.frame(gene=rownames(expr), loading=w1), file="@OUT@", row.names=FALSE)
    cat(length(w1))
    '''
    r_script = r_script.replace('@EXPR@', expr_csv).replace('@YFILE@', y_csv).replace('@OUT@', load_csv)
    r_file = os.path.join(tmp_dir, '_pls.R')
    with open(r_file, 'w', encoding='utf-8') as f:
        f.write(r_script)
    r = subprocess.run([RSCRIPT, r_file], capture_output=True, text=True)
    if r.returncode != 0:
        raise RuntimeError(f'R pls 失败: {r.stderr[-500:]}')
    load = pd.read_csv(load_csv)
    return load.set_index('gene')['loading'], len(common)


def main():
    sal = pd.read_csv(SAL_CSV).set_index('region')['salience']

    # ---- 1. 全供体（现有产物，避免重跑）+ 自洽性验证 ----
    full = pd.read_csv(FULL_EXPR, index_col=0)
    print(f'[1] 全供体表达: {full.shape}', flush=True)

    # ---- 2. 6 个 LODO 子集（严格：每次完整 abagen pipeline）----
    lodo_matrices = {}
    for i, d in enumerate(DONORS):
        subset = [x for x in DONORS if x != d]
        tag = f'LODO-{d}'
        m = run_abagen(subset, tag)
        m.to_csv(f'{OUT}/PSF_ahba_LODO_{d}.csv', encoding='utf-8-sig')
        lodo_matrices[f'lodo_{d}'] = m
        print(f'  [{tag}] {m.shape}', flush=True)

    # ---- 3. 7 组 PLS 载荷（6 LODO + 全供体）----
    tmp = tempfile.mkdtemp(prefix='lodo_pls_')
    loadings = {}
    nreg = {}
    for d in DONORS:
        w, n = pls_loadings(lodo_matrices[f'lodo_{d}'], sal, tmp)
        loadings[f'lodo_{d}'] = w
        nreg[f'lodo_{d}'] = n
        print(f'  [PLS lodo_{d}] 载荷 {len(w)} 基因 | 有效区 {n}', flush=True)
    w_full, n_full = pls_loadings(full, sal, tmp)
    loadings['full'] = w_full
    nreg['full'] = n_full
    print(f'  [PLS full] 载荷 {len(w_full)} 基因 | 有效区 {n_full}', flush=True)

    # ---- 4. 保存载荷 ----
    load_df = pd.DataFrame(loadings)
    load_df.to_csv(f'{OUT}/PSF_lodo_pls_loadings.csv', encoding='utf-8-sig')

    # ---- 5. 相关：LODO vs 全供体（Seidlitz 金标准）----
    common_genes = load_df.dropna().index
    print(f'[5] 公共基因 {len(common_genes)}（全供体与 LODO 载荷均有的基因）', flush=True)
    rs = {}
    for d in DONORS:
        key = f'lodo_{d}'
        r = float(np.corrcoef(load_df.loc[common_genes, 'full'],
                              load_df.loc[common_genes, key])[0, 1])
        rs[key] = r
        print(f'  corr(full, {key}) = {r:+.4f}', flush=True)

    # ---- 6. 自洽性：全供体载荷 vs 主分析 ranking（应完全一致）----
    main_rk = pd.read_csv(MAIN_RANK)
    w_main = main_rk.set_index('gene')['loading']
    common2 = w_full.index.intersection(w_main.index)
    r_self = float(np.corrcoef(w_full.loc[common2], w_main.loc[common2])[0, 1])
    maxdiff = float(np.abs(w_full.loc[common2] - w_main.loc[common2]).max())
    print(f'[6] 自洽: corr(全供体载荷, 主分析ranking载荷) = {r_self:+.6f} | max|diff| = {maxdiff:.2e}', flush=True)

    # ---- 7. 报告 ----
    rvals = np.array(list(rs.values()))
    lines = [
        '=' * 78,
        'LODO 留一供体 PLS 分析报告（2026-09-02）',
        '=' * 78,
        '方法: 6 次完整 abagen pipeline（每次排除 1 供体, 全默认参数与 gen_ahba_100 一致）',
        '      → 皮层 100 区表达 → 每基因 Z-score → R pls SIMPLS ncomp=1',
        '      Y = 断连-疲劳 salience（仅残差 age+sex, 主分析版）',
        '文献背书: Seidlitz 2020 Nat Commun 11:3358（LODO 原文）',
        '          Markello 2021 eLife 10:e72129（abagen 官方推荐 LODO）',
        '-' * 78,
        '有效区数（表达非全NaN ∩ salience 有效）:',
    ]
    for d in DONORS:
        lines.append(f'  lodo_{d}: {nreg[f"lodo_{d}"]} 区')
    lines.append(f'  full  : {n_full} 区')
    lines.append('-' * 78)
    lines.append('LODO 载荷 vs 全供体载荷 Pearson r（Seidlitz 金标准 0.84-0.98）:')
    for d in DONORS:
        lines.append(f'  lodo_{d}: r = {rs[f"lodo_{d}"]:+.4f}')
    lines.append(f'  范围: [{rvals.min():.4f}, {rvals.max():.4f}] | 均值: {rvals.mean():.4f}')
    ok = (rvals.min() >= 0.84)
    lines.append(f'  判定: {"✓ 达到 Seidlitz 0.84-0.98 金标准" if ok else "⚠️ 未达金标准（需讨论）"}')
    lines.append('-' * 78)
    lines.append(f'自洽性: corr(全供体载荷, 主分析ranking载荷) = {r_self:+.6f} | max|diff| = {maxdiff:.2e}')
    lines.append(f'  （≈1 说明全供体表达与主分析一致, LODO 流程无偏差）')
    lines.append('-' * 78)
    lines.append('文件:')
    lines.append('  PSF_ahba_LODO_{donor}.csv × 6（LODO 表达矩阵）')
    lines.append('  PSF_lodo_pls_loadings.csv（7 组载荷）')
    with open(f'{OUT}/PSF_lodo_report.txt', 'w', encoding='utf-8') as f:
        f.write('\n'.join(lines))
    print('\n' + '\n'.join(lines))
    print(f'\n完成。报告: {OUT}/PSF_lodo_report.txt')


if __name__ == '__main__':
    main()
