# -*- coding: utf-8 -*-
"""
LODO 敏感性：公共有效区重算（2026-09-02）
=====================================================================
背景：主 LODO 中 lodo_9861(r=0.747)/lodo_10021(r=0.829) 偏低，
      诊断根因=这两个是 AHBA 仅有的双侧供体（abagen 文档：9861/10021
      have samples from both hemispheres; all others left hemisphere），
      排除它们后右半球 4-5 区无样本覆盖 → 全 NaN → 有效区 94/95。
      若低 r 只是"区覆盖缺失"的伪影，则所有矩阵限制在公共有效区
      （94 区）重跑 PLS 后，相关应大幅回升（>0.9）。

方法：7 个矩阵（6 LODO + 全供体）一律用公共有效区（94 区）
      → 每基因 Z-score → R pls SIMPLS（Y=salience）→ 载荷
      → LODO vs 全供体 Pearson r（与主 LODO 同法）
输出: PSF_ahba_100/lodo/PSF_lodo_sensitivity_common94_report.txt

口径注记（2026-10-10）：公共区数按数据动态取交；终稿结果 = 90 个公共皮层区，
见 derived_data/lodo_common90_r_by_donor.csv（文件名 common94 系初期 94 区口径沿用）。
"""
import os
import numpy as np
import pandas as pd
import subprocess
import tempfile

BASE = 'F:/dcmtonii/run62_crosssectional'
OUT = BASE + '/PSF_ahba_100/lodo_newmasks_20260915'
RSCRIPT = 'I:/R/bin/Rscript.exe'
SAL_CSV = BASE + '/PSF_pls_behavior/PSF_pls_discon_salience_residY_age_sex_cortex100_newmasks_20260915.csv'
FULL_EXPR = BASE + '/PSF_ahba_100/PSF_ahba_100_expression.csv'
DONORS = ['9861', '10021', '12876', '14380', '15496', '15697']


def pls_loadings(expr_df, sal, tmp_dir):
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
    '''
    r_script = r_script.replace('@EXPR@', expr_csv).replace('@YFILE@', y_csv).replace('@OUT@', load_csv)
    r_file = os.path.join(tmp_dir, '_pls.R')
    with open(r_file, 'w', encoding='utf-8') as f:
        f.write(r_script)
    r = subprocess.run([RSCRIPT, r_file], capture_output=True, text=True)
    if r.returncode != 0:
        raise RuntimeError(r.stderr[-400:])
    load = pd.read_csv(load_csv)
    return load.set_index('gene')['loading'], len(common)


def main():
    sal = pd.read_csv(SAL_CSV).set_index('region')['salience']
    full = pd.read_csv(FULL_EXPR, index_col=0)
    mats = {'full': full}
    for d in DONORS:
        mats[f'lodo_{d}'] = pd.read_csv(f'{OUT}/PSF_ahba_LODO_{d}.csv', index_col=0)

    # 公共有效区：全供体有效 ∩ 所有 LODO 有效
    valid_sets = [set(r for r in m.index if not m.loc[r].isna().all()) for m in mats.values()]
    common_regions = sorted(set.intersection(*valid_sets))
    print(f'公共有效区: {len(common_regions)} 个', flush=True)
    print('  ', common_regions, flush=True)

    tmp = tempfile.mkdtemp(prefix='lodo_sens_')
    loadings = {}
    for tag, m in mats.items():
        w, n = pls_loadings(m.loc[common_regions], sal, tmp)
        loadings[tag] = w
        print(f'  [PLS {tag}] {len(w)} 基因 | 区 {n}', flush=True)

    load_df = pd.DataFrame(loadings)
    load_df.to_csv(f'{OUT}/PSF_lodo_sensitivity_common94_loadings.csv', encoding='utf-8-sig')
    cg = load_df.dropna().index
    rs = {d: float(np.corrcoef(load_df.loc[cg, 'full'], load_df.loc[cg, f'lodo_{d}'])[0, 1])
          for d in DONORS}
    rvals = np.array(list(rs.values()))
    lines = [
        'LODO 敏感性（公共 94 区重算）——验证低 r 是否区覆盖伪影',
        f'公共有效区: {len(common_regions)}（{common_regions}）',
        'r（LODO vs 全供体, 公共区约束）:',
    ]
    for d in DONORS:
        lines.append(f'  lodo_{d}: r = {rs[d]:+.4f}')
    lines.append(f'  范围: [{rvals.min():.4f}, {rvals.max():.4f}] | 均值: {rvals.mean():.4f}')
    lines.append(f'  判定: {"✓ 全部≥0.90 → 主LODO低r为右半球覆盖缺失伪影" if rvals.min() >= 0.90 else "⚠️ 仍有低于0.90 → 需进一步讨论"}')
    print('\n' + '\n'.join(lines), flush=True)
    with open(f'{OUT}/PSF_lodo_sensitivity_common94_report.txt', 'w', encoding='utf-8') as f:
        f.write('\n'.join(lines))


if __name__ == '__main__':
    main()
