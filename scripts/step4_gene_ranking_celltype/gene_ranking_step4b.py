# -*- coding: utf-8 -*-
"""
Step 4'（皮层 100 区·主分析版）：用新 salience（仅残差Y）重算基因排序 + 生成 Metascape 输入
=====================================================================
背景：主分析定稿为"仅残差化 Y(FSS9) 对[年龄,性别]"（Petrican 2022 先例），
      产生新 salience（PSF_pls_discon_salience_residY_age_sex_cortex100_newmasks_20260915.csv）。
      为保证全链路口径一致，Metascape 基因列表必须基于新 salience 重算，
      弃用旧的（基线未控 salience 版）pos300/neg300。

方法（与 Step 4 完全一致，Seidlitz 2020 原法）：
  X = 皮层表达（99 区有值 × 15633 基因，每基因跨区 Z-score）
  Y = 新 salience（99 区）
  PLS 回归 = R pls 包 SIMPLS, ncomp=1（Seidlitz 原文指定）
  排序 = 载荷 w1 大正→大负
Metascape 输入 = 正共变 top300 / 负共变 bottom300（每行一个基因符号）

输出：
  PSF_ahba_100/PSF_pls_gene_ranking_cortex100_residY_newmasks_20260915.csv
  PSF_ahba_100/PSF_pls_gene_ranking_cortex100_residY_newmasks_20260915_report.txt
  PSF_ahba_100/metascape_input_pos300_residY_newmasks_20260915.txt
  PSF_ahba_100/metascape_input_neg300_residY_newmasks_20260915.txt
=====================================================================
"""
import numpy as np
import pandas as pd
import subprocess
import os
import tempfile

BASE = 'F:/dcmtonii/run62_crosssectional'
SAL = BASE + '/PSF_pls_behavior/PSF_pls_discon_salience_residY_age_sex_cortex100_newmasks_20260915.csv'
EXPR = BASE + '/PSF_ahba_100/PSF_ahba_100_expression.csv'
OUT = BASE + '/PSF_ahba_100'
RSCRIPT = 'I:/R/bin/Rscript.exe'
TOP_N = 300

os.makedirs(OUT, exist_ok=True)

# 1. 加载（新 salience 列名 = region/salience）
sal = pd.read_csv(SAL)
expr = pd.read_csv(EXPR, index_col=0)
print(f'[1] salience(residY): {sal.shape[0]} 区 | expr: {expr.shape[0]} 区 × {expr.shape[1]} 基因', flush=True)

# 2. 对齐：salience ∩ 表达（排除表达全 NaN 区 RH_SomMot_2）
sal_map = dict(zip(sal['region'], sal['salience']))
expr_valid = ~expr.isna().all(axis=1)
common = [r for r in expr.index if r in sal_map and not np.isnan(sal_map[r]) and expr_valid[r]]
common_expr = expr.loc[common].astype(float)
common_sal = np.array([sal_map[r] for r in common])
print(f'[2] 交集区: {len(common)}（排除表达全 NaN 区）', flush=True)

# 3. 每基因跨区 Z-score（Seidlitz 原文 "each gene's regional values were normalized (Z-scored)"）
gmean = common_expr.mean(0)
gstd = common_expr.std(0, ddof=1)
expr_z = (common_expr - gmean) / (gstd + 1e-12)
print(f'[3] Z-score 完成: {expr_z.shape}', flush=True)

# 4. PLS 回归（R pls 包 SIMPLS, ncomp=1）
TMP = tempfile.mkdtemp(prefix='pls_rank100rY_')
expr_csv = os.path.join(TMP, 'expr_z.csv').replace('\\', '/')
y_csv = os.path.join(TMP, 'y.csv').replace('\\', '/')
load_csv = os.path.join(TMP, 'loadings.csv').replace('\\', '/')
expr_z.T.to_csv(expr_csv, index=True)          # 15633 × 99
pd.DataFrame({'y': common_sal}).to_csv(y_csv, index=False)

r_script = r'''
suppressMessages(library(pls))
expr <- as.matrix(read.csv("@EXPR@", row.names=1))   # 15633 × 99
y <- read.csv("@YFILE@")[,1]                          # 99
fit <- plsr(y ~ t(expr), ncomp=1, method="simpls", center=TRUE, scale=FALSE)
w1 <- as.vector(coef(fit))
write.csv(data.frame(gene=rownames(expr), loading=w1), file="@OUT@", row.names=FALSE)
cat("PLS done:", length(w1), "loadings\n")
'''
r_script = r_script.replace('@EXPR@', expr_csv).replace('@YFILE@', y_csv).replace('@OUT@', load_csv)
r_file = os.path.join(TMP, '_pls.R')
with open(r_file, 'w', encoding='utf-8') as f:
    f.write(r_script)
r = subprocess.run([RSCRIPT, r_file], capture_output=True, text=True)
if r.returncode != 0:
    print('R 失败:', r.stderr[-800:]); raise SystemExit(1)
print('[4] R pls SIMPLS 完成:', r.stdout.strip()[-40:], flush=True)

# 5. 读回 + 排序（rank 1 = 最大正载荷）
load_df = pd.read_csv(load_csv)
w1 = load_df['loading'].values.astype(float)
genes = load_df['gene'].astype(str).tolist()
order = np.argsort(-w1)
rank_of = np.empty(len(w1), dtype=int)
rank_of[order] = np.arange(1, len(w1) + 1)

out_df = pd.DataFrame({'gene': genes, 'loading': w1, 'rank': rank_of})
out_df.to_csv(OUT + '/PSF_pls_gene_ranking_cortex100_residY_newmasks_20260915.csv', index=False, encoding='utf-8-sig')

# 6. Metascape 输入（正 top300 / 负 bottom300）
pos = out_df.sort_values('rank').head(TOP_N)['gene'].tolist()
neg = out_df.sort_values('rank').tail(TOP_N)['gene'].tolist()
with open(OUT + '/metascape_input_pos300_residY_newmasks_20260915.txt', 'w', encoding='utf-8') as f:
    f.write('\n'.join(pos))
with open(OUT + '/metascape_input_neg300_residY_newmasks_20260915.txt', 'w', encoding='utf-8') as f:
    f.write('\n'.join(neg))

# 7. 报告
report = [
    '=' * 70,
    "Step 4'（皮层 100 区·主分析版）：基因排序（新 salience = 仅残差Y对年龄+性别）",
    f'样本区: {len(common)} 个皮层区（salience ∩ 表达，排除 RH_SomMot_2）',
    f'基因数: {len(w1)} | Metascape 输入: 正 top{TOP_N} / 负 bottom{TOP_N}',
    '方法: PLS 回归（R pls 包 SIMPLS, ncomp=1, Seidlitz 2020 原法）',
    '-' * 70,
    '正共变 top 15（与疲劳相关断连模式共变的基因）:',
]
for _, r in out_df.sort_values('rank').head(15).iterrows():
    report.append(f'  {r["rank"]:>4}. {r["gene"]:<15s} loading={r["loading"]:+.5f}')
report.append('-' * 70)
report.append('负共变 bottom 15:')
for _, r in out_df.sort_values('rank').tail(15).iterrows():
    report.append(f'  {r["rank"]:>4}. {r["gene"]:<15s} loading={r["loading"]:+.5f}')
report.append('-' * 70)
report.append('Metascape 输入文件:')
report.append(f'  {OUT}/metascape_input_pos300_residY_newmasks_20260915.txt（正共变, 建议先跑）')
report.append(f'  {OUT}/metascape_input_neg300_residY_newmasks_20260915.txt（负共变, 对照）')
with open(OUT + '/PSF_pls_gene_ranking_cortex100_residY_newmasks_20260915_report.txt', 'w', encoding='utf-8') as f:
    f.write('\n'.join(report))
print()
print('\n'.join(report[:20]))
print(f'\n[5] 已保存:\n  {OUT}/PSF_pls_gene_ranking_cortex100_residY_newmasks_20260915.csv')
print(f'  {OUT}/metascape_input_pos300_residY_newmasks_20260915.txt')
print(f'  {OUT}/metascape_input_neg300_residY_newmasks_20260915.txt')
