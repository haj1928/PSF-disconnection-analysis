# -*- coding: utf-8 -*-
"""
补算 1 & 2：PLS 奇异值谱 + 逐例 Lx/Ly（A1 新掩膜口径）

用途：Fig. 4「共变模式的稳健性」的 a（碎石图）与 b（脑-行为散点）面板。
同时补算「无协变量基线」并把结果落盘 —— 该值（p = 0.0170）此前只存在于
`A1执行记录_20260915.md`，**没有原始输出文件**，属于溯源缺口。

口径（与主分析完全一致）
  X = 皮层 100 区断连矩阵（124×100，**不残差化**）
  Y = FSS9（124×9）对 [年龄, 性别] 残差化
  pyls.behavioral_pls(n_perm=5000, n_boot=5000, seed=0)

输出（均于 figures_20260915/）
  fig3_pls_singvals.csv       —— LV 序号 / 奇异值 / 方差解释% / 累计%
  fig3_subject_scores.csv     —— 逐例 Lx、Ly、FSS-9 总分、病灶体积
  fig3_pls_baseline_check.csv —— 主分析 vs 无协变量基线 vs 已知值的核对
"""
import os
import sys

import numpy as np
import pandas as pd

if not hasattr(np, 'row_stack'):
    np.row_stack = np.vstack
import pyls

SRC = 'F:/dcmtonii'
BASE = SRC + '/run62_crosssectional'
HERE = os.path.dirname(os.path.abspath(__file__))

# ---------------------------------------------------------------- 读数据
mt = pd.read_csv(BASE + '/PSF_master_table_baseline_acute.csv')
mt['sub'] = [f'sub_{x}' if not str(x).startswith('sub_') else str(x) for x in mt['subject_id']]

pd_ = pd.read_csv(SRC + '/PSF_SDC_parcel_discon_N124_cortex100_newmasks_20260915.csv')
ids = [str(x) for x in pd_['ID']]
X = pd_.drop(columns=['ID']).values.astype(float)
pnames = [c for c in pd_.columns if c != 'ID']

Y9 = pd.read_csv(BASE + '/FSS9_baseline_124.csv', index_col=0)
Y = Y9.loc[ids].values.astype(float)
n = len(ids)
print('N=%d | X %s | Y %s' % (n, X.shape, Y.shape), flush=True)

mt_i = mt.set_index('sub').loc[ids]
age = mt_i['年龄'].values.astype(float)
sex = mt_i['性别_男1女0'].values.astype(float)


def residualize(Xa, C):
    """X 每列对 C 线性回归取残差（含截距）。C 可为 (n,) / (n,k) / [v1, v2, ...]"""
    if isinstance(C, (list, tuple)):
        C = np.column_stack([np.asarray(c, float).ravel() for c in C])
    C = np.asarray(C, float)
    if C.ndim == 1:
        C = C[:, None]
    Xa = np.asarray(Xa, float)
    assert C.shape[0] == Xa.shape[0], 'C%r 与 X%s 行数不符' % (C.shape, Xa.shape)
    Cw = np.column_stack([np.ones(C.shape[0]), C])
    beta, *_ = np.linalg.lstsq(Cw, Xa, rcond=None)
    return Xa - Cw @ beta


def run(Xa, Ya, tag):
    res = pyls.behavioral_pls(Xa, Ya, n_perm=5000, n_boot=5000, seed=0,
                              verbose=False, test_split=None)
    sv = np.asarray(getattr(res, 'singvals', None), dtype=float) if hasattr(res, 'singvals') else None
    if sv is None:
        for a in ('s', 'S', 'singular_values'):
            if hasattr(res, a):
                sv = np.asarray(getattr(res, a), dtype=float)
                break
    lx = np.asarray(res.x_scores[:, 0], dtype=float)
    ly = np.asarray(res.y_scores[:, 0], dtype=float)
    p = float(res.permres['pvals'][0])
    r = float(np.corrcoef(lx, ly)[0, 1])
    ratio = res.bootres['x_weights_normed'][:, 0]
    nsig = int(np.sum(np.abs(ratio) > 2))
    print('  [%s] perm_p=%.4f corr=%+.4f n_sig=%d | singvals=%s'
          % (tag, p, r, nsig, None if sv is None else np.round(sv[:5], 3).tolist()), flush=True)
    return dict(tag=tag, perm_p=p, corr=r, n_sig=nsig), sv, lx, ly, res


print('\n=== 主分析（Y 残差化 [年龄,性别]，X 不动）===', flush=True)
main, sv, lx, ly, res = run(X, residualize(Y, [age, sex]), 'main')

print('\n=== 无协变量基线（X、Y 均不残差化）===', flush=True)
base, sv_b, lx_b, ly_b, _ = run(X, Y, 'baseline_no_covariate')

# ---------------------------------------------------------------- 奇异值谱
if sv is None:
    raise SystemExit('pyls 未返回奇异值 —— 请检查属性名')
varexp = sv ** 2 / np.sum(sv ** 2) * 100.0
cum = np.cumsum(varexp)
df_sv = pd.DataFrame(dict(LV=np.arange(1, len(sv) + 1), singval=sv,
                          varexp_pct=varexp, cum_varexp_pct=cum))
out1 = os.path.join(HERE, 'fig3_pls_singvals.csv')
df_sv.to_csv(out1, index=False, encoding='utf-8-sig')
print('\nLV1 方差解释 = %.4f%%（正文口径 94.0%%）' % varexp[0])
print('已写出 %s（%d 个 LV）' % (out1, len(df_sv)))

# ---------------------------------------------------------------- 逐例得分
vol = pd.read_csv(BASE + '/PSF_pls_behavior/断连稿写作_20260907/修订_revision_20260910/'
                         'PSF_lesion_volume_N124_newmasks_20260915.csv', encoding='utf-8-sig')
vol_i = vol.set_index('subject').loc[ids]
df_s = pd.DataFrame(dict(subject=ids, Lx=lx, Ly=ly,
                         FSS9_total=Y.sum(1),
                         lesion_volume_mL=vol_i['volume_mL'].values.astype(float)))
out2 = os.path.join(HERE, 'fig3_subject_scores.csv')
df_s.to_csv(out2, index=False, encoding='utf-8-sig')
print('已写出 %s' % out2)

# ---------------------------------------------------------------- 核对
checks = [
    ('main', 'perm_p', 0.0234), ('main', 'corr', 0.3174), ('main', 'n_sig', 24),
    ('baseline', 'perm_p', 0.0170),
]
rows = []
rows.append(dict(item='主分析 perm_p', got=main['perm_p'], expected=0.0234,
                 ok=abs(main['perm_p'] - 0.0234) < 5e-4))
rows.append(dict(item='主分析 corr', got=main['corr'], expected=0.3174,
                 ok=abs(main['corr'] - 0.3174) < 5e-4))
rows.append(dict(item='主分析 n_sig', got=main['n_sig'], expected=24, ok=main['n_sig'] == 24))
rows.append(dict(item='主分析 varexp LV1', got=round(float(varexp[0]), 4), expected=94.0,
                 ok=abs(varexp[0] - 94.0) < 0.05))
rows.append(dict(item='无协变量基线 perm_p', got=base['perm_p'], expected=0.0170,
                 ok=abs(base['perm_p'] - 0.0170) < 5e-4))
out3 = os.path.join(HERE, 'fig3_pls_baseline_check.csv')
pd.DataFrame(rows).to_csv(out3, index=False, encoding='utf-8-sig')
print('\n=== 核对（复现主分析口径）===')
for r in rows:
    print('  %-22s 本脚本 %-10s 已知 %-10s %s'
          % (r['item'], r['got'] if not isinstance(r['got'], float) else round(r['got'], 4),
             r['expected'], 'OK' if r['ok'] else '**FAIL**'))
print('已写出 %s' % out3)
