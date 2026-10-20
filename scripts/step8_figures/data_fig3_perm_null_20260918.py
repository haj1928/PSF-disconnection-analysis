# -*- coding: utf-8 -*-
"""
补算 5：LV1 奇异值的置换空分布（供 Fig. 3a 新版使用）
=====================================================

背景
  原 Fig. 4a 是碎石图（各 LV 的方差解释），信息量低（LV1 一根柱 + 8 根几乎为零的柱）。
  郝医生 2026-09-18 决定改为**置换空分布**：把"观测 LV1 奇异值 vs 5,000 次置换的零分布"
  画出来，直接可视化 Methods 2.5 的奇异值置换检验（正文置换 p = 0.023）。

为什么要单独跑这一遍
  pyls 只返回 p 值（`permres['pvals']`）与置换索引（`permres['permsamples']`），
  **不返回置换后的奇异值矩阵 `d_perm`**——而画空分布正需要它。
  故此处以 monkeypatch 包住 `pyls.compute.perm_sig`，把传给它的 `perm`（L×P）截获下来。
  这是"只取内部量、不改 pyls 行为"的做法：spy 函数原样调用原函数并返回其结果。

管道与数值口径
  完全沿用 `data_fig3_sensitivity_20260916.py`（同一份 X / Y / 协变量 / seed / n_perm），
  并断言重算得到的置换 p 与正文所用值（0.0234）一致——不一致即中止。

输出
  fig3_perm_null_singvals.csv  —— 9 列 × 5000 行（各 LV 的置换空分布奇异值）
  fig3_perm_null_summary.csv   —— 每个 LV 的：观测奇异值 / 空分布均值·SD·最大 / 置换 p
"""
import os
import sys

import numpy as np
import pandas as pd

if not hasattr(np, 'row_stack'):
    np.row_stack = np.vstack
import pyls
import pyls.compute

HERE = os.path.dirname(os.path.abspath(__file__))
BASE = 'F:/dcmtonii/run62_crosssectional'
SRC = 'F:/dcmtonii'
REV = BASE + '/PSF_pls_behavior/断连稿写作_20260907/修订_revision_20260910'
N_PERM = N_BOOT = 5000
SEED = 0
P_EXPECTED = 0.023395320935812838        # 正文 §3.2 主分析置换 p（= 117/5001）

# ---------------------------------------------------------------- 数据（与 data_fig3_sensitivity 完全相同）
mt = pd.read_csv(BASE + '/PSF_master_table_baseline_acute.csv')
mt['sub'] = [f'sub_{x}' if not str(x).startswith('sub_') else str(x) for x in mt['subject_id']]

pd_ = pd.read_csv(SRC + '/PSF_SDC_parcel_discon_N124_cortex100_newmasks_20260915.csv')
ids = [str(x) for x in pd_['ID']]
X = pd_.drop(columns=['ID']).values.astype(float)
Y = pd.read_csv(BASE + '/FSS9_baseline_124.csv', index_col=0).loc[ids].values.astype(float)

mt_i = mt.set_index('sub').loc[ids].copy()
for c in ['年龄', '性别_男1女0']:
    if mt_i[c].isna().any():
        mt_i[c] = mt_i[c].fillna(mt_i[c].median())
age = mt_i['年龄'].values.astype(float)
sex = mt_i['性别_男1女0'].values.astype(float)


def residualize(Xa, C):
    C = np.column_stack([np.asarray(c, float).ravel() for c in C])
    Cw = np.column_stack([np.ones(C.shape[0]), C])
    Xa = np.asarray(Xa, float)
    beta, *_ = np.linalg.lstsq(Cw, Xa, rcond=None)
    return Xa - Cw @ beta


# ---------------------------------------------------------------- 截获 d_perm
_BOX = {}


def _spy(orig_singvals, perm):
    _BOX['orig'] = np.asarray(orig_singvals, float).copy()
    _BOX['perm'] = np.asarray(perm, float).copy()
    return _ORIG_PERM_SIG(orig_singvals, perm)


_ORIG_PERM_SIG = pyls.compute.perm_sig
pyls.compute.perm_sig = _spy

print('运行主分析 PLS（%d 次置换，seed=%d）…' % (N_PERM, SEED), flush=True)
res = pyls.behavioral_pls(X, residualize(Y, [age, sex]), n_perm=N_PERM, n_boot=N_BOOT,
                          seed=SEED, verbose=False, test_split=None)
pyls.compute.perm_sig = _ORIG_PERM_SIG          # 复原

assert 'perm' in _BOX, '未截获 d_perm —— pyls 内部调用路径可能变了，需重查'
d_perm = _BOX['perm']                            # (L, P)
# 注意：pyls 内部把 `singvals` 以 **(L, L) 对角阵**形式传给 perm_sig，
# 而结果对象 `res.singvals` 是**一维 (L,)** —— 两者取法不同，勿混。
obs_diag = np.diag(_BOX['orig'])                 # (L,) ← 与 perm_sig 内部判据一致
obs_sv = np.asarray(res.singvals).ravel()        # (L,)
assert obs_sv.ndim == 1 and obs_sv.size == obs_diag.size, \
    'singvals 形状异常：结果对象 %s vs 截获 %s' % (obs_sv.shape, obs_diag.shape)
assert np.allclose(obs_diag, obs_sv), '截获的观测奇异值与结果对象不一致'
assert d_perm.shape == (len(obs_sv), N_PERM), '空分布形状异常：%s' % (d_perm.shape,)

pvals = (np.sum(d_perm > obs_sv[:, None], axis=1) + 1) / (N_PERM + 1)
assert abs(pvals[0] - P_EXPECTED) < 1e-12, \
    'LV1 置换 p = %.10f，期望 %.10f —— 数据或 seed 口径已变，须重查' % (pvals[0], P_EXPECTED)
assert np.allclose(pvals, np.asarray(res.permres['pvals'])), '与 pyls 返回的 pvals 不一致'

print('LV1 观测奇异值 = %.6f ｜ 空分布：均值 %.4f，SD %.4f，最大 %.4f，超观测者 %d/%d'
      % (obs_sv[0], d_perm[0].mean(), d_perm[0].std(ddof=1), d_perm[0].max(),
         int((d_perm[0] > obs_sv[0]).sum()), N_PERM))
print('LV1 置换 p = %.6f（正文 0.023）｜ 超观测次数 + 1 = %d'
      % (pvals[0], int((d_perm[0] > obs_sv[0]).sum()) + 1))

# ---------------------------------------------------------------- 落盘
n_lv = d_perm.shape[0]
pd.DataFrame({('LV%d' % (i + 1)): d_perm[i] for i in range(n_lv)}).to_csv(
    os.path.join(HERE, 'fig3_perm_null_singvals.csv'), index=False, encoding='utf-8-sig')

summary = pd.DataFrame(dict(
    LV=np.arange(1, n_lv + 1),
    observed_singval=obs_sv,
    null_mean=d_perm.mean(axis=1),
    null_sd=d_perm.std(axis=1, ddof=1),
    null_max=d_perm.max(axis=1),
    n_exceed=np.sum(d_perm > obs_sv[:, None], axis=1),
    perm_p=pvals,
))
summary.to_csv(os.path.join(HERE, 'fig3_perm_null_summary.csv'), index=False, encoding='utf-8-sig')

print('\n=== 各 LV 摘要 ===')
print(summary.to_string(index=False, float_format=lambda v: '%.4f' % v))
print('\n已写出 fig3_perm_null_singvals.csv（%d×%d）与 fig3_perm_null_summary.csv' % (N_PERM, n_lv))
