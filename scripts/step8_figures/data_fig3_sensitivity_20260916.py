# -*- coding: utf-8 -*-
"""
补算 4：多口径敏感性总表（A1 新掩膜口径，单一管线重算）

背景
  该项目的敏感性结果散落在多个 txt 报告里，且**口径混杂**：
  `_sens_lesion_volume_report_newmasks_20260915.txt` 名为 newmasks，
  实际仍读旧体积文件（范围 [0.10, 81.68] mL、中位 2.51）与旧断连矩阵；
  其数值与定稿口径不一致。
  → 本脚本用**一条管线**重算全部口径，落成单一 CSV，供 Fig. 3 出图。

统一设置
  X = 皮层 100 区断连矩阵（124×100，`..._newmasks_20260915.csv`）
  Y = FSS9（124×9）；体积 = `PSF_lesion_volume_N124_newmasks_20260915.csv`
  ⚠️ 涉及 log10 体积的三个口径**不设 expected 校验**：原 A1 脚本
     `downstream_newmasks_20260915/_sens_lesion_volume_log_cortex100_newmasks_20260915.py`
     第 38 行仍读**旧**体积文件 `PSF_lesion_volume_N124.csv`（63 例数值不同，
     其中最大病灶一例为 81.682 而非 82.768 mL），故其记录的 0.0326 / 0.0094
     属旧体积口径；本脚本用 A1 正确体积重算，值见输出。
  pyls.behavioral_pls(n_perm=5000, n_boot=5000, seed=0)
  协变量缺失（NIHSS 1 例、HAMD 3 例）→ 中位数插补（与主分析一致）

输出
  fig3_sensitivity_specs.csv   —— 各口径 perm_p / corr / n_sig / 变体类型
  fig3_loading_main_vs_rank.csv —— 主载荷 vs 秩变换载荷（100 区，供散点）

变更（2026-10-10）
  原「口径 3：杠杆点 LOO（单例）」与「口径 4：LOO 阶梯」（按体积剔 A / 按 FSS 剔 B）
  移出本脚本（该检验已于 2026-09-22 从稿件删除）。
"""
import os
import sys

import numpy as np
import pandas as pd
from scipy.stats import rankdata

if not hasattr(np, 'row_stack'):
    np.row_stack = np.vstack
import pyls

SRC = 'F:/dcmtonii'
BASE = SRC + '/run62_crosssectional'
REV = BASE + '/PSF_pls_behavior/断连稿写作_20260907/修订_revision_20260910'
HERE = os.path.dirname(os.path.abspath(__file__))
N_PERM = N_BOOT = 5000
SEED = 0

# ---------------------------------------------------------------- 数据
mt = pd.read_csv(BASE + '/PSF_master_table_baseline_acute.csv')
mt['sub'] = [f'sub_{x}' if not str(x).startswith('sub_') else str(x) for x in mt['subject_id']]

pd_ = pd.read_csv(SRC + '/PSF_SDC_parcel_discon_N124_cortex100_newmasks_20260915.csv')
ids = [str(x) for x in pd_['ID']]
X = pd_.drop(columns=['ID']).values.astype(float)
parcels = [c for c in pd_.columns if c != 'ID']

Y = pd.read_csv(BASE + '/FSS9_baseline_124.csv', index_col=0).loc[ids].values.astype(float)

vol_df = pd.read_csv(REV + '/PSF_lesion_volume_N124_newmasks_20260915.csv', encoding='utf-8-sig')
vol = vol_df.set_index('subject').loc[ids]['volume_mL'].values.astype(float)

mt_i = mt.set_index('sub').loc[ids].copy()
for c in ['年龄', '性别_男1女0', 'HAMD_total']:
    if mt_i[c].isna().any():
        mt_i[c] = mt_i[c].fillna(mt_i[c].median())
age, sex = mt_i['年龄'].values.astype(float), mt_i['性别_男1女0'].values.astype(float)
hamd = mt_i['HAMD_total'].values.astype(float)
logvol = np.log10(vol)
print('N=%d | 体积 %.3f–%.3f mL（最大者 %s）'
      % (len(ids), vol.min(), vol.max(), ids[int(np.argmax(vol))]), flush=True)


def residualize(Xa, C):
    if isinstance(C, (list, tuple)):
        C = np.column_stack([np.asarray(c, float).ravel() for c in C])
    C = np.asarray(C, float)
    if C.ndim == 1:
        C = C[:, None]
    Xa = np.asarray(Xa, float)
    Cw = np.column_stack([np.ones(C.shape[0]), C])
    beta, *_ = np.linalg.lstsq(Cw, Xa, rcond=None)
    return Xa - Cw @ beta


def run(Xa, Ya):
    res = pyls.behavioral_pls(Xa, Ya, n_perm=N_PERM, n_boot=N_BOOT, seed=SEED,
                              verbose=False, test_split=None)
    lx, ly = res.x_scores[:, 0], res.y_scores[:, 0]
    return (float(res.permres['pvals'][0]), float(np.corrcoef(lx, ly)[0, 1]),
            int(np.sum(np.abs(res.bootres['x_weights_normed'][:, 0]) > 2)),
            np.asarray(res.x_weights[:, 0], float))


SPECS = []
def add(name, group, Xs, Ys, expect=None):
    p, r, ns, w = run(Xs, Ys)
    tag = 'OK' if (expect is None or abs(p - expect) < 5e-4) else '**期望%s**' % expect
    print('  %-38s p=%.4f corr=%+.4f n_sig=%3d  %s' % (name, p, r, ns, tag), flush=True)
    SPECS.append(dict(spec=name, group=group, perm_p=p, corr=r, n_sig=ns, expected=expect))
    return w


print('\n=== 口径 1：协变量谱 ===', flush=True)
w_main = add('主分析（Y~年龄+性别）', 'covariate', X, residualize(Y, [age, sex]), 0.0234)
add('无协变量', 'covariate', X, Y, 0.0170)
add('+log10 体积（仅 Y 残差）', 'covariate', X, residualize(Y, [age, sex, logvol]))
# 2026-09-21 移除以下三个口径（与 fig4_robustness 的 SPEC_DROPPED 一致）：
#   +log10 体积（双侧残差）／+NIHSS+mRS（双侧残差）／+NIHSS+mRS+HAMD（双侧残差）
#   理由：正文不披露卒中严重度（NIHSS/mRS），且这三个口径在 Fig. 4b 中不展示。
#   移除后 SPECS 为 4 个口径；另 2 个杠杆口径已于 2026-10-10 移出，
#   现与 derived_data/sensitivity_specs.csv 的 4 行一一对应。

print('\n=== 口径 2：单调秩变换 ===', flush=True)
tot = Y.sum(1)
Yrank = rankdata(tot).reshape(-1, 1)
w_rank = run(X, residualize(Yrank, [age, sex]))
SPECS.append(dict(spec='FSS-9 总分单调秩变换', group='rank', perm_p=w_rank[0],
                  corr=w_rank[1], n_sig=w_rank[2], expected=0.0444))
print('  %-38s p=%.4f corr=%+.4f n_sig=%3d' % ('FSS-9 总分单调秩变换', *w_rank[:3]), flush=True)
r_load = float(np.corrcoef(w_main / np.linalg.norm(w_main),
                           w_rank[3] / np.linalg.norm(w_rank[3]))[0, 1])
print('  主载荷 vs 秩变换载荷 Pearson r = %+.4f（正文口径 +0.981）' % r_load, flush=True)
pd.DataFrame(dict(region=parcels, loading_main=w_main, loading_rank=w_rank[3])).to_csv(
    os.path.join(HERE, 'fig3_loading_main_vs_rank.csv'), index=False, encoding='utf-8-sig')

# ---------------------------------------------------------------- 落盘 + 核对
df = pd.DataFrame(SPECS)
df.to_csv(os.path.join(HERE, 'fig3_sensitivity_specs.csv'), index=False, encoding='utf-8-sig')
bad = [s for s in SPECS if s['expected'] is not None and abs(s['perm_p'] - s['expected']) > 5e-4]
print('\n=== 与正文/记录值核对 ===')
for s in SPECS:
    if s['expected'] is None:
        continue
    ok = 'OK' if abs(s['perm_p'] - s['expected']) < 5e-4 else '**不符**'
    print('  %-38s 本脚本 %.4f | 记录 %.4f  %s' % (s['spec'], s['perm_p'], s['expected'], ok))
print('\n已写出 fig3_sensitivity_specs.csv / fig3_loading_main_vs_rank.csv')
print('不符项 %d 个' % len(bad))
