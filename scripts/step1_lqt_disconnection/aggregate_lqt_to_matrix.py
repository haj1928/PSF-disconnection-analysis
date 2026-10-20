# -*- coding: utf-8 -*-
"""
code_aggregate_lqt_to_matrix_20260915.py

把 LQT 逐被试输出汇总为 124×135 区域断连矩阵（口径与既有主矩阵完全一致），
并支持"61 例复用旧输出 + 63 例用新掩膜重跑输出"的混合汇总。

汇总口径（取自 lqt_matlab_full_N124.m 注释 + verify_matlab_vs_r_lqt.py，已逐值验证过）：
    pat_raw       = pct_sdc_matrix / 100 .* atlas_con          （逐元素）
    parc_discon[i]= sum_j pat_raw[i,j] / sum_j atlas_con[i,j]  （行内占比）
其中 pct_sdc_matrix 来自 <root>/<sub>/Parcel_Disconnection/<sub>_Yeo7100_percent_parcel_SDC.mat，
atlas_con 来自      <root>/Atlas/atlas_Yeo7100.*.count.end.connectivity.mat（变量 connectivity）。

模式
----
--verify : 用全部旧输出（LQT_MATLAB_Full）复现既有主矩阵 PSF_SDC_parcel_discon_N124.csv，
           逐值比对以证明本脚本口径正确。
--build  : 61 例（未变掩膜）取旧输出 + 63 例取 LQT_MATLAB_Full_v2_20260915，
           写出新矩阵与皮层 100 版。
"""
import argparse
import io
import os
import glob
import warnings

import numpy as np
import pandas as pd
from scipy import io

R = 'F:/dcmtonii'
RUN = os.path.join(R, 'run62_crosssectional')
OUT_DIR = os.path.join(RUN, 'PSF_pls_behavior', '断连稿写作_20260907', '修订_revision_20260910')
MAIN_FULL = os.path.join(R, 'PSF_SDC_parcel_discon_N124.csv')
MAIN_C100 = os.path.join(R, 'PSF_SDC_parcel_discon_N124_cortex100.csv')
OLD_ROOT = os.path.join(RUN, 'LQT_MATLAB_Full')
NEW_ROOT = os.path.join(RUN, 'LQT_MATLAB_Full_v2_20260915')
PLAN = os.path.join(OUT_DIR, 'mask_switch_plan_20260915.csv')

warnings.filterwarnings('ignore')


def load_atlas(root):
    f = glob.glob(os.path.join(root, 'Atlas', '*connectivity.mat'))
    assert len(f) == 1, 'Atlas connectivity 文件不唯一: %s' % f
    m = io.loadmat(f[0])
    return np.asarray(m['connectivity'], dtype=float)


def load_pct(root, sub):
    f = glob.glob(os.path.join(root, sub, 'Parcel_Disconnection', '*percent_parcel_SDC.mat'))
    if not f:
        return None
    m = io.loadmat(f[0])
    k = [x for x in m if not x.startswith('__')][0]
    return np.asarray(m[k], dtype=float)


def parc_discon(pct, atlas_con):
    """逐元素重建 raw 断连，再取行内占比（与 R 版 pdis 同口径）"""
    pat_raw = pct / 100.0 * atlas_con
    num = pat_raw.sum(axis=1)
    den = atlas_con.sum(axis=1)
    with np.errstate(invalid='ignore', divide='ignore'):
        out = num / den
    out[den == 0] = np.nan
    return out


def aggregate(plan, root_of):
    """plan: DataFrame(subject,...)；root_of: 函数 subject -> root 路径"""
    cols = list(pd.read_csv(MAIN_FULL, nrows=0).columns)      # ['ID', 135 区域名]
    regions = cols[1:]
    atlas_cache = {}
    rows, missing = [], []
    for s in plan.subject:
        root = root_of(s)
        if root not in atlas_cache:
            atlas_cache[root] = load_atlas(root)
        pct = load_pct(root, s)
        if pct is None:
            missing.append(s)
            rows.append([np.nan] * len(regions))
            continue
        rows.append(parc_discon(pct, atlas_cache[root]))
    X = np.vstack([np.asarray(r, dtype=float) if np.ndim(r) else np.full(len(regions), np.nan) for r in rows])
    return pd.DataFrame(np.c_[np.array(plan.subject, dtype=object), X], columns=cols), missing


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--verify', action='store_true')
    ap.add_argument('--build', action='store_true')
    a = ap.parse_args()

    plan = pd.read_csv(PLAN)
    if a.verify:
        df, miss = aggregate(plan, lambda s: OLD_ROOT)
        ref = pd.read_csv(MAIN_FULL)
        assert list(df.columns) == list(ref.columns)
        assert list(df['ID'].astype(str)) == list(ref['ID'].astype(str))
        A = df.drop(columns=['ID']).values.astype(float)
        B = ref.drop(columns=['ID']).values.astype(float)
        m = ~(np.isnan(A) | np.isnan(B))
        print('[验证] 复用位置一致的旧输出复现主矩阵')
        print('  shape %s vs %s | 可比元素 %d | NaN 位置一致 %s'
              % (A.shape, B.shape, int(m.sum()), bool(np.array_equal(np.isnan(A), np.isnan(B)))))
        print('  max|diff| = %.3e | 相关 r = %.10f'
              % (np.abs(A[m] - B[m]).max(),
                 np.corrcoef(A[m], B[m])[0, 1]))
        if miss:
            print('  缺输出被试:', miss)
        ok = np.abs(A[m] - B[m]).max() < 1e-10
        print('  → 口径复现%s' % ('成功 ✓' if ok else '失败 ✗'))
        return

    if a.build:
        def root_of(s):
            reuse = bool(plan.set_index('subject').loc[s, 'reuse_lqt'])
            return OLD_ROOT if reuse else NEW_ROOT
        df, miss = aggregate(plan, root_of)
        assert not miss, '以下被试在对应 LQT 目录缺输出: %s' % miss
        out_full = os.path.join(R, 'PSF_SDC_parcel_discon_N124_newmasks_20260915.csv')
        df.to_csv(out_full, index=False, encoding='utf-8-sig')
        cort = [c for c in df.columns if c.startswith('LH_') or c.startswith('RH_')]
        out_c100 = os.path.join(R, 'PSF_SDC_parcel_discon_N124_cortex100_newmasks_20260915.csv')
        df[['ID'] + cort].to_csv(out_c100, index=False, encoding='utf-8-sig')
        print('[构建] 新矩阵已写出')
        print('  全脑版 : %s  shape %s' % (out_full, df.shape))
        print('  皮层100: %s  shape %s' % (out_c100, (len(df), len(cort) + 1)))
        # 与主矩阵对比
        ref = pd.read_csv(MAIN_FULL).drop(columns=['ID']).values.astype(float)
        A = df.drop(columns=['ID']).values.astype(float)
        m = ~(np.isnan(A) | np.isnan(ref))
        r_row = [np.corrcoef(A[i][m[i]], ref[i][m[i]])[0, 1] for i in range(len(A))]
        print()
        print('  与主矩阵逐被试相关：中位 %.4f | 最小 %.4f | >0.999 的比例 %.1f%%'
              % (np.median(r_row), np.min(r_row), 100 * np.mean(np.array(r_row) > 0.999)))
        return

    ap.print_help()


if __name__ == '__main__':
    main()
