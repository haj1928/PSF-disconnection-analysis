# -*- coding: utf-8 -*-
"""
细胞层表达图与疲劳模式的关联及其特异性（新掩膜口径）
==================================================
判据：表达图的特异性以**偏相关**评估 —— 控制区域级断连负荷图后，重新检验该表达图与
    LV1 loading 的关联；其 spin p 由同一旋转零模型给出（与 PET 受体层的
    run_spin_partial_receptor.py 同口径）。

方法要点（与 run_sal_specificity_spin_Dataset5_*.py 完全同源）：
  - 表达图 = 细胞类基因集 median regional expression score（基因跨区 Z 后按区取 median）
  - 旋转 = Schaefer100 fsa5 球面旋转 1000×（seed=0, rotate_parcellation）
  - 同向单侧 p（按观测方向计数 exceed）；偏相关列为双向平均（方向A=旋转模式图、方向B=旋转细胞图，各取 p 后平均）
  - 同时报告各表达图与断连负荷图（特异性对照图）的空间相关

历史注：2026-09 曾一并计算「与两图空间相关之差」的配对置换检验；
    因两张模式图本身存在空间相关（r = -0.203），差值零分布的方差被放大约 1.7 倍、
    分辨力不足，该检验未用于稿件结论，已于 2026-09-22 自本脚本移除。

变更（2026-10-09）：① 坐标获取修正为行索引对齐（原 zip 写法系缺陷；官方语义见 enigmatoolbox centroid_extraction_sphere）；② 口径统一为同向单侧；③ r_partial 改为双向平均（方向A=旋转模式图；方向B=旋转细胞图；burden 均固定为实测图）；④ 方向B 对缺失区采用 0-填充后再旋转（修复 NaN 导致的 p 低估）。

输出: downstream_newmasks_20260915/celltype_specificity_contrast_20260915.csv
"""
import json
import sys

import nibabel as nb
import numpy as np
import pandas as pd

BASE = 'F:/dcmtonii/run62_crosssectional/'
PET = BASE + 'PSF_repro_PET'
EXPR = BASE + 'PSF_ahba_100/PSF_ahba_100_expression.csv'
SAL = BASE + 'PSF_pls_behavior/PSF_pls_discon_salience_residY_age_sex_cortex100_newmasks_20260915.csv'
CTRL = BASE + 'PSF_ahba_100/PSF_discon_control_vectors_newmasks_20260915.csv'
PANEL = BASE + 'PSF_references/Seidlitz2020/seidlitz2020_Dataset5_7class_ahba_aligned.json'
OUT = ('F:/dcmtonii/run62_crosssectional/PSF_pls_behavior/断连稿写作_20260907/修订_revision_20260910/'
       'downstream_newmasks_20260915/celltype_specificity_contrast_20260915.csv')
N_ROT, SEED = 1000, 0

# ---------- 1. 数据与旋转 ----------
expr = pd.read_csv(EXPR, index_col=0)
sal_df = pd.read_csv(SAL)
ctrl = pd.read_csv(CTRL).set_index('region')
PANEL_D = json.load(open(PANEL, encoding='utf-8'))
parcels = sal_df['region'].tolist()

sys.path.insert(0, PET)
from spin_test_numpy import rotate_parcellation  # noqa: E402

sphere_dir = f'{PET}/assets'
sp_lh = nb.load(f'{sphere_dir}/fsa5_sphere_lh.gii').darrays[0].data
sp_rh = nb.load(f'{sphere_dir}/fsa5_sphere_rh.gii').darrays[0].data


def coords_by_order(annot_f, sphere):
    labels, ctab, names = nb.freesurfer.io.read_annot(annot_f)
    names = [nm.decode() for nm in names]
    out = {}
    for i, nm in enumerate(names):
        if nm in ('unknown', 'medialwall'):
            continue
        m = labels == i
        if m.sum() == 0:
            continue
        out[nm.replace('7Networks_', '')] = sphere[m].mean(0)
    return out


coords_l = coords_by_order(f'{sphere_dir}/schaefer100_fsa5/lh.Schaefer2018_100Parcels_7Networks_order.annot', sp_lh)
coords_r = coords_by_order(f'{sphere_dir}/schaefer100_fsa5/rh.Schaefer2018_100Parcels_7Networks_order.annot', sp_rh)
cl, cr = [], []
for p in parcels:
    src = coords_l if p.startswith('LH') else coords_r
    (cl if p.startswith('LH') else cr).append(src[p])
perm_id = rotate_parcellation(np.array(cl), np.array(cr), nrot=N_ROT, seed=SEED).astype(int)

ez = (expr - expr.mean(0)) / (expr.std(0, ddof=1) + 1e-12)
keep = [i for i in range(100) if not expr.loc[parcels[i]].isna().all()]

sal_100 = sal_df.set_index('region')['salience'].reindex(parcels).values.astype(float)
bur_100 = ctrl['discon_mean'].reindex(parcels).values.astype(float)


def p_one(null, obs):
    null = np.asarray(null, dtype=float)
    return ((1 + int(np.sum(null >= obs))) / (1 + len(null)) if obs >= 0
            else (1 + int(np.sum(null <= obs))) / (1 + len(null)))


# ---------- 2. 两个模式本身的相关（混杂程度） ----------
sal_k, bur_k = sal_100[keep], bur_100[keep]
r_modes = float(np.corrcoef(sal_k, bur_k)[0, 1])
null_modes = np.array([np.corrcoef(sal_100[perm_id[:, k]][keep], bur_k)[0, 1] for k in range(N_ROT)])
p_modes = p_one(null_modes, r_modes)
print('=' * 96)
print('[0] 两个模式向量本身的相关（决定两检验的混杂程度）')
print('    r(LV1 loading, disconnection burden) = %+.4f | spin p = %.4f' % (r_modes, p_modes))
print('=' * 96)

# ---------- 2b. 工具：残差化与偏相关（双向平均口径，与 Table 3 / 表 S3 同源） ----------
def resid_cols(Y, Zcols):
    Zc = np.column_stack([np.ones(len(Zcols[0]))] + [z for z in Zcols])
    beta, *_ = np.linalg.lstsq(Zc, Y, rcond=None)
    return Y - Zc @ beta


def resid_rows(M, Zcols):
    Zc = np.column_stack([np.ones(len(Zcols[0]))] + [z for z in Zcols])
    beta, *_ = np.linalg.lstsq(Zc, M.T, rcond=None)
    return (M.T - Zc @ beta).T


def unitize(X):
    Xc = X - X.mean(1, keepdims=True)
    return Xc / (np.linalg.norm(Xc, axis=1, keepdims=True) + 1e-12)


# ---------- 3. 逐细胞检验（同向单侧；偏相关双向平均） ----------
KEEP = np.array(keep)
rows = []
for ct in ['Micro', 'Oligodendrocytes', 'Astro', 'Endo', 'Excitatory_Neurons', 'Inhibitory_Neurons', 'OPC']:
    gs = [g for g in PANEL_D.get(ct, []) if g in ez.columns]
    if not gs:
        continue
    cm = np.array([np.median(ez.loc[parcels[i], gs]) for i in range(100)])

    r_sal = float(np.corrcoef(cm[keep], sal_k)[0, 1])
    r_bur = float(np.corrcoef(cm[keep], bur_k)[0, 1])

    null_sal, null_bur = [], []
    for k in range(N_ROT):
        s_k = sal_100[perm_id[:, k]][keep]
        b_k = bur_100[perm_id[:, k]][keep]
        null_sal.append(np.corrcoef(cm[keep], s_k)[0, 1])
        null_bur.append(np.corrcoef(cm[keep], b_k)[0, 1])
    null_sal = np.array(null_sal)
    null_bur = np.array(null_bur)

    p_sal = p_one(null_sal, r_sal)
    p_bur = p_one(null_bur, r_bur)

    # 偏相关（双向平均 × 同向单侧）：方向A=旋转模式图（cm 固定）；方向B=旋转细胞图（模式图固定）
    tr = resid_cols(cm[KEEP][:, None], [bur_k])[:, 0]
    Sr = resid_cols(sal_k[:, None], [bur_k])[:, 0]
    part_obs = float(np.corrcoef(tr, Sr)[0, 1])
    Msal = sal_100[perm_id].T[:, KEEP]
    cm_f = cm.copy()
    cm_f[np.isnan(cm_f)] = 0.0          # 缺失区（RH_SomMot_2）0-填充：修复方向B的 NaN 污染（2026-10-09）
    Trot = cm_f[perm_id].T[:, KEEP]
    A = unitize(resid_rows(Msal, [bur_k])) @ unitize(tr[None, :])[0]
    Bv = unitize(resid_rows(Trot, [bur_k])) @ unitize(Sr[None, :])[0]
    p_part = (p_one(A, part_obs) + p_one(Bv, part_obs)) / 2

    rows.append(dict(cell_type=ct, n_genes=len(gs),
                     r_salience=round(r_sal, 4), p_salience=round(p_sal, 4),
                     r_burden=round(r_bur, 4), p_burden=round(p_bur, 4),
                     r_partial=round(float(part_obs), 4), p_partial=round(p_part, 4)))
    print('%-22s n=%-4d  r_sal=%+.4f(p=%.4f)  r_bur=%+.4f(p=%.4f)  '
          'r_partial=%+.4f(p=%.4f)'
          % (ct, len(gs), r_sal, p_sal, r_bur, p_bur, part_obs, p_part))

df = pd.DataFrame(rows)
df.to_csv(OUT, index=False, encoding='utf-8-sig')
print('\n已保存:', OUT)
print('\n判读规则：p_partial < 0.05 表示控制断连负荷图后，该表达图与 LV1 loading 的关联仍成立'
      '（表达图特异性的判据）。')
