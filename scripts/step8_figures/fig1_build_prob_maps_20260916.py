# -*- coding: utf-8 -*-
"""
fig1_build_prob_maps_20260916.py —— Fig. 1 数据层：按病灶体积三分位构建组内病灶重叠概率图。

输入（数据，不硬编码）
  * 病灶体积表：修订_revision_20260910/PSF_lesion_volume_N124_newmasks_20260915.csv
    （含 subject / volume_mL / path，路径即 A1 新口径重跑后的 MNI 掩码）
输出
  * fig1_prob_tertile{1,2,3}.nii.gz —— 各三分位的组内重叠概率（0–1）
  * fig1_tertile_summary.csv       —— 分层切点、n、最大概率、非零体素统计（供出图脚本读取）
"""
import csv
import io
import os

import nibabel as nib
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = 'F:/dcmtonii/run62_crosssectional/'
VOL = (ROOT + 'PSF_pls_behavior/断连稿写作_20260907/修订_revision_20260910/'
       'PSF_lesion_volume_N124_newmasks_20260915.csv')

# ---------------------------------------------------------------- 1. 读数据
rows = list(csv.DictReader(io.open(VOL, encoding='utf-8-sig')))
vols = np.array([float(x['volume_mL']) for x in rows], dtype=float)
n_all = len(rows)
assert n_all == 124, '样本量异常：%d' % n_all

# 等例数三分位（41 / 42 / 41），与「约 33%/33%/33%」一致
order = np.argsort(vols, kind='stable')
grp = np.empty(n_all, dtype=int)
grp[order[:41]] = 0
grp[order[41:83]] = 1
grp[order[83:]] = 2

# ---------------------------------------------------------------- 2. 叠加
stack = {g: None for g in (0, 1, 2)}
ref_img = None
for i, x in enumerate(rows):
    img = nib.load(x['path'].replace('\\', '/'))
    if ref_img is None:
        ref_img = img
    k = (np.asarray(img.dataobj) > 0.5).astype(np.float32)
    g = grp[i]
    stack[g] = k if stack[g] is None else stack[g] + k

# ---------------------------------------------------------------- 3. 概率 + 落盘
summary = []
gmaps = {}
for g, (name_cn, name_en) in enumerate([('小', 'small'), ('中', 'medium'), ('大', 'large')]):
    n_g = int((grp == g).sum())
    m = stack[g]
    prob = m / float(n_g)
    gmaps[g] = prob
    out = os.path.join(HERE, 'fig1_prob_tertile%d.nii.gz' % (g + 1))
    nib.save(nib.Nifti1Image(prob.astype(np.float32), ref_img.affine, ref_img.header), out)
    nz = prob[prob > 0]
    summary.append(dict(
        tertile=g + 1, group=name_en, label_cn=name_cn, n=n_g,
        vol_min_mL='%.4f' % vols[order[[0, 41, 83][g]]],
        vol_max_mL='%.4f' % vols[order[[40, 82, 123][g]]],
        vol_median_mL='%.4f' % np.median(vols[grp == g]),
        max_overlap_n=int(m.max()),
        max_prob='%.6f' % prob.max(),
        nonzero_voxels=int(nz.size),
        prob_median_nonzero='%.6f' % np.median(nz),
        prob_p90_nonzero='%.6f' % np.percentile(nz, 90),
        prob_p99_nonzero='%.6f' % np.percentile(nz, 99),
        out_nii=os.path.basename(out),
    ))
    print('tertile %d (%s) n=%d  体积 %.3f–%.3f mL  最大重叠 %d 例 → 概率 %.4f  非零体素 %d'
          % (g + 1, name_en, n_g, vols[grp == g].min(), vols[grp == g].max(),
             int(m.max()), prob.max(), nz.size))

# 全局最大概率（决定共享色标上限）
gmax = max(gmaps[g].max() for g in gmaps)
print('\n全局最大概率 = %.6f' % gmax)

# 层内概率图（各层自身最大）与全局的关系
for g in gmaps:
    print('  tertile %d 峰值 %.4f  = 全局峰值的 %.1f%%' % (g + 1, gmaps[g].max(), 100 * gmaps[g].max() / gmax))

# 掩码与模板几何
aff = np.asarray(ref_img.affine)
print('\n掩码几何: shape=%s  diag=%s  trans=%s'
      % (ref_img.shape, np.diag(aff)[:3].round(3).tolist(), aff[:3, 3].round(3).tolist()))

# 脑范围（用全体掩码并集近似不可行，改用各层面解剖由模板决定）——输出病灶包围盒
miss = np.zeros(ref_img.shape, dtype=bool)
for g in gmaps:
    miss |= gmaps[g] > 0
idx = np.where(miss)
print('病灶包围盒（体素索引）x %d–%d  y %d–%d  z %d–%d'
      % (idx[0].min(), idx[0].max(), idx[1].min(), idx[1].max(), idx[2].min(), idx[2].max()))
print('病灶包围盒（MNI mm）  x %.0f–%.0f  y %.0f–%.0f  z %.0f–%.0f'
      % (90 - idx[0].max(), 90 - idx[0].min(),
         -126 + idx[1].min(), -126 + idx[1].max(),
         -72 + idx[2].min(), -72 + idx[2].max()))

out_csv = os.path.join(HERE, 'fig1_tertile_summary.csv')
with io.open(out_csv, 'w', encoding='utf-8-sig', newline='') as f:
    w = csv.DictWriter(f, fieldnames=list(summary[0].keys()))
    w.writeheader()
    w.writerows(summary)
print('\n已写出 %s' % out_csv)

# 逐例分层归属（可溯源）
assign = os.path.join(HERE, 'fig1_tertile_assignment.csv')
with io.open(assign, 'w', encoding='utf-8-sig', newline='') as f:
    w = csv.writer(f)
    w.writerow(['subject', 'volume_mL', 'tertile', 'group'])
    names = ['small', 'medium', 'large']
    for i, x in enumerate(rows):
        w.writerow([x['subject'], x['volume_mL'], grp[i] + 1, names[grp[i]]])
print('已写出 %s' % assign)
print('已写出 fig1_prob_tertile1/2/3.nii.gz')
