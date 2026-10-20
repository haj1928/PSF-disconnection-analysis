# -*- coding: utf-8 -*-
"""
补算 3：21 张 PET 受体/转运体图谱的 Schaefer100 皮层分区值矩阵（21 × 100）

用途
  * 【已移出稿件】原补充图「21 张受体图谱的皮层分布」（7×3 表面图网格）
  * Fig. 8a 需要图谱名与靶点元数据
方法
  与主分析 `run_receptor_discon_cortex100_residY_newmasks_20260915.py` **同一套**
  图谱清单、同一 parcellate_receptor 实现（含 affine 对齐 + 区内 nanmean），
  仅**不跑 spin**、只把分区值落盘。
输出
  unused_receptor_parcel_values.csv  —— 21 行（图谱）× 100 列（区）+ 元数据列
"""
import os

import nibabel as nib
import numpy as np
import pandas as pd
from scipy.ndimage import affine_transform

BASE = 'F:/dcmtonii/run62_crosssectional'
PET = BASE + '/PSF_repro_PET'
ATLAS = BASE + '/extdata/Schaefer_Yeo_Plus_Subcort/100Parcels7Networks.nii.gz'
NEUROMAPS = PET + '/data/neuromaps_annotations'
SAL_CSV = (BASE + '/PSF_pls_behavior/PSF_pls_discon_salience_residY_age_sex_cortex100'
           '_newmasks_20260915.csv')
HERE = os.path.dirname(os.path.abspath(__file__))

sal = pd.read_csv(SAL_CSV)
parcels = sal['region'].tolist()
assert len(parcels) == 100

atlas_img = nib.load(ATLAS)
atlas_data = atlas_img.get_fdata()
atlas_affine = atlas_img.affine


def parcellate_receptor(nii_path):
    img = nib.load(nii_path)
    rec = img.get_fdata().squeeze()
    if rec.shape != atlas_data.shape or not np.allclose(img.affine, atlas_affine):
        comb = np.linalg.inv(atlas_affine) @ img.affine
        atlas_r = affine_transform(atlas_data, comb[:3, :3], offset=comb[:3, 3],
                                   output_shape=rec.shape, order=0).astype(int)
    else:
        atlas_r = atlas_data.astype(int)
    vals = np.full(100, np.nan)
    for lab in range(1, 101):
        m = (atlas_r == lab) & ~np.isnan(rec)
        if m.sum() > 0:
            vals[lab - 1] = np.nanmean(rec[m])
    return vals


RECEPTOR_LIST = [
    ('beliveau2017', 'dasb', 'Serotonin', '5-HTT', '1mm'),
    ('beliveau2017', 'cumi101', 'Serotonin', '5-HT1A', '1mm'),
    ('beliveau2017', 'az10419369', 'Serotonin', '5-HT1B', '1mm'),
    ('beliveau2017', 'cimbi36', 'Serotonin', '5-HT2A', '1mm'),
    ('beliveau2017', 'sb207145', 'Serotonin', '5-HT4', '1mm'),
    ('savli2012', 'way100635', 'Serotonin', '5-HT1A(2)', '3mm'),
    ('savli2012', 'altanserin', 'Serotonin', '5-HT2A(2)', '3mm'),
    ('kaller2017', 'sch23390', 'Dopamine', 'D1', '3mm'),
    ('jaworska2020', 'fallypride', 'Dopamine', 'D2/D3', '1mm'),
    ('dukart2018', 'fpcit', 'Dopamine', 'DAT', '3mm'),
    ('sandiego2015', 'flb457', 'Dopamine', 'D2/D3(2)', '1mm'),
    ('norgaard2021', 'flumazenil', 'GABA', 'GABA-A', '1mm'),
    ('lukow2022', 'ro154513', 'GABA', 'GABA-A-a5', '2mm'),
    ('dubois2015', 'abp688', 'Glutamate', 'mGluR5', '1mm'),
    ('galovic2021', 'ge179', 'Glutamate', 'NMDA', '1mm'),
    ('aghourian2017', 'feobv', 'Acetylcholine', 'VAChT', '1mm'),
    ('hillmer2016', 'flubatine', 'Acetylcholine', 'a4b2-nACh', '1mm'),
    ('naganawa2020', 'lsn3172176', 'Acetylcholine', 'M1 mAChR', '1mm'),
    ('hesse2017', 'methylreboxetine', 'Noradrenaline', 'NET', '3mm'),
    ('kantonen2020', 'carfentanil', 'Opioid', 'mu-opioid', '3mm'),
    ('finnema2016', 'ucbj', 'Synaptic', 'SV2A', '1mm'),
]

rows = []
for src, tracer, nt, target, res in RECEPTOR_LIST:
    p = os.path.join(NEUROMAPS, src, tracer, 'MNI152',
                     'source-%s_desc-%s_space-MNI152_res-%s_feature.nii.gz'
                     % (src, tracer, res))
    if not os.path.exists(p):
        print('  [err] 缺文件: %s' % p)
        continue
    vec = parcellate_receptor(p)
    n_nan = int(np.isnan(vec).sum())
    if n_nan:
        vec = np.where(np.isnan(vec), np.nanmean(vec), vec)
    row = {'atlas': '%s_%s' % (src, tracer), 'neurotransmitter': nt, 'target': target,
           'res': res, 'n_nan_filled': n_nan}
    row.update({'%s' % pc: vec[i] for i, pc in enumerate(parcels)})
    rows.append(row)
    print('  [ok] %-26s %-14s %-12s NaN 填充 %d' % (row['atlas'], nt, target, n_nan))

df = pd.DataFrame(rows)
out = os.path.join(HERE, 'unused_receptor_parcel_values.csv')
df.to_csv(out, index=False, encoding='utf-8-sig')
print('\n图谱 %d / 21 | 形状 %s' % (len(df), df.shape))
print('已写出 %s' % out)

# 自检：与主分析输出表（r 值）核对图谱名集合是否一致
main = pd.read_csv(BASE + '/PSF_ahba_100/'
                   'PSF_discon_receptor_mapping_cortex100_residY_newmasks_20260915.csv')
same = set(df['atlas']) == set(main['receptor'])
print('图谱名与主分析输出一致：%s' % ('PASS' if same else 'FAIL'))
if not same:
    print('  仅本表有：', set(df['atlas']) - set(main['receptor']))
    print('  仅主表有：', set(main['receptor']) - set(df['atlas']))
# z 标准化（供表面图使用）后再存一份
z = df.copy()
cols = [c for c in parcels]
z[cols] = (df[cols] - df[cols].mean(0)) / df[cols].std(0)
z.to_csv(os.path.join(HERE, 'unused_receptor_parcel_values_z.csv'), index=False,
         encoding='utf-8-sig')
print('已写出 unused_receptor_parcel_values_z.csv（区内 z 标准化，供表面渲染）')
