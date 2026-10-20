# -*- coding: utf-8 -*-
"""
salience 空间置换检验（2026-09-02）——用户方案：固定表达矩阵，打乱/旋转 salience 权重
=====================================================================
背景：用户提出（正确）——salience 是 100 个脑区的权重值，把权重打乱或取平均，
      重跑 AHBA PLS 排序 + 细胞富集，若结果与真实 salience 不同 → 疲劳模式特异。
      这正是 Fulcher et al. 2021 (Nat Commun 12:2669) 检验的
      "random neural phenotypes" 思路（标准 GCEA 假阳性>500倍），
      以及 Alexander-Bloch et al. 2018 (NeuroImage 178:540) spin test
      （保留空间自相关的球面旋转置换）。

本脚本（第 1 步）：生成三类 salience 置换
  A. spin 置换（1000×, seed=0）——球面旋转保留空间自相关（最严格，标准做法）
  B. 简单置换（1000×, seed=0）——完全打乱区标签（更宽松的 null）
  C. 均匀权重（1×）——salience 全为均值（"权重取平均"= 无空间模式）
输出: PSF_ahba_100/permutation_newmasks_20260915/PSF_sal_perm_spin.csv (1000×100)
      PSF_ahba_100/permutation_newmasks_20260915/PSF_sal_perm_simple.csv (1000×100)
      PSF_ahba_100/permutation_newmasks_20260915/PSF_sal_uniform.csv (1×100)
"""
import os
import sys
import numpy as np
import pandas as pd
import nibabel as nb

BASE = 'F:/dcmtonii/run62_crosssectional'
PET = BASE + '/PSF_repro_PET'
SAL_CSV = BASE + '/PSF_pls_behavior/PSF_pls_discon_salience_residY_age_sex_cortex100_newmasks_20260915.csv'
OUT = BASE + '/PSF_ahba_100/permutation_newmasks_20260915'
N_PERM = 1000
SEED = 0
os.makedirs(OUT, exist_ok=True)

# ---- 1. salience（100 区, 顺序=id 1-100）----
sal = pd.read_csv(SAL_CSV)
salience = sal['salience'].values.astype(float)
parcels = sal['region'].tolist()
print(f'[1] salience 100 区 OK | 范围 [{salience.min():.4f}, {salience.max():.4f}]', flush=True)

# ---- 2. spin 置换（复用受体分析同一套 Schaefer100 fsa5 旋转）----
sys.path.insert(0, PET)
from spin_test_numpy import rotate_parcellation

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
coord_l_list, coord_r_list = [], []
for p in parcels:
    hemi = p.split('_')[0]
    src = coords_l if hemi == 'LH' else coords_r
    assert p in src, f'皮层区 {p} 质心未对齐'
    (coord_l_list if hemi == 'LH' else coord_r_list).append(src[p])
perm_id = rotate_parcellation(np.array(coord_l_list), np.array(coord_r_list),
                              nrot=N_PERM, seed=SEED)
assert perm_id.shape == (100, N_PERM)
perm_id = perm_id.astype(int)   # rotate_parcellation 可能返回 float 索引
print(f'[2] spin perm_id 就绪: {perm_id.shape} (seed={SEED})', flush=True)

# ---- 3. 应用置换到 salience ----
sal_perm_spin = salience[perm_id].T          # 1000 × 100（每行=一次旋转后的 salience）
sal_perm_spin_df = pd.DataFrame(sal_perm_spin, columns=parcels)
sal_perm_spin_df.to_csv(f'{OUT}/PSF_sal_perm_spin.csv', index=False, encoding='utf-8-sig')

# 简单置换（完全打乱区标签）
rng = np.random.default_rng(SEED)
sal_perm_simple = np.array([salience[rng.permutation(100)] for _ in range(N_PERM)])
pd.DataFrame(sal_perm_simple, columns=parcels).to_csv(
    f'{OUT}/PSF_sal_perm_simple.csv', index=False, encoding='utf-8-sig')

# 均匀权重（取平均 = 所有区相同权重 = 无空间模式）
sal_uniform = np.full(100, salience.mean())
pd.DataFrame([sal_uniform], columns=parcels).to_csv(
    f'{OUT}/PSF_sal_uniform.csv', index=False, encoding='utf-8-sig')

# ---- 4. 质量检查 ----
print(f'[3] spin 置换: 每行均值 {sal_perm_spin.mean(1).mean():.4f} (应≈{salience.mean():.4f})', flush=True)
r_obs = np.corrcoef(salience, sal_perm_spin[0])[0, 1]
print(f'    第 1 次 spin 与真实 salience 的 corr = {r_obs:+.3f} (应≈0, 置换有效)', flush=True)
r_obs2 = np.corrcoef(salience, sal_perm_simple[0])[0, 1]
print(f'    第 1 次 simple 与真实 salience 的 corr = {r_obs2:+.3f} (应≈0)', flush=True)
print(f'    均匀权重: 全 = {sal_uniform[0]:.4f} (无空间模式)', flush=True)
print('完成。')
