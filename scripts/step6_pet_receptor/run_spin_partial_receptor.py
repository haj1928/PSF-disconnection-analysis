# -*- coding: utf-8 -*-
"""
R1：受体层偏相关的空间零模型校准（回应第二轮审稿意见 R1）
=====================================================================
问题：正文报 alpha4beta2 偏相关 r = +0.289, p = 0.0042。该 p 按常规自由度
      （df ~ 97）计算，隐含"100 个皮层区彼此空间独立"，而受体图谱与行为脑图
      均有强空间自相关，有效自由度远低于 97，故该 p 值系统性偏小。

做法：用与本分析主链完全相同的球面旋转零模型（Schaefer100 fsa5，
      1000 次旋转，seed = 0），重算三张显著图谱的偏相关零分布与 p 值。
      - 观测统计量：partial corr(atlas, LV1 loading | 其余两张显著图谱)
      - 零分布：分别旋转"焦点图谱"与"LV1 loading"，保持对照图谱固定
        （与 enigmatoolbox perm_sphere_p 的双向取平均约定一致；计数为含入式
        (1 + k)/(N + 1)，与论文数值口径一致）
      - 另报告：简单相关 r 的 spin p（复现原结果，作为管线校验）
      - 另报告：控制全部 21 张图谱第一主成分后的偏相关（回应 R7 共享成分）

输入：
  PSF_ahba_100/PSF_discon_receptor_mapping_cortex100_residY_newmasks_20260915.csv（原结果，校验用）
  PSF_pls_behavior/PSF_pls_discon_salience_residY_age_sex_cortex100_newmasks_20260915.csv（LV1 loading）
  PSF_repro_PET/data/neuromaps_annotations/...（21 张 PET 图谱）
  extdata/Schaefer_Yeo_Plus_Subcort/100Parcels7Networks.nii.gz（分区）
输出：
  report_spin_partial_receptor_newmasks_20260915.txt
  spin_partial_receptor_result_newmasks_20260915.csv
"""
import os
import sys
import numpy as np
import pandas as pd
import nibabel as nb
from scipy.ndimage import affine_transform

BASE = 'F:/dcmtonii/run62_crosssectional'
PET = BASE + '/PSF_repro_PET'
SAL_CSV = BASE + '/PSF_pls_behavior/PSF_pls_discon_salience_residY_age_sex_cortex100_newmasks_20260915.csv'
ATLAS = BASE + '/extdata/Schaefer_Yeo_Plus_Subcort/100Parcels7Networks.nii.gz'
NEUROMAPS = PET + '/data/neuromaps_annotations'
ORIG_CSV = BASE + '/PSF_ahba_100/PSF_discon_receptor_mapping_cortex100_residY_newmasks_20260915.csv'
OUT = os.path.dirname(os.path.abspath(__file__))
N_PERM = 1000
SEED = 0

log = []


def W(s=''):
    print(s, flush=True)
    log.append(s)


# ---------------- 1. LV1 loading ----------------
sal = pd.read_csv(SAL_CSV)
assert len(sal) == 100, f'须 100 皮层区, 实际 {len(sal)}'
lv1 = sal['salience'].values.astype(float)
parcels = sal['region'].tolist()
W(f'[1] LV1 loading（主分析版，皮层 100 区）范围 [{lv1.min():.4f}, {lv1.max():.4f}]')

# ---------------- 2. 受体图谱 -> Schaefer100 ----------------
atlas_img = nb.load(ATLAS)
atlas_data = atlas_img.get_fdata()
atlas_affine = atlas_img.affine


def parcellate_receptor(nii_path):
    img = nb.load(nii_path)
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

receptor_maps = {}
for src, tracer, nt, target, res in RECEPTOR_LIST:
    p = os.path.join(NEUROMAPS, src, tracer, 'MNI152',
                     f'source-{src}_desc-{tracer}_space-MNI152_res-{res}_feature.nii.gz')
    vec = parcellate_receptor(p)
    vec = np.where(np.isnan(vec), np.nanmean(vec), vec)
    receptor_maps[f'{src}_{tracer}'] = (vec, nt, target)
W(f'[2] 受体图谱加载 {len(receptor_maps)}/21')

# 管线校验：复算简单相关，应与已保存结果一致
orig = pd.read_csv(ORIG_CSV, encoding='utf-8-sig')
W('')
W('[3] 管线校验（复算简单相关 vs 已保存结果）')
chk = []
for _, row in orig.iterrows():
    vec = receptor_maps[row['receptor']][0]
    r_new = float(np.corrcoef(lv1, vec)[0, 1])
    chk.append(abs(r_new - row['r']))
W(f'    21 张 |r 复算 - 原值| 最大值 = {max(chk):.6f}  ->  {"一致" if max(chk) < 1e-3 else "不一致！"}')

# ---------------- 3. spin 旋转 ----------------
sys.path.insert(0, PET)
from spin_test_numpy import rotate_parcellation, perm_sphere_p

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
    hemi = p.split('_')[0]
    src = coords_l if hemi == 'LH' else coords_r
    (cl if hemi == 'LH' else cr).append(src[p])
perm_id = rotate_parcellation(np.array(cl), np.array(cr), nrot=N_PERM, seed=SEED)
assert perm_id.shape == (100, N_PERM)
W(f'[4] spin perm_id 就绪 {perm_id.shape} (seed = {SEED}, nrot = {N_PERM})')


def rotate(v, rr):
    return v[perm_id[:, rr].astype(int)]


# ---------------- 4. 偏相关工具 ----------------
def part_corr(x, y, Z):
    """x, y 对 [1, Z] 回归取残差后的 Pearson 相关。"""
    Zc = np.column_stack([np.ones(len(x)), Z]) if Z is not None else np.ones((len(x), 1))
    def res(v):
        b, _, _, _ = np.linalg.lstsq(Zc, v, rcond=None)
        return v - Zc @ b
    return float(np.corrcoef(res(x), res(y))[0, 1])


def p_bidir_inclusive(null_xy, null_yx, obs, n_perm):
    """含入式双向单侧 p：两方向各取 (1 + k)/(N + 1) 后平均（与论文数值口径一致）。

    官方 perm_sphere_p 的计数为严格式 k/N，两者相差约一格（≈0.001）。
    """
    def _one(v):
        if obs >= 0:
            return (1 + int((np.asarray(v) >= obs).sum())) / (n_perm + 1)
        return (1 + int((np.asarray(v) <= obs).sum())) / (n_perm + 1)

    return (_one(null_xy) + _one(null_yx)) / 2


def spin_p_partial(x, y, Z, obs, n_perm=N_PERM):
    """双向旋转零分布（旋转焦点图谱 / 旋转 LV1）；计数为含入式 (1 + k)/(N + 1) 双向平均。"""
    null_xy = np.array([part_corr(rotate(x, rr), y, Z) for rr in range(n_perm)])
    null_yx = np.array([part_corr(x, rotate(y, rr), Z) for rr in range(n_perm)])
    return p_bidir_inclusive(null_xy, null_yx, obs, n_perm), null_xy, null_yx


# ---------------- 5. 三张显著图谱的偏相关 ----------------
SIG = ['hillmer2016_flubatine', 'kantonen2020_carfentanil', 'norgaard2021_flumazenil']
LBL = {'hillmer2016_flubatine': 'a4b2-nACh', 'kantonen2020_carfentanil': 'mu-opioid',
       'norgaard2021_flumazenil': 'GABA-A'}

W('')
W('=' * 78)
W('[5] 三张显著图谱：简单相关（spin）与偏相关（常规 df / spin 校准）')
W('=' * 78)
W('')
W('| 图谱 | 简单 r | 简单 r spin p | 偏 r（控其余两张） | 偏 r 常规 p (df=97) | 偏 r spin p | 判读 |')
W('|---|---|---|---|---|---|---|')

rows = []
for name in SIG:
    x = receptor_maps[name][0]
    others = [receptor_maps[o][0] for o in SIG if o != name]
    Z = np.column_stack(others)

    r_simple = float(np.corrcoef(lv1, x)[0, 1])
    _, _nsp = perm_sphere_p(lv1, x, perm_id, corr_type='pearson', null_dist=True)
    p_simple_spin = p_bidir_inclusive(_nsp[:N_PERM], _nsp[N_PERM:], r_simple, N_PERM)

    obs = part_corr(x, lv1, Z)
    n = len(x)
    k = Z.shape[1] + 1  # 控制变量数 + 焦点变量
    df = n - k - 1
    t = obs * np.sqrt(df / max(1e-12, 1 - obs ** 2))
    from scipy import stats as st
    p_df = float(2 * (1 - st.t.cdf(abs(t), df)))

    p_spin, nxy, nyx = spin_p_partial(x, lv1, Z, obs)

    if p_spin < 0.05:
        verdict = '独立成立'
    elif obs == 0 or p_spin > 0.5:
        verdict = '不独立（被其余图谱吸收）'
    else:
        verdict = '边缘 / 不独立'
    rows.append({'atlas': LBL[name], 'receptor_key': name,
                 'r_simple': round(r_simple, 4), 'p_simple_spin': round(p_simple_spin, 4),
                 'partial_r': round(obs, 4), 'p_partial_parametric_df97': round(p_df, 6),
                 'p_partial_spin': round(p_spin, 4),
                 'null_xy_mean': round(float(nxy.mean()), 4), 'null_yx_mean': round(float(nyx.mean()), 4),
                 'verdict': verdict})
    W(f'| {LBL[name]} | {r_simple:+.3f} | {p_simple_spin:.4f} | **{obs:+.3f}** | {p_df:.4f} | '
      f'**{p_spin:.4f}** | {verdict} |')

# ---------------- 6. 控制全 21 张共享成分（R7）----------------
W('')
W('[6] 附加：控制 21 张图谱第一主成分后的偏相关（回应 R7「共享成分」）')
all_maps = np.column_stack([receptor_maps[k][0] for k in receptor_maps])
# 标准化后再做 SVD：与正文"21 张图谱 PC1 解释 41.7% 方差"（基于相关矩阵）口径一致
Zstd = (all_maps - all_maps.mean(0)) / all_maps.std(0, ddof=1)
U, S, Vt = np.linalg.svd(Zstd, full_matrices=False)
pc1 = U[:, 0] * S[0]
var1 = S[0] ** 2 / np.sum(S ** 2)
W(f'    21 张图谱（标准化后）PC1 解释 {100*var1:.1f}% 方差')
W('')
W('| 图谱 | 简单 r | 偏 r（控 PC1） | spin p |')
W('|---|---|---|---|')
for name in SIG:
    x = receptor_maps[name][0]
    obs = part_corr(x, lv1, pc1.reshape(-1, 1))
    p_spin, _, _ = spin_p_partial(x, lv1, pc1.reshape(-1, 1), obs)
    W(f'| {LBL[name]} | {float(np.corrcoef(lv1, x)[0,1]):+.3f} | {obs:+.3f} | {p_spin:.4f} |')
    rows.append({'atlas': LBL[name], 'receptor_key': name, 'r_simple': None, 'p_simple_spin': None,
                 'partial_r': round(obs, 4), 'p_partial_parametric_df97': None,
                 'p_partial_spin': round(p_spin, 4), 'null_xy_mean': None, 'null_yx_mean': None,
                 'verdict': f'控21张PC1（{100*var1:.1f}%方差）'})

pd.DataFrame(rows).to_csv(os.path.join(OUT, 'spin_partial_receptor_result_newmasks_20260915.csv'),
                          index=False, encoding='utf-8-sig')

with open(os.path.join(OUT, 'report_spin_partial_receptor_newmasks_20260915.txt'), 'w', encoding='utf-8') as f:
    f.write('\n'.join(log))
W('')
W(f'输出: report_spin_partial_receptor_newmasks_20260915.txt / spin_partial_receptor_result_newmasks_20260915.csv')
