# -*- coding: utf-8 -*-
"""
Step 8'（皮层 100 区·主分析版）：断连 salience × PET 受体图谱（21 张）空间相关
=====================================================================
背景（2026-09-01）：主分析定稿为"仅残差 Y(age+sex)"皮层 100 区 salience。
      Step8（run_receptor_discon.py）已是 Schaefer100 皮层 100 区口径，
      仅需把 salience 输入从旧版（PSF_pls_discon_parcel_salience.csv）换成
      新主分析版（PSF_pls_discon_salience_residY_age_sex_cortex100_newmasks_20260915.csv）。

方法（与旧版完全一致，仅换 salience）：
  - salience_100 : 新主分析版（region 顺序 = id 1-100, LH 1-50 / RH 51-100）
  - receptor_100 : 21 张 neuromaps PET 图谱（MNI152, 1/2/3mm）→ Schaefer100 分区 nanmean
  - r : Pearson | spin : Schaefer100 fsa5 球面旋转 1000× seed=0 双向平均，含入式计数 (1+k)/(N+1) | FDR : 21 张组内 BH
输出: PSF_ahba_100/PSF_discon_receptor_mapping_cortex100_residY_newmasks_20260915.csv + report
"""
import os, sys
import numpy as np, pandas as pd
import nibabel as nb
from scipy.ndimage import affine_transform

BASE = 'F:/dcmtonii/run62_crosssectional'
PET = BASE + '/PSF_repro_PET'
SAL_CSV = BASE + '/PSF_pls_behavior/PSF_pls_discon_salience_residY_age_sex_cortex100_newmasks_20260915.csv'
ATLAS = BASE + '/extdata/Schaefer_Yeo_Plus_Subcort/100Parcels7Networks.nii.gz'
NEUROMAPS = PET + '/data/neuromaps_annotations'
OUT = BASE + '/PSF_ahba_100'
N_PERM = 1000
SEED = 0

# ============ 1. salience（新主分析版：100 行顺序 = id 1-100）============
sal = pd.read_csv(SAL_CSV)
assert len(sal) == 100, f'须 100 皮层区, 实际 {len(sal)}'
salience = sal['salience'].values.astype(float)
parcels = sal['region'].tolist()
print(f'[1] salience(residY) 皮层 100 区 OK | 范围 [{salience.min():.4f}, {salience.max():.4f}]', flush=True)

# ============ 2. 受体图谱 → Schaefer100 分区 ============
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
    if not os.path.exists(p):
        print(f'  [err] 缺文件: {p}'); continue
    vec = parcellate_receptor(p)
    n_nan = int(np.isnan(vec).sum())
    vec = np.where(np.isnan(vec), np.nanmean(vec), vec)
    receptor_maps[f'{src}_{tracer}'] = vec
    print(f'  [ok] {src}_{tracer:24s} {nt}/{target}  NaN区={n_nan}', flush=True)
print(f'[2] 受体图谱加载 {len(receptor_maps)}/{len(RECEPTOR_LIST)}', flush=True)

# ============ 3. spin 置换（Schaefer100 fsa5, seed=0）============
sys.path.insert(0, PET)
from spin_test_numpy import rotate_parcellation, perm_sphere_p


def p_bidir_inclusive(null_both, obs, n_perm):
    """含入式双向单侧 p：两方向各取 (1 + k)/(N + 1) 后平均（k = 同向更极端的旋转数）。

    与论文数值口径一致（观测计入置换分布）；官方 perm_sphere_p 的计数为严格式 k/N，
    两者相差约一格（≈0.001，如 0.1670 严格 vs 0.1678 含入）。零模型仍由官方函数生成
    （null_dist=True 依序返回 [方向1 的 N 个 null, 方向2 的 N 个 null]）。
    """
    a = np.asarray(null_both[:n_perm], dtype=float)
    b = np.asarray(null_both[n_perm:], dtype=float)

    def _one(v):
        if obs >= 0:
            return (1 + int((v >= obs).sum())) / (n_perm + 1)
        return (1 + int((v <= obs).sum())) / (n_perm + 1)

    return (_one(a) + _one(b)) / 2

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
print(f'[3] spin perm_id 就绪: {perm_id.shape} (seed={SEED})', flush=True)

# ============ 4. r + spin p ============
results = []
for name, vec in receptor_maps.items():
    r = float(np.corrcoef(salience, vec)[0, 1])
    _, null = perm_sphere_p(salience, vec, perm_id, corr_type='pearson', null_dist=True)
    p = p_bidir_inclusive(null, r, N_PERM)   # 含入式双向单侧（与论文口径一致）
    meta = [m for m in RECEPTOR_LIST if f'{m[0]}_{m[1]}' == name][0]
    results.append({'receptor': name, 'neurotransmitter': meta[2], 'target': meta[3],
                    'r': round(r, 4), 'spin_p': round(float(p), 4)})
    sig = '***' if p < 0.001 else ('**' if p < 0.01 else ('*' if p < 0.05 else 'ns'))
    print(f'  {meta[3]:15s} r={r:+.4f} spin_p={p:.4f} {sig}', flush=True)

# ============ 5. BH-FDR（21 张组内）============
df = pd.DataFrame(results)
pv = df['spin_p'].values
order = np.argsort(pv)
fdr = pv[order] * len(pv) / (np.arange(1, len(pv) + 1))
fdr = np.minimum.accumulate(fdr[::-1])[::-1]
df['fdr'] = np.nan
df.loc[df.index[order], 'fdr'] = np.clip(fdr, 0, 1)
df['fdr'] = df['fdr'].round(4)
df = df.sort_values('spin_p').reset_index(drop=True)
df.to_csv(f'{OUT}/PSF_discon_receptor_mapping_cortex100_residY_newmasks_20260915.csv', index=False, encoding='utf-8-sig')

# ============ 6. 报告 ============
with open(f'{OUT}/PSF_discon_receptor_report_cortex100_residY_newmasks_20260915.txt', 'w', encoding='utf-8') as f:
    f.write('Step8\' 断连 salience(residY) × PET 受体/转运体图谱（21 张, Schaefer100 皮层 100 区）\n')
    f.write(f'salience: 新主分析版（仅残差 age+sex）| spin: Schaefer100 fsa5 球面旋转 {N_PERM}× seed={SEED} | BH-FDR 组内\n\n')
    for _, r in df.iterrows():
        mark = '★' if r['fdr'] < 0.05 else ('·' if r['spin_p'] < 0.05 else '')
        f.write(f"  {r['target']:15s} {r['neurotransmitter']:<12s} r={r['r']:+.4f} spin_p={r['spin_p']:.4f} FDR={r['fdr']:.4f} {mark}\n")
    n_sig_spin = int((df['spin_p'] < 0.05).sum())
    n_sig_fdr = int((df['fdr'] < 0.05).sum())
    f.write(f'\nspin p<0.05: {n_sig_spin}/21 | FDR<0.05: {n_sig_fdr}/21\n')
print(f'\n输出: PSF_discon_receptor_mapping_cortex100_residY_newmasks_20260915.csv (+report)')
print(f'spin p<0.05: {(df["spin_p"]<0.05).sum()}/21 | FDR<0.05: {(df["fdr"]<0.05).sum()}/21')
