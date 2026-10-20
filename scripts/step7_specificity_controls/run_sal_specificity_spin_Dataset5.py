# -*- coding: utf-8 -*-
"""
细胞表达图 × {salience, disconnection burden} 空间相关 + spin（Dataset5 官方面板版）
=====================================================================
背景（2026-09-02）：Step7 表 A 原为 63 marker 面板时代产物（Micro +0.312/p=0.002 等）。
      细胞富集主口径切换为 Seidlitz Dataset5 面板后，空间层检验须按新面板重算，
      并明确列出 Exci/Inhi（Dataset5 下 PRAND 显著，需空间层交叉验证）。
方法（与 run_sal_specificity_spin_v2.py 同款，仅面板不同）：
  - 表达图 = 细胞类基因集 median regional expression score（Seidlitz 2020 用法；
    基因跨区 Z 后按脑区取基因集 median）
  - 相关 = Pearson(表达图, 模式向量)（100 区，99 有效区）
  - spin = Schaefer100 fsa5 球面旋转 1000×（seed=0, rotate_parcellation,
    spin_test_numpy.py 与受体层同一套）→ P_SPIN（双侧）
  - 模式 = salience（疲劳, 主分析）+ disconnection burden（断连负荷对照）
输入: salience_residY(100区) + PSF_discon_control_vectors_newmasks_20260915.csv + Dataset5 面板 JSON
输出: PSF_ahba_100/permutation_newmasks_20260915/PSF_sal_specificity_spin_Dataset5_result_newmasks_20260915.csv
"""
import pandas as pd, numpy as np, json, sys
import nibabel as nb

BASE = 'F:/dcmtonii/run62_crosssectional/'
PET = BASE + 'PSF_repro_PET'
EXPR = BASE + 'PSF_ahba_100/PSF_ahba_100_expression.csv'
SAL = BASE + 'PSF_pls_behavior/PSF_pls_discon_salience_residY_age_sex_cortex100_newmasks_20260915.csv'
CTRL = BASE + 'PSF_ahba_100/PSF_discon_control_vectors_newmasks_20260915.csv'
PANEL = BASE + 'PSF_references/Seidlitz2020/seidlitz2020_Dataset5_7class_ahba_aligned.json'
OUT = BASE + 'PSF_ahba_100/permutation_newmasks_20260915/PSF_sal_specificity_spin_Dataset5_result_newmasks_20260915.csv'
N_ROT = 1000
SEED = 0

# ---- 1. 数据 ----
expr = pd.read_csv(EXPR, index_col=0)
sal_df = pd.read_csv(SAL)
ctrl = pd.read_csv(CTRL).set_index('region')
PANEL_D = json.load(open(PANEL, encoding='utf-8'))
parcels = sal_df['region'].tolist()

# spin perm_id（Schaefer100 fsa5, seed=0, 与受体层同一套）
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
        if nm in ('unknown', 'medialwall'): continue
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
print(f'[1] spin perm_id: {perm_id.shape} (seed={SEED}, nrot={N_ROT})', flush=True)

# ---- 2. 表达图（median 聚合）+ 模式 ----
ez = (expr - expr.mean(0)) / (expr.std(0, ddof=1) + 1e-12)
keep = [i for i in range(100) if not expr.loc[parcels[i]].isna().all()]
print(f'[2] 有效区 {len(keep)}（100 - RH_SomMot_2）', flush=True)

modes = {
    'salience': sal_df.set_index('region')['salience'].reindex(parcels).values.astype(float),
    'discon_mean': ctrl['discon_mean'].reindex(parcels).values.astype(float),
}

# ---- 3. 逐细胞 spin ----
rows = []
for ct in ['Micro', 'Oligodendrocytes', 'Astro', 'Excitatory_Neurons', 'Inhibitory_Neurons', 'Endo', 'OPC']:
    gs = [g for g in PANEL_D.get(ct, []) if g in ez.columns]
    if not gs:
        continue
    cmap_100 = np.array([np.median(ez.loc[parcels[i], gs]) for i in range(100)])
    for mode_name, mode_100 in modes.items():
        cm = cmap_100[keep]; mv = mode_100[keep]
        r_obs = np.corrcoef(cm, mv)[0, 1]
        null = np.array([np.corrcoef(cm, mode_100[perm_id[:, k]][keep])[0, 1] for k in range(N_ROT)])
        p = (1 + np.sum(np.abs(null) >= abs(r_obs))) / (1 + N_ROT)
        rows.append({'feature': f'Cell_{ct}', 'mode': mode_name, 'n_genes': len(gs),
                     'r': round(r_obs, 4), 'P_SPIN': p})
        print(f'  Cell_{ct:<22} {mode_name:<11} n={len(gs):>4} r={r_obs:+.4f} P_SPIN={p:.4f}', flush=True)

# ---- 4. 保存 ----
df = pd.DataFrame(rows)
df.to_csv(OUT, index=False, encoding='utf-8-sig')
print(f'\n已保存: {OUT}')
