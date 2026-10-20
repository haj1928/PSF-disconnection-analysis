# -*- coding: utf-8 -*-
"""
蛋白稳态（核糖体/蛋白酶体）通路表达图 spin 补跑（2026-09-08 S2c 审计补全）
=====================================================================
背景：3.5 空间层验证原覆盖能量代谢（OXPHOS/糖酵解）与免疫（MHC I），
      缺蛋白稳态（第一段主发现两类核心之一）。本脚本补 KEGG Ribosome/Proteasome。
方法：与 run_sal_specificity_spin_v2.py 的 Path_* 行完全同法：
      通路基因集（KEGG_2021_Human.gmt）→ median regional expression score 表达图
      （逐基因跨区 Z → 每皮层区取基因集中位数，99 有效区）
      → × {salience, discon_mean} Pearson 相关 → Schaefer100 fsa5 spin 1000× seed=0 双侧
输入: PSF_ahba_100_expression.csv / salience CSV / PSF_discon_control_vectors_newmasks_20260915.csv / KEGG GMT
输出: PSF_ahba_100/permutation_newmasks_20260915/PSF_sal_specificity_spin_proteostasis_newmasks_20260915.csv
结果(2026-09-08): Ribosome salience r=+0.3405 P_SPIN=0.0020 (burden −0.1136 ns)
                   Proteasome salience r=+0.3989 P_SPIN=0.0010 (burden −0.0561 ns)
"""
import pandas as pd, numpy as np, nibabel as nb, sys, os
BASE = 'F:/dcmtonii/run62_crosssectional/'
PET = BASE + 'PSF_repro_PET'
EXPR = BASE + 'PSF_ahba_100/PSF_ahba_100_expression.csv'
SAL = BASE + 'PSF_pls_behavior/PSF_pls_discon_salience_residY_age_sex_cortex100_newmasks_20260915.csv'
CTRL = BASE + 'PSF_ahba_100/PSF_discon_control_vectors_newmasks_20260915.csv'
GMT = BASE + 'PSF_references/gene_sets/KEGG_2021_Human.gmt'
OUT = BASE + 'PSF_ahba_100/permutation_newmasks_20260915/PSF_sal_specificity_spin_proteostasis_newmasks_20260915.csv'
N_ROT, SEED = 1000, 0

expr = pd.read_csv(EXPR, index_col=0)
sal_df = pd.read_csv(SAL)
ctrl = pd.read_csv(CTRL).set_index('region')
parcels = sal_df['region'].tolist()

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

ez = (expr - expr.mean(0)) / (expr.std(0, ddof=1) + 1e-12)
keep = [i for i in range(100) if not expr.loc[parcels[i]].isna().all()]
modes = {'salience': sal_df.set_index('region')['salience'].reindex(parcels).values.astype(float),
         'discon_mean': ctrl['discon_mean'].reindex(parcels).values.astype(float)}

def load_gmt(fn):
    out = {}
    for line in open(fn, encoding='utf-8'):
        parts = line.rstrip('\n').split('\t')
        if len(parts) >= 3:
            out[parts[0]] = parts[2:]
    return out
kegg = load_gmt(GMT)

rows = []
for name in ['Ribosome', 'Proteasome']:
    gs = [g for g in kegg.get(name, []) if g in ez.columns]
    cmap = np.array([np.median(ez.loc[parcels[i], gs]) for i in range(100)])
    for mname, mode in modes.items():
        cm, mv = cmap[keep], mode[keep]
        r_obs = np.corrcoef(cm, mv)[0, 1]
        null = np.array([np.corrcoef(cm, mode[perm_id[:, k]][keep])[0, 1] for k in range(N_ROT)])
        p = (1 + np.sum(np.abs(null) >= abs(r_obs))) / (1 + N_ROT)
        # 2026-10-09：新增同向单侧（论文定稿口径；报告 44 §0）
        p_one = ((1 + np.sum(null >= r_obs)) / (1 + N_ROT)) if r_obs >= 0 else ((1 + np.sum(null <= r_obs)) / (1 + N_ROT))
        rows.append({'feature': f'Path_KEGG_{name}', 'mode': mname, 'n_genes': len(gs),
                     'r': round(r_obs, 4), 'P_SPIN': p, 'P_SPIN_one': p_one})
        print(f'  Path_KEGG_{name:<11} {mname:<11} n={len(gs):>4} r={r_obs:+.4f} P_SPIN={p:.4f} P_one={p_one:.4f}')
pd.DataFrame(rows).to_csv(OUT, index=False, encoding='utf-8-sig')
print('saved:', OUT)
