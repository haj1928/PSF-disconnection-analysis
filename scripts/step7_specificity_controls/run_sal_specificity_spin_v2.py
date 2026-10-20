# -*- coding: utf-8 -*-
"""
salience 特异性检验 v2（2026-09-02）——与已验证单测同构
=====================================================================
聚焦目标：细胞/通路表达图与 salience 的共变是 salience 特异的，
  而非断连总量共有。
方法（Seidlitz 2020 median + Alexander-Bloch spin）：
  gene set map = 基因集内 median Z（100 区）
  P_SPIN = corr(map, 模式) vs 模式 1000× spin 旋转 null（双侧）
输出: PSF_ahba_100/permutation/PSF_sal_specificity_spin_result_newmasks_20260915.csv + report
"""
import numpy as np
import pandas as pd
import sys
import nibabel as nb

BASE = 'F:/dcmtonii/run62_crosssectional'
PET = BASE + '/PSF_repro_PET'
OUT = BASE + '/PSF_ahba_100/permutation_newmasks_20260915'
GMT = BASE + '/PSF_references/gene_sets'
N_SPIN = 1000

expr = pd.read_csv(BASE + '/PSF_ahba_100/PSF_ahba_100_expression.csv', index_col=0)
sal_df = pd.read_csv(BASE + '/PSF_pls_behavior/PSF_pls_discon_salience_residY_age_sex_cortex100_newmasks_20260915.csv')
ctrl = pd.read_csv(BASE + '/PSF_ahba_100/PSF_discon_control_vectors_newmasks_20260915.csv').set_index('region')
parcels = sal_df['region'].tolist()
assert list(expr.index) == parcels

ez = (expr - expr.mean(0)) / (expr.std(0, ddof=1) + 1e-12)     # 100 × G
mode_sal = sal_df.set_index('region')['salience'].reindex(parcels).values
mode_mean = ctrl['discon_mean'].reindex(parcels).values
assert not np.isnan(mode_sal).any() and not np.isnan(mode_mean).any()
valid_idx = [i for i, p in enumerate(parcels) if not expr.loc[p].isna().all()]
print(f'[0] 有效区 {len(valid_idx)} | ez {ez.shape}', flush=True)

# ---- spin perm_id ----
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
perm_id = rotate_parcellation(np.array(cl), np.array(cr), nrot=N_SPIN, seed=0).astype(int)
print(f'[1] perm_id {perm_id.shape}', flush=True)

def spin_p(feat_map_100, mode_arr):
    cm = feat_map_100[valid_idx]
    mv = mode_arr[valid_idx]
    r_obs = float(np.corrcoef(cm, mv)[0, 1])
    null = np.array([np.corrcoef(cm, mode_arr[perm_id[:, k]][valid_idx])[0, 1]
                     for k in range(N_SPIN)])
    p = (1 + int((np.abs(null) >= abs(r_obs)).sum())) / (1 + N_SPIN)
    # 2026-10-09：新增同向单侧（论文定稿口径；报告 44 §0）
    p_one = ((1 + int((null >= r_obs).sum())) / (1 + N_SPIN)) if r_obs >= 0 else ((1 + int((null <= r_obs).sum())) / (1 + N_SPIN))
    return r_obs, p, p_one, float(np.median(null))

CELL_SETS = {
    'Micro': ['CX3CR1','CSF1R','AIF1','TMEM119','PTPRC','ITGAM','CD68','TYROBP','P2RY12','TREM2'],
    'Astro': ['GFAP','AQP4','SLC1A2','SLC1A3','ALDH1L1','S100B','GJA1','SOX9','ALDOC'],
    'Oligo': ['MBP','PLP1','MOG','MAG','CNP','MOBP','CLDN11','OPALIN','ST18'],
    'Exci': ['SLC17A7','SATB2','NEFL','NEFM','CAMK2A','RORB','CUX2','FEZF2','TBR1','NRGN'],
    'Inhi': ['GAD1','GAD2','SLC32A1','PVALB','SST','VIP','LAMP5','RELN','CALB1','CALB2'],
    'OPC': ['PDGFRA','CSPG4','OLIG1','OLIG2','SOX10','GPR17','LHFPL3'],
    'Endo': ['CLDN5','FLT1','PECAM1','VWF','ABCB1','SLC2A1','CDH5','TEK'],
}
def load_gmt(fn):
    sets = {}
    with open(fn, encoding='utf-8', errors='ignore') as f:
        for line in f:
            parts = line.rstrip('\n').split('\t')
            if len(parts) >= 3:
                sets[parts[0]] = parts[2:]
    return sets
go = load_gmt(f'{GMT}/GO_Biological_Process_2021.gmt')
kegg = load_gmt(f'{GMT}/KEGG_2021_Human.gmt')
PATH_SETS = {
    'KEGG_OXPHOS': kegg.get('Oxidative phosphorylation'),
    'KEGG_Glycolysis': kegg.get('Glycolysis / Gluconeogenesis'),
    'GO_glycolysis': go.get('glycolytic process (GO:0006096)'),
    'GO_aerobic_ETC': go.get('aerobic electron transport chain (GO:0019646)'),
    'GO_synapse_pruning': go.get('synapse pruning (GO:0098883)'),
    # ⚠ 2026-09-16 修正：原用 next(k.startswith('antigen processing and presentation
    # of exogenous peptide antigen')) 取 GMT 中**第一个**同前缀条目，实际命中
    # GO:0002478（103 基因，非 MHC class I），而 Table 3 该行声明的是
    # GO:0002479（TAP 依赖的 MHC class I 抗原提呈，73 基因）——两者不同，
    # 造成「同一行混用两个基因集」（NES/spin/基因集大小 来自 GO:0002479，
    # r 来自 GO:0002478）。此处改为**显式取 GO:0002479**，并把特征名写清。
    'GO_MHC_I_antigen_TAP': go.get(
        'antigen processing and presentation of exogenous peptide antigen '
        'via MHC class I, TAP-dependent (GO:0002479)'),
}

rows = []
for tag, genes in {**{'Cell_' + k: v for k, v in CELL_SETS.items()},
                   **{'Path_' + k: v for k, v in PATH_SETS.items()}}.items():
    if not genes:
        print(f'[跳过] {tag}: 空集'); continue
    gs = [g for g in genes if g in ez.columns]
    if not gs:
        print(f'[跳过] {tag}: 无基因命中'); continue
    feat_map = ez[gs].median(axis=1).reindex(parcels).values   # median 聚合（Seidlitz 用法）
    for mode_name, mode_arr in [('salience', mode_sal), ('discon_mean', mode_mean)]:
        r_obs, p, p_one, null_m = spin_p(feat_map, mode_arr)
        rows.append({'feature': tag, 'mode': mode_name, 'n_genes': len(gs),
                     'r': round(r_obs, 4), 'spin_null_median': round(null_m, 4),
                     'P_SPIN': p, 'P_SPIN_one': p_one})
        print(f'  {tag:<24} × {mode_name:<12} r={r_obs:+.3f} P_SPIN={p:.4f}', flush=True)

res = pd.DataFrame(rows)
res.to_csv(f'{OUT}/PSF_sal_specificity_spin_result_newmasks_20260915.csv', index=False, encoding='utf-8-sig')
# 报告还需 spin_null_median（原漏出，导致报告 null 列为 nan）
wide = res.pivot_table(index='feature', columns='mode',
                       values=['r', 'P_SPIN', 'P_SPIN_one', 'spin_null_median'])
# 2026-09-16 修正：原写作 f'{b}_{a}' 得到 'salience_r'，而下方循环读的是
# 'r_salience' → 报告表格全 NaN（CSV 本身正确）。改为 f'{a}_{b}'。
wide.columns = [f'{a}_{b}' for a, b in wide.columns]
wide = wide.reset_index()
rep = ['=' * 102,
       'salience 特异性检验（细胞+通路表达图 × {salience, 断连总量} spin, median 聚合）',
       '判读: salience P_SPIN_one<0.05 ★ + discon_mean P_SPIN_one≥0.05 → salience 特异（同向单侧）',
       '-' * 102,
       f"{'特征':<24}{'r_sal':>9}{'P_sal':>9}{'P1_sal':>9}{'null_sal':>10}{'r_mean':>9}{'P_mean':>9}{'P1_mean':>9}{'null_mean':>11}   判读"]
for _, r in wide.iterrows():
    # 2026-10-09：判读改用同向单侧（论文口径）
    ps_ = r.get('P_SPIN_one_salience', np.nan); pm_ = r.get('P_SPIN_one_discon_mean', np.nan)
    if np.isnan(ps_): interp = 'NA'; mark = ''
    elif ps_ < 0.05 and (np.isnan(pm_) or pm_ >= 0.05): interp = 'salience 特异'; mark = ' ★'
    elif ps_ < 0.05 and pm_ < 0.05: interp = '两模式均显著'; mark = ''
    else: interp = '不显著'; mark = ''
    rep.append(f"{r['feature']:<24}{r.get('r_salience',np.nan):>+9.3f}{r.get('P_SPIN_salience',np.nan):>9.4f}{ps_:>9.4f}{r.get('spin_null_median_salience',np.nan):>10.3f}{r.get('r_discon_mean',np.nan):>+9.3f}{r.get('P_SPIN_discon_mean',np.nan):>9.4f}{pm_:>9.4f}{r.get('spin_null_median_discon_mean',np.nan):>11.3f}   {interp}{mark}")
with open(f'{OUT}/PSF_sal_specificity_spin_report_newmasks_20260915.txt', 'w', encoding='utf-8') as f:
    f.write('\n'.join(rep))
print('\n' + '\n'.join(rep))
print(f'\n完成: {OUT}/PSF_sal_specificity_spin_result_newmasks_20260915.csv')
