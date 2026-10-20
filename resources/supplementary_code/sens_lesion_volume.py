# -*- coding: utf-8 -*-
"""
病灶体积敏感性：断连×疲劳 PLS 是否超越病灶大小（皮层100 主口径同参）
=====================================================================
主分析核对：仅残差 Y(FSS9) 对 [年龄,性别]            → 期望 perm_p≈0.0232
S1 (仅Y侧加体积): Y 对 [年龄,性别,病灶体积] 残差, X 不动
S2 (双侧控体积):  X 每列对 [年龄,性别,病灶体积] 残差, Y 同 S1  → 探索性
                （双残差+简单置换的膨胀注意: Winkler 2020, 显著性解读从保守）
全部: n_perm=5000, n_boot=5000, seed=0（与主分析同参）
病灶体积: F:/dcmtonii/PSF_lesion_volume_N124.csv (subject, volume_mL)
"""
import numpy as np
import pandas as pd
from scipy.stats import spearmanr

if not hasattr(np, 'row_stack'):
    np.row_stack = np.vstack
import pyls

SRC = 'F:/dcmtonii'
BASE = SRC + '/run62_crosssectional'
OUT = BASE + '/PSF_pls_behavior'

def sid(x):
    return f'sub_{x}' if not str(x).startswith('sub_') else str(x)

mt = pd.read_csv(BASE + '/PSF_master_table_baseline_acute.csv')
mt['sub'] = mt['subject_id'].map(sid)
pd_ = pd.read_csv(SRC + '/PSF_SDC_parcel_discon_N124_cortex100_newmasks_20260915.csv')
pd_['sub'] = pd_['ID'].map(sid)
X = pd_.drop(columns=['ID', 'sub']).values.astype(float)
ids = [str(x) for x in pd_['ID']]
pnames = [c for c in pd_.columns if c != 'ID']
print(f'断连: {X.shape} | 区数 {len(pnames)}', flush=True)

Y9 = pd.read_csv(BASE + '/FSS9_baseline_124.csv', index_col=0)
Y = Y9.loc[ids].values.astype(float)

lv = pd.read_csv(SRC + '/PSF_lesion_volume_N124.csv')
lv['sub'] = lv['subject'].map(sid)
vol = lv.set_index('sub')['volume_mL'].loc[ids].values.astype(float)
print(f'病灶体积: n={len(vol)} 缺失 {np.isnan(vol).sum()} | '
      f'范围 [{np.nanmin(vol):.2f}, {np.nanmax(vol):.2f}] mL | '
      f'中位 {np.nanmedian(vol):.2f} mL', flush=True)

mt_i = mt.set_index('sub').loc[ids]
age = mt_i['年龄'].values.astype(float)
sex = mt_i['性别_男1女0'].values.astype(float)

def residualize(M, C):
    M = np.asarray(M, float)
    Cw = np.column_stack([np.ones(len(C)), C])
    beta, *_ = np.linalg.lstsq(Cw, M, rcond=None)
    return M - Cw @ beta

# ---- 描述：病灶体积与其他变量（Spearman）----
FSSsum = Y.sum(axis=1)
desc = []
for nm, v in [('FSS9总分', FSSsum), ('断连总量', X.mean(axis=1)),
              ('年龄', age), ('HAMD_total', mt_i['HAMD_total'].values.astype(float))]:
    m = ~(np.isnan(vol) | np.isnan(v))
    rho, p = spearmanr(vol[m], v[m])
    desc.append(f'  rho(体积, {nm}) = {rho:+.3f}  p={p:.4f}  (n={m.sum()})')

def run_pls(Xin, Yin, tag):
    res = pyls.behavioral_pls(Xin, Yin, n_perm=5000, n_boot=5000, seed=0,
                              verbose=False, test_split=None)
    v1 = res.x_weights[:, 0]
    ratio = res.bootres['x_weights_normed'][:, 0]
    sig = int(np.sum(np.abs(ratio) > 2))
    perm_p = float(res.permres['pvals'][0])
    r = float(np.corrcoef(res.x_scores[:, 0], res.y_scores[:, 0])[0, 1])
    ve = float(res.singvals[0] ** 2 / np.sum(res.singvals ** 2) * 100)
    print(f'[{tag}] perm_p={perm_p:.4f} | corr(Lx,Ly)={r:+.4f} | varexp={ve:.2f}% | BSR>2 {sig}/{len(pnames)}', flush=True)
    return dict(tag=tag, perm_p=round(perm_p, 4), corr_Lx_Ly=round(r, 4),
                varexp_pct=round(ve, 2), n_sig=sig)

# ---- 主分析核对（仅 Y 残差 [年龄,性别]）----
Yr0 = residualize(Y, np.column_stack([age, sex]))
r_main = run_pls(X, Yr0, '主分析核对(年龄+性别)')

# ---- S1：仅 Y 侧加病灶体积 ----
C1 = np.column_stack([age, sex, vol])
Yr1 = residualize(Y, C1)
r_s1 = run_pls(X, Yr1, 'S1_仅Y残差[年龄,性别,体积]')

# ---- S2：双侧控 [年龄,性别,体积]（探索）----
Xr2 = residualize(X, C1)
r_s2 = run_pls(Xr2, Yr1, 'S2_双侧残差[年龄,性别,体积]')

lines = []
lines.append('=' * 78)
lines.append('病灶体积敏感性：断连×疲劳 PLS（皮层100, 仅残差Y 口径, N=124）')
lines.append('病灶体积文件: F:/dcmtonii/PSF_lesion_volume_N124.csv')
lines.append('-' * 78)
lines.append('体积描述: 范围 [%.2f, %.2f] mL | 中位 %.2f mL | 缺失 %d' %
             (np.nanmin(vol), np.nanmax(vol), np.nanmedian(vol), np.isnan(vol).sum()))
lines.append('\n[描述] 病灶体积 × 其他变量 (Spearman):')
lines += desc
lines.append('\n[PLS 同参对照] n_perm=5000 n_boot=5000 seed=0:')
for r in [r_main, r_s1, r_s2]:
    lines.append('  %-28s perm_p=%.4f  corr=%+.4f  varexp=%.2f%%  BSR>2=%d/100' %
                 (r['tag'], r['perm_p'], r['corr_Lx_Ly'], r['varexp_pct'], r['n_sig']))
lines.append('''
[解读]
- S1 vs 主分析: Y 侧加体积后仍显著 → 断连模式关联超越病灶大小（仅控人口学层面的检验）
- S2 (探索): 双侧控体积后是否存活 → "断连模式在体积/年龄/性别之外的独立贡献"。
  注意双残差+简单置换可能膨胀错误率(Winkler 2020 NeuroImage 220:117065), S2 不显著不能当证伪,
  显著需复算确认; S1 与主口径方法链一致, 无此争议。
- 若体积与 FSS9 相关弱 → 体积本就不构成强混杂, 控与不控差异小（看[描述]行）。
''')
lines.append('=' * 78)
txt = '\n'.join(lines)
print('\n' + txt, flush=True)
with open(OUT + '/_sens_lesion_volume_report_newmasks_20260915.txt', 'w', encoding='utf-8') as f:
    f.write(txt)

df = pd.DataFrame([r_main, r_s1, r_s2])
df.to_csv(OUT + '/_sens_lesion_volume_summary_newmasks_20260915.csv', index=False, encoding='utf-8-sig')
print('\n[已存] PSF_pls_behavior/_sens_lesion_volume_report_newmasks_20260915.txt / _summary.csv')
