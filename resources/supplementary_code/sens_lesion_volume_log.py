# -*- coding: utf-8 -*-
"""
病灶体积敏感性——稳健性核查
疑点：S1 仅Y侧加体积后 perm_p 0.0232→0.1954，但 rho(体积,FSS9)=0.024(ns)——
      若体积与疲劳无关，残差化几乎不该动 Y。（右偏分布）
      （体积右偏 0.10~81.68 mL, 中位 2.51）。
核查：
  1) 体积 × FSS9 逐条目 Spearman（确认条目级是否也无关）
  2) log10(体积) 后重跑 S1/S2（同参 n_perm=5000 n_boot=5000 seed=0）
  3) 【已于 2026-09-22 删除】最大病灶个体的杠杆点检查
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

Y9 = pd.read_csv(BASE + '/FSS9_baseline_124.csv', index_col=0)
Y = Y9.loc[ids].values.astype(float)

lv = pd.read_csv(SRC + '/PSF_lesion_volume_N124.csv')
vol = lv.set_index('subject')['volume_mL'].loc[ids].values.astype(float)
logvol = np.log10(vol)

mt_i = mt.set_index('sub').loc[ids]
age = mt_i['年龄'].values.astype(float)
sex = mt_i['性别_男1女0'].values.astype(float)

def residualize(M, C):
    M = np.asarray(M, float)
    Cw = np.column_stack([np.ones(len(C)), C])
    beta, *_ = np.linalg.lstsq(Cw, M, rcond=None)
    return M - Cw @ beta

def run_pls(Xin, Yin, tag, perm=5000):
    res = pyls.behavioral_pls(Xin, Yin, n_perm=perm, n_boot=5000, seed=0,
                              verbose=False, test_split=None)
    ratio = res.bootres['x_weights_normed'][:, 0]
    sig = int(np.sum(np.abs(ratio) > 2))
    perm_p = float(res.permres['pvals'][0])
    r = float(np.corrcoef(res.x_scores[:, 0], res.y_scores[:, 0])[0, 1])
    print(f'[{tag}] perm_p={perm_p:.4f} | corr={r:+.4f} | BSR>2 {sig}/100', flush=True)
    return dict(tag=tag, perm_p=round(perm_p, 4), corr=round(r, 4), n_sig=sig)

out = []
out.append('=' * 78)
out.append('病灶体积稳健性核查（右偏→log 变换）')

# 1) 逐条目
out.append('\n[1] 体积 × FSS9 逐条目 Spearman:')
for i in range(9):
    rho, p = spearmanr(vol, Y[:, i])
    out.append('  rho(体积, FSS_%d) = %+.3f  p=%.4f' % (i + 1, rho, p))
rho_t, p_t = spearmanr(vol, Y.sum(1))
out.append('  rho(体积, FSS总分) = %+.3f  p=%.4f' % (rho_t, p_t))
rho_l, p_l = spearmanr(logvol, Y.sum(1))
out.append('  rho(log10体积, FSS总分) = %+.3f  p=%.4f' % (rho_l, p_l))

# 2) log 版 S1/S2
C = np.column_stack([age, sex, logvol])
Yr0 = residualize(Y, np.column_stack([age, sex]))
r_main = run_pls(X, Yr0, '主分析核对(年龄+性别)')
Yr1 = residualize(Y, C)
r_s1 = run_pls(X, Yr1, 'S1_log: 仅Y残差[年龄,性别,log体积]')
Xr2 = residualize(X, C)
r_s2 = run_pls(Xr2, Yr1, 'S2_log: 双侧残差[年龄,性别,log体积]')

# 3) 杠杆点检查已于 2026-09-22 随稿件删除（该检验不再报告）

out.append('\n[汇总]')
for r in [r_main, r_s1, r_s2]:
    out.append('  %-34s perm_p=%.4f  corr=%+.4f  BSR>2=%d/100' %
               (r['tag'], r['perm_p'], r['corr'], r['n_sig']))
out.append('''
[解读]
- 若 log 版 S1 恢复显著 → 原 S1 不显著=右偏杠杆伪影, 体积非真混杂（其与 FSS 本就无关）
- S2 log 版: 双侧控体积后是否存活 → 断连模式超越体积的独立贡献（探索性, Winkler 2020 注意）
''')
out.append('=' * 78)
txt = '\n'.join(out)
print('\n' + txt, flush=True)
with open(OUT + '/_sens_lesion_volume_log_report_newmasks_20260915.txt', 'w', encoding='utf-8') as f:
    f.write(txt)
print('[已存] PSF_pls_behavior/_sens_lesion_volume_log_report_newmasks_20260915.txt')
