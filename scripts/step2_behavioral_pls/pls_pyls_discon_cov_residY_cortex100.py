# -*- coding: utf-8 -*-
"""
Step 2b'（皮层 100 区主分析版）：仅残差化 Y(FSS9) 对 [年龄,性别]，X 不动
=====================================================================
背景（2026-09-01 用户拍板）：
  - 主分析协变量 = 年龄 + 性别（与主研究 §2.5 主检验同口径）。
  - NIHSS/mRS/HAMD 从主分析移除（NIHSS/mRS 是断连下游中介，控制属
    overadjustment；HAMD 与疲劳高度共线，控制会吸收真实方差）。
  - 不做双残差：Winkler et al. 2020 (NeuroImage 220:117065) 证明残差化后
    简单置换检验会膨胀错误率（残差化引入观测间依赖、违反可交换性）；
    Petrican et al. 2022 (NeuroImage 255:119209) 明确仅残差化非成像
    （行为）侧并引 Winkler 2020。本脚本遵循该先例：只残差 Y(FSS9)。
  - 数学上单侧残差对线性混杂的排除与双侧等效：Yr ⊥ C 时，
    <P_C X, Yr> = <X, P_C Yr> = 0，另一侧协变量线性成分的共变自动归零。

做法：
  Y_r = 每列对 [1, 年龄, 性别] 线性回归取残差（含截距）
  X    = 皮层 100 区断连矩阵（124×100，不残差）
  PLS  = pyls.behavioral_pls(X, Y_r, n_perm=5000, n_boot=5000, seed=0)

输出：
  PSF_pls_discon_residY_age_sex_cortex100.csv        （汇总：perm_p/corr/varexp/n_sig）
  PSF_pls_discon_salience_residY_age_sex_cortex100.csv（salience + BSR，供下游）
  PSF_pls_discon_residY_age_sex_cortex100_report.txt  （完整报告）
=====================================================================
"""
import numpy as np
import pandas as pd

# ---- pyls 兼容补丁（numpy 2.5 移除 row_stack）----
if not hasattr(np, 'row_stack'):
    np.row_stack = np.vstack

import pyls

SRC = 'F:/dcmtonii'
BASE = SRC + '/run62_crosssectional'
OUT = BASE + '/PSF_pls_behavior'

# ---------------- 1. 读数据 ----------------
mt = pd.read_csv(BASE + '/PSF_master_table_baseline_acute.csv')
subs = [f'sub_{x}' if not str(x).startswith('sub_') else str(x) for x in mt['subject_id']]
mt['sub'] = subs

pd_ = pd.read_csv(SRC + '/PSF_SDC_parcel_discon_N124_cortex100.csv')   # 皮层 100 区
ids = [str(x) for x in pd_['ID']]
Xp = pd_.drop(columns=['ID']).values.astype(float)                     # 124×100
pnames = [c for c in pd_.columns if c != 'ID']
n = len(ids)
print(f'断连(皮层100): {Xp.shape} | 区数 {len(pnames)} | N={n}', flush=True)

Y9 = pd.read_csv(BASE + '/FSS9_baseline_124.csv', index_col=0)
Y = Y9.loc[ids].values.astype(float)                                   # 124×9
print(f'FSS9: {Y.shape}', flush=True)

# ---------------- 2. 协变量：仅年龄 + 性别 ----------------
mt_i = mt.set_index('sub').loc[ids]
for c in ['年龄', '性别_男1女0']:
    n_miss = int(mt_i[c].isna().sum())
    if n_miss:
        print(f'  [警告] 协变量 {c} 缺失 {n_miss} 例', flush=True)
age = mt_i['年龄'].values.astype(float)
sex = mt_i['性别_男1女0'].values.astype(float)   # 男=1/女=0
C = np.column_stack([age, sex])
print(f'协变量 C=[年龄,性别] {C.shape} | 缺失: 年龄 {int(mt_i["年龄"].isna().sum())} 性别 {int(mt_i["性别_男1女0"].isna().sum())}', flush=True)


# ---------------- 3. 仅残差化 Y（FSS9），X 不动 ----------------
def residualize(M, C):
    """M 每列对 [1, C] 线性回归取残差（含截距）"""
    M = np.asarray(M, float)
    Cw = np.column_stack([np.ones(len(C)), C])
    beta, _, _, _ = np.linalg.lstsq(Cw, M, rcond=None)
    return M - Cw @ beta

Yr = residualize(Y, C)   # 仅 Y 残差化
print(f'Y(FSS9) 对 [年龄,性别] 残差化完成 | X 未动', flush=True)

# ---------------- 4. 行为 PLS ----------------
print('跑 pyls.behavioral_pls(seed=0, n_perm=5000, n_boot=5000) ...', flush=True)
res = pyls.behavioral_pls(Xp, Yr, n_perm=5000, n_boot=5000, seed=0,
                          verbose=False, test_split=None)
v1 = res.x_weights[:, 0]                                   # 脑 salience（LV1）
ratio = res.bootres['x_weights_normed'][:, 0]              # bootstrap ratio
sig_idx = np.where(np.abs(ratio) > 2)[0]                   # 显著区（|BSR|>2 ≈ 95%CI）
perm_p = float(res.permres['pvals'][0])
r_ly = float(np.corrcoef(res.x_scores[:, 0], res.y_scores[:, 0])[0, 1])
varexp = float(res.singvals[0] ** 2 / np.sum(res.singvals ** 2) * 100)
print(f'perm_p={perm_p:.4f} | corr(Lx,Ly)={r_ly:+.4f} | varexp={varexp:.2f}% | 显著区 {len(sig_idx)}/{len(pnames)}', flush=True)

# ---------------- 5. 输出 ----------------
df_sum = pd.DataFrame([dict(step='Step2b_主分析_仅残差Y',
                            n=n, n_regions=len(pnames), perm_p=round(perm_p, 4),
                            corr_Lx_Ly=round(r_ly, 4), varexp_pct=round(varexp, 2),
                            n_sig=int(len(sig_idx)),
                            covariate='年龄+性别', resid_side='Y(FSS9) only',
                            X_resid=False, Y_resid=True,
                            NIHSS_mRS_HAMD='excluded')])
df_sum.to_csv(OUT + '/PSF_pls_discon_residY_age_sex_cortex100.csv',
              index=False, encoding='utf-8-sig')

df_sal = pd.DataFrame(dict(region=pnames, salience=np.round(v1, 6),
                           BSR=np.round(ratio, 4),
                           sig=ratio > 2))
df_sal.to_csv(OUT + '/PSF_pls_discon_salience_residY_age_sex_cortex100.csv',
              index=False, encoding='utf-8-sig')

rep = [
    '=' * 74,
    'Step 2b\' 主分析：仅残差化 Y(FSS9) 对 [年龄,性别]，X 不动（皮层 100 区）',
    f'数据: N={n} 急性期 | X=断连 124×100 皮层(未残差) | Y=FSS9 124×9(残差化)',
    f'协变量: 年龄 + 性别 | NIHSS/mRS/HAMD 已从主分析移除',
    '方法: Petrican et al. 2022 (NeuroImage 255:119209) 先例——仅残差化行为侧,',
    '      引 Winkler et al. 2020 (NeuroImage 220:117065) 规避双残差+置换偏差',
    '-' * 74,
    f'  perm_p = {perm_p:.4f}',
    f'  corr(Lx,Ly) = {r_ly:+.4f}',
    f'  varexp(LV1) = {varexp:.2f}%',
    f'  显著区(|BSR|>2) = {len(sig_idx)}/{len(pnames)}',
    '-' * 74,
    '判读: 主分析仅控年龄+性别；断连-疲劳共变经人口学控制后的稳健性由此值体现。',
    '      注意: 未控 NIHSS/mRS(断连下游中介, 控之属 overadjustment, Schisterman 2009)',
    '      与 HAMD(与疲劳共线 rho≈0.47, 控之吸收真实方差)。',
    '=' * 74,
]
with open(OUT + '/PSF_pls_discon_residY_age_sex_cortex100_report.txt', 'w', encoding='utf-8') as f:
    f.write('\n'.join(rep))
print('\n'.join(rep), flush=True)
print(f'\n[输出]', flush=True)
print(f'  {OUT}/PSF_pls_discon_residY_age_sex_cortex100.csv', flush=True)
print(f'  {OUT}/PSF_pls_discon_salience_residY_age_sex_cortex100.csv', flush=True)
print(f'  {OUT}/PSF_pls_discon_residY_age_sex_cortex100_report.txt', flush=True)
