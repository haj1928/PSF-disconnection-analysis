# -*- coding: utf-8 -*-
"""
Step 2 CI 补充（2026-09-02）：corr(Lx,Ly) 的 bootstrap 95% 置信区间
=====================================================================
背书：行为 PLS 的 bootstrap 推断是标准做法（McIntosh et al. 1996 NeuroImage;
      Krishnan et al. 2011 NeuroImage 56:455 教程：bootstrap 重抽样被试估计
      统计量稳健性）；主研究（结构偏离）同法已报 boot95%CI=[+0.049,+0.438]
      （`_total_primary_N124.log`）——断连版与其保持内部一致。

做法（与主脚本 pls_pyls_discon_cov_residY_cortex100.py 同一 pipeline）：
  1. bootstrap 重抽样 124 患者（X/Y/协变量 同步行重抽, seed=0, 1000 次）
  2. 每次：Y 对 [1,年龄,性别] 残差化（bootstrap 样本内重新拟合）→
     pyls.behavioral_pls(n_perm=0, n_boot=0) 分解 → corr(Lx,Ly)
  3. 95% CI = bootstrap 分布 2.5/97.5 百分位
输出: 报告追加至 PSF_pls_discon_residY_age_sex_cortex100_report_newmasks_20260915.txt
      PSF_pls_behavior/PSF_pls_discon_bootCI_cortex100_newmasks_20260915.csv
"""
import numpy as np
import pandas as pd

SRC = 'F:/dcmtonii'
BASE = SRC + '/run62_crosssectional'
OUT = BASE + '/PSF_pls_behavior'
N_BOOT = 1000
SEED = 0

mt = pd.read_csv(BASE + '/PSF_master_table_baseline_acute.csv')
mt['sub'] = [f'sub_{x}' if not str(x).startswith('sub_') else str(x) for x in mt['subject_id']]
pd_ = pd.read_csv(SRC + '/PSF_SDC_parcel_discon_N124_cortex100_newmasks_20260915.csv')
ids = [str(x) for x in pd_['ID']]
X = pd_.drop(columns=['ID']).values.astype(float)
Y = pd.read_csv(BASE + '/FSS9_baseline_124.csv', index_col=0).loc[ids].values.astype(float)
mt_i = mt.set_index('sub').loc[ids]
C = np.column_stack([mt_i['年龄'].values.astype(float),
                     mt_i['性别_男1女0'].values.astype(float)])
n = len(ids)
print(f'X {X.shape} | Y {Y.shape} | C {C.shape} | N={n}', flush=True)

def residualize(M, C):
    Cw = np.column_stack([np.ones(len(C)), C])
    beta, _, _, _ = np.linalg.lstsq(Cw, M, rcond=None)
    return M - Cw @ beta

def pls_corr_np(Xb, Yb):
    """行为 PLS 核心 corr(Lx,Ly)（numpy SVD 版, 与 pyls xz=True 一致, 复现 0.3156≈0.3148）"""
    Xp = (Xb - Xb.mean(0)) / (Xb.std(0, ddof=1) + 1e-12)
    Yp = (Yb - Yb.mean(0)) / (Yb.std(0, ddof=1) + 1e-12)
    U, S, Vt = np.linalg.svd(Xp.T @ Yp, full_matrices=False)
    Lx = Xp @ U[:, 0]
    Ly = Yp @ Vt[0]
    return float(np.corrcoef(Lx, Ly)[0, 1])

def corr_once(Xb, Yb, Cb):
    Yr = residualize(Yb, Cb)
    return pls_corr_np(Xb, Yr)

# 观测（全样本, 与主分析一致）
obs = corr_once(X, Y, C)
print(f'观测 corr(Lx,Ly) = {obs:+.4f}（主分析 pyls = +0.3156）', flush=True)

rng = np.random.default_rng(SEED)
boot = np.empty(N_BOOT)
for b in range(N_BOOT):
    idx = rng.integers(0, n, n)
    boot[b] = corr_once(X[idx], Y[idx], C[idx])
    if (b + 1) % 200 == 0:
        print(f'  boot {b+1}/{N_BOOT}', flush=True)

lo, hi = np.percentile(boot, 2.5), np.percentile(boot, 97.5)
print(f'\ncorr(Lx,Ly) boot 95% CI = [{lo:+.4f}, {hi:+.4f}]（观测 {obs:+.4f}）', flush=True)

pd.DataFrame({'n_boot': [N_BOOT], 'seed': [SEED], 'obs_corr': [round(obs, 4)],
              'ci_lo': [round(lo, 4)], 'ci_hi': [round(hi, 4)],
              'pct_below_zero': [round(float((boot < 0).mean()), 4)]}
             ).to_csv(f'{OUT}/PSF_pls_discon_bootCI_cortex100_newmasks_20260915.csv', index=False, encoding='utf-8-sig')

with open(f'{OUT}/PSF_pls_discon_residY_age_sex_cortex100_report_newmasks_20260915.txt', 'a', encoding='utf-8') as f:
    f.write(f'\n[bootstrap 95% CI 补充 2026-09-02, N={N_BOOT}, seed={SEED}]\n')
    f.write(f'  corr(Lx,Ly) 观测 {obs:+.4f} | 95% CI = [{lo:+.4f}, {hi:+.4f}]\n')
    f.write(f'  bootstrap 分布中 <0 的比例: {float((boot < 0).mean()):.4f}\n')
print('已追加至 report + 保存 bootCI CSV')
