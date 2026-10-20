# -*- coding: utf-8 -*-
"""
Step1 复现：从 LQT 全脑断连矩阵裁剪皮层 100 区 + 一致性验证
=====================================================================
背景：新口径主分析使用皮层 100 区断连矩阵（PSF_SDC_parcel_discon_N124_cortex100.csv）。
      该矩阵【不是重新跑 LQT】得到，而是从 LQT 全脑版
      （PSF_SDC_parcel_discon_N124.csv, 124×135 区）直接裁剪 LH_/RH_ 皮层 100 列
      （与全脑版皮层列 max|diff|=0，已单独验证）。
本脚本：①复现裁剪；②验证裁剪结果与已发布 cortex100 完全一致。
LQT 全脑版计算本身由 psf_sdc_lqt_n124_v2.R 完成（本目录内，Pan 2022 同款工具
      LQT/Griffis 2021 + HCP-842 规范连接组；详见 PSF_SDC_LQT操作清单.md）。
输入: F:/dcmtonii/PSF_SDC_parcel_discon_N124.csv（LQT 全脑版, 124×136）
输出: F:/dcmtonii/PSF_SDC_parcel_discon_N124_cortex100.csv（验证后, 124×101）
"""
import pandas as pd
import numpy as np

SRC_FULL = 'F:/dcmtonii/PSF_SDC_parcel_discon_N124.csv'          # LQT 全脑版（R 脚本产出）
SRC_C100 = 'F:/dcmtonii/PSF_SDC_parcel_discon_N124_cortex100.csv' # 已发布 cortex100（裁剪源）
OUT = 'F:/dcmtonii/PSF_SDC_parcel_discon_N124_cortex100.csv'

full = pd.read_csv(SRC_FULL)
c100_pub = pd.read_csv(SRC_C100)

# 1. 裁剪：皮层列 = LH_/RH_ 前缀（Schaefer100 命名）
cort_cols = [c for c in full.columns if c.startswith('LH_') or c.startswith('RH_')]
print(f'[1] 全脑版: {full.shape} | 皮层列(LH_/RH_): {len(cort_cols)}')

id_col = full.columns[0]
c100 = full[[id_col] + cort_cols].copy()
print(f'[2] 裁剪后: {c100.shape}')

# 2. 验证：与已发布 cortex100 逐值一致（max|diff|=0）
common = [c for c in cort_cols if c in c100_pub.columns]
diff = (c100[common].values - c100_pub[common].values)
max_diff = np.nanmax(np.abs(diff)) if diff.size else float('nan')
print(f'[3] 与已发布版比较（{len(common)} 列）: max|diff| = {max_diff}')
assert max_diff < 1e-12, f'裁剪与已发布版不一致 (max|diff|={max_diff})'
print('[4] ✓ 一致（浮点容差内）— 裁剪逻辑复现成功')

# 3. 写出（若路径与已发布相同则覆盖为同一文件；如只想验证不改源可改 OUT）
c100.to_csv(OUT, index=False)
print(f'[5] 已写出: {OUT}')
