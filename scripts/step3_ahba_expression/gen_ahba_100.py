# =============================================================================
# gen_ahba_100.py —— 用 abagen 生成 皮层 100 区 AHBA 基因表达矩阵
# =============================================================================
# 输入：
#   atlas      = LQT 图谱 100Parcels7Networks.nii.gz（MNI 182³，135 区）
#   atlas_info = 100Parcels7Networks_abagen_info.csv 的【皮层子集】（structure=='cortex'，id 1-100）
#   data_dir   = F:/dcmtonii/abagen-data（6 个供体 zip 已就位并 CRC 验证）
# 参数：全部 abagen 默认（Markello 2021 eLife 推荐：ibf_threshold=0.5,
#   probe_selection=diff_stability, region_agg=donors, norm=...）
# 输出：
#   PSF_ahba_100_expression.csv（100 皮层区 × 基因）
#   PSF_ahba_100_expression_report.txt（缺失区报告）
# =============================================================================
# 说明：只用皮层 100 区（跟随 Hansen 2022 / Seidlitz 2020 权威先例）；
#       abagen 只对 atlas_info 中列出的 region 返回数据，非皮层区自然不进入结果。
# =============================================================================

import pandas as pd
import numpy as np
import abagen

BASE = 'F:/dcmtonii/run62_crosssectional'
ATLAS = BASE + '/extdata/Schaefer_Yeo_Plus_Subcort/100Parcels7Networks.nii.gz'
INFO  = BASE + '/extdata/Schaefer_Yeo_Plus_Subcort/100Parcels7Networks_abagen_info.csv'
DATA  = 'F:/dcmtonii/abagen-data'
OUT   = BASE + '/PSF_ahba_100'

import os
os.makedirs(OUT, exist_ok=True)


def main():
    # Windows multiprocessing 铁律：主逻辑放在函数 + __main__ 保护内
    print(f'图谱: {ATLAS}', flush=True)
    print(f'数据: {DATA}', flush=True)

    # 1. atlas_info 必须覆盖图谱全部 labels（abagen 校验要求），但结果只取皮层
    info_full = pd.read_csv(INFO)
    print(f'atlas_info 全量: {len(info_full)} 行（abagen 校验要求覆盖全部 135 区）', flush=True)
    cortex_ids = info_full[info_full['structure'] == 'cortex']['id'].tolist()
    print(f'目标皮层 100 区 id: {len(cortex_ids)} 个（{cortex_ids[0]}-{cortex_ids[-1]}）', flush=True)

    # 2. abagen 映射（atlas_info 用全量，满足校验）
    print('开始 abagen.get_expression_data（6 供体，单进程）...', flush=True)
    expr, report = abagen.get_expression_data(
        ATLAS,
        atlas_info=info_full,
        data_dir=DATA,
        donors='all',
        return_report=True,
        verbose=1,
        n_proc=1,
    )

    # 3. 只保留皮层 100 区（abagen 返回 index = 区 id 整数，直接按 id 过滤）
    expr_cortex = expr.loc[[i for i in cortex_ids if i in expr.index]]
    print(f'皮层提取: {expr.shape[0]} 区 → {expr_cortex.shape[0]} 区（按 id 1-100）', flush=True)

    # ---- 保存表达矩阵（皮层 100 区，index 转为 label 便于与 salience 对齐）----
    id2label = info_full.set_index('id')['label'].to_dict()
    expr_cortex.index = [id2label.get(int(i), i) for i in expr_cortex.index]
    expr_cortex.index.name = 'region'
    expr_cortex.to_csv(OUT + '/PSF_ahba_100_expression.csv', encoding='utf-8-sig')
    print(f'表达矩阵已存: {expr_cortex.shape[0]} 区 × {expr_cortex.shape[1]} 基因（index=label）', flush=True)

    # ---- 缺失报告 ----
    nan_rows = expr_cortex[expr_cortex.isna().all(axis=1)]
    nan_partial = expr_cortex[expr_cortex.isna().any(axis=1) & ~expr_cortex.isna().all(axis=1)]
    with open(OUT + '/PSF_ahba_100_expression_report.txt', 'w', encoding='utf-8') as f:
        f.write(f'皮层 100 区表达矩阵: {expr_cortex.shape[0]} 区 × {expr_cortex.shape[1]} 基因\n')
        f.write(f'全 NaN 区（AHBA 无采样）: {len(nan_rows)} 个\n')
        f.write('  区: ' + ', '.join(nan_rows.index.astype(str)) + '\n\n')
        f.write(f'部分 NaN 区: {len(nan_partial)} 个\n')
        f.write('  区: ' + ', '.join(nan_partial.index.astype(str)) + '\n\n')
    print('缺失报告已存: PSF_ahba_100_expression_report.txt', flush=True)

    # ---- 与 135 区矩阵的皮层部分对比（验证一致性；135 版 index=id，映射后比）----
    try:
        expr135 = pd.read_csv(BASE + '/PSF_ahba_135/PSF_ahba_135_expression.csv', index_col=0)
        # 135 版 index 是 id → 映射到 label
        expr135_labels = expr135.copy()
        expr135_labels.index = [id2label.get(int(i), i) for i in expr135_labels.index]
        common = set(expr_cortex.index) & set(expr135_labels.index)
        print(f'\n对比验证: 与 135 区矩阵的公共皮层区 {len(common)} 个', flush=True)
        diff_max = 0.0
        for r in list(common)[:20]:
            d = (expr_cortex.loc[r].values - expr135_labels.loc[r].values)
            d = np.nanmax(np.abs(d)) if d.size else 0
            diff_max = max(diff_max, float(d))
        print(f'  前 20 公共区 max|diff| = {diff_max:.2e}（应≈0，同一 abagen 管线）', flush=True)
    except Exception as e:
        print('对比验证跳过:', e, flush=True)

    print('完成。', flush=True)


if __name__ == '__main__':
    main()
