# -*- coding: utf-8 -*-
"""
审稿自查 C：细胞层（小胶质）与基因程序层（线粒体能量代谢）是否用的是同一批基因？
（不改动任何原分析，仅对既有排序产物做集合运算）

输入: PSF_ahba_100/PSF_pls_gene_ranking_cortex100_residY_newmasks_20260915.csv（15,633 基因 + PLS 载荷排序）
      标记面板 = run_celltype_enrich_seidlitz_prand_cortex100_residY.py 内的 66 marker（7 类）
输出: report_selfcheck_gene_overlap_newmasks_20260915.txt
"""
import numpy as np
import pandas as pd

BASE = 'F:/dcmtonii/run62_crosssectional'
OUT = BASE + '/PSF_pls_behavior/断连稿写作_20260907/修订_revision_20260910'

MARKERS = {
    'Excitatory_Neurons': ['SLC17A7', 'SATB2', 'NEFL', 'NEFM', 'CAMK2A', 'RORB', 'CUX2', 'FEZF2', 'TBR1', 'NRGN'],
    'Inhibitory_Neurons': ['GAD1', 'GAD2', 'SLC32A1', 'PVALB', 'SST', 'VIP', 'LAMP5', 'RELN', 'CALB1', 'CALB2'],
    'Astrocytes': ['GFAP', 'AQP4', 'SLC1A2', 'SLC1A3', 'ALDH1L1', 'S100B', 'GJA1', 'SOX9', 'ALDOC'],
    'Oligodendrocytes': ['MBP', 'PLP1', 'MOG', 'MAG', 'CNP', 'MOBP', 'CLDN11', 'OPALIN', 'ST18'],
    'OPC': ['PDGFRA', 'CSPG4', 'OLIG1', 'OLIG2', 'SOX10', 'GPR17', 'LHFPL3'],
    'Microglia': ['CX3CR1', 'CSF1R', 'AIF1', 'TMEM119', 'PTPRC', 'ITGAM', 'CD68', 'TYROBP', 'P2RY12', 'TREM2'],
    'Endothelial': ['CLDN5', 'FLT1', 'PECAM1', 'VWF', 'ABCB1', 'SLC2A1', 'CDH5', 'TEK'],
}
MITO_CORE = ['AMDHD1', 'FAM96B', 'GALNT1', 'HMOX2', 'CRYBB1', 'DNAJC15', 'PRDX5', 'FIS1', 'RIT2',
             'NOSIP', 'NDUFA13', 'RPS19BP1', 'MRPL53', 'SSNA1', 'FABP6']

rk = pd.read_csv(BASE + '/PSF_ahba_100/PSF_pls_gene_ranking_cortex100_residY_newmasks_20260915.csv')
rk['gene'] = rk['gene'].astype(str)
rank_of = dict(zip(rk['gene'], rk['rank']))
N = len(rk)
top300 = set(rk.nsmallest(300, 'rank')['gene'])
bot300 = set(rk.nlargest(300, 'rank')['gene'])

L = []
W = L.append
W('=' * 78)
W('审稿自查 C：细胞层（小胶质）与代谢程序层是否共用同一批基因？')
W('=' * 78)
W(f'基因库 N = {N}（abagen 99 皮层区 × 15,633 基因；PLS w1 排序）')
W('')

W('[C1] 7 类细胞标记面板的内部重叠（66 个 marker 是否互斥）')
allm = [g for v in MARKERS.values() for g in v]
W(f'  marker 总数 {len(allm)} | 去重后 {len(set(allm))} -> {"互斥" if len(allm) == len(set(allm)) else "存在重叠"}')
ov = []
ks = list(MARKERS)
for i in range(len(ks)):
    for j in range(i + 1, len(ks)):
        s = set(MARKERS[ks[i]]) & set(MARKERS[ks[j]])
        if s:
            ov.append((ks[i], ks[j], s))
W(f'  两两交集: {ov if ov else "无"}')
W('')

W('[C2] 小胶质 10 个 marker 在排序中的位次')
micro = MARKERS['Microglia']
for g in sorted(micro, key=lambda x: rank_of.get(x, 10 ** 9)):
    W(f'  {g:10s} rank = {rank_of.get(g, "未在库中")}')
W(f'  落在正共变 top300 的: {sorted(set(micro) & top300)}')
W(f'  落在负共变 bottom300 的: {sorted(set(micro) & bot300)}')
W('')

W('[C3] 代谢/线粒体核心基因（top15）在排序中的位次')
for g in MITO_CORE:
    W(f'  {g:10s} rank = {rank_of.get(g, "未在库中")}')
W('')

W('[C4] 两套基因的直接交集')
inter = set(micro) & set(MITO_CORE)
W(f'  小胶质 marker(10) ∩ 代谢核心(top15) = {inter if inter else "空集"}')
inter300 = set(micro) & top300
W(f'  小胶质 marker(10) ∩ 正共变 top300 = {sorted(inter300)}（占 10 个 marker 的 {100*len(inter300)/10:.0f}%）')
W(f'  代谢核心 top15 ∩ 正共变 top300 = {len(set(MITO_CORE) & top300)} / 15')
W('')

W('[C5] 判读')
W('  - 若两套基因交集为空、且小胶质 marker 与代谢核心基因位次明显错开')
W('    -> 细胞层与基因程序层由不同基因承载，属两条独立证据，')
W('       可驳回"代谢核心只是小胶质信号的另一种呈现"的质疑')
W('  - 仍需在正文补一句完整性声明：两层使用彼此独立的基因集（给出 marker 清单出处）')
W('=' * 78)

rep = '\n'.join(L)
print(rep)
with open(OUT + '/report_selfcheck_gene_overlap_newmasks_20260915.txt', 'w', encoding='utf-8') as f:
    f.write(rep + '\n')
print('\n[写出] report_selfcheck_gene_overlap_newmasks_20260915.txt')
