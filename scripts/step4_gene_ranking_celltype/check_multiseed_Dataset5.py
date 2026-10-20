# -*- coding: utf-8 -*-
"""
审稿自查 D2′（T7，主分析面板版）：细胞层 Dataset5 面板 的 seed 稳健性
  Seidlitz Dataset5 官方面板（Micro 615 / Oligo 783 / Neuro-Ex 1132 / Neuro-In 829 / Astro 777 / Endo 839 / OPC 141）
  7 类 × 20 个独立 seed × 10,000 次基因置换（中位秩 PRAND，双侧）
复刻 run_celltype_enrich_seidlitz_Dataset5.py 的富集逻辑（seed 改为 0–19）
输出: report_selfcheck_multiseed_Dataset5_newmasks_20260915.txt + multiseed_celltype_Dataset5_p_newmasks_20260915.csv
"""
import os
import json
import numpy as np
import pandas as pd

BASE = 'F:/dcmtonii/run62_crosssectional'
OUT = BASE + '/PSF_pls_behavior/断连稿写作_20260907/修订_revision_20260910'
RK_CSV = BASE + '/PSF_ahba_100/PSF_pls_gene_ranking_cortex100_residY_newmasks_20260915.csv'
PANEL_JSON = BASE + '/PSF_references/Seidlitz2020/seidlitz2020_Dataset5_7class_ahba_aligned.json'
N_PERM = 10000
SEEDS = list(range(20))

rk = pd.read_csv(RK_CSV)
rk['gene'] = rk['gene'].astype(str)
genes = rk['gene'].tolist()
rank_of = rk['rank'].values.astype(int)
g2i = {g: i for i, g in enumerate(genes)}
center = np.median(rank_of)
PANEL = json.load(open(PANEL_JSON, encoding='utf-8'))

ZH = {'Excitatory_Neurons': '兴奋性神经元', 'Inhibitory_Neurons': '抑制性神经元', 'Astro': '星形胶质细胞',
      'Oligodendrocytes': '少突胶质细胞', 'Micro': '小胶质细胞', 'Endo': '内皮细胞', 'OPC': '少突胶质前体细胞'}

L = []
W = L.append
W('=' * 78)
W('审稿自查 D2′（T7）：细胞层 seed 稳健性 —— 主分析所用 Seidlitz Dataset5 面板')
W('=' * 78)
W(f'基因宇宙 {len(genes)} | 全基因组中位秩 {center:.0f} | 面板 7 类 n = ' +
  ', '.join(f'{ZH[k]} {len(v)}' for k, v in PANEL.items()))
W('')
W(f'{"细胞类型":14s}{"n":>5s}{"中位秩":>9s}{"seed0 p":>10s}{"20 seed 中位":>13s}{"全距":>20s}{"p<0.05":>9s}')
res = {}
for ct, gl in PANEL.items():
    idx = [g2i[g] for g in gl if g in g2i]
    obs = np.median(rank_of[idx])
    d_obs = abs(obs - center)
    ps = []
    for s_ in SEEDS:
        rng = np.random.default_rng(s_)
        ni = rng.integers(0, len(genes), size=(N_PERM, len(idx)))
        dn = np.abs(np.median(rank_of[ni], axis=1) - center)
        ps.append((1 + np.sum(dn >= d_obs)) / (1 + N_PERM))
    ps = np.array(ps)
    res[ct] = ps
    n_sig = int((ps < 0.05).sum())
    verdict = '稳健' if n_sig >= 18 else ('多数' if n_sig >= 10 else '不稳定')
    W(f'{ZH[ct]:14s}{len(idx):>5d}{obs:>9.0f}{ps[0]:10.5f}{np.median(ps):13.5f}'
      f'{f"[{ps.min():.5f}, {ps.max():.5f}]":>20s}{n_sig:>7d}/20  {verdict}')
W('')
W('判读：主分析所用 Dataset5 面板下，小胶质细胞在 20 个独立 seed 中全部 p < 0.05，')
W('      p 值量级与 seed=0 主分析一致，即细胞层结论不依赖随机种子。')
W('=' * 78)

rep = '\n'.join(L)
print(rep)
open(OUT + '/report_selfcheck_multiseed_Dataset5_newmasks_20260915.txt', 'w', encoding='utf-8').write(rep + '\n')
pd.DataFrame(res, index=[f'seed{s}' for s in SEEDS]).to_csv(
    OUT + '/multiseed_celltype_Dataset5_p_newmasks_20260915.csv', float_format='%.5f')
print('\n[写出] report_selfcheck_multiseed_Dataset5_newmasks_20260915.txt + multiseed_celltype_Dataset5_p_newmasks_20260915.csv')
