# -*- coding: utf-8 -*-
"""Fig. 5 数据准备：通路富集。

产出（figures_20260915/）：
  fig5_es_curves.csv          5 条代表通路的 running ES 曲线（按 rank 抽样）
  fig5_volcano.csv            GO + KEGG 全部条目的 NES / FDR q（火山图）
  fig5_spin.csv               184 条被检验通路的 NES / FDR q / p_spin（存活标记）
  fig5_pathway_maps.csv       5 条代表通路的皮层表达图（区内中位 Z，99 区）
  fig5_pathway_map_r.csv      表达图 × LV1 loading 的空间 r（与正文对账）
  fig5_representatives.csv    代表通路汇总（NES / FDR / spin p / r）

对账基准（正文 / Table 3）：
  OXPHOS NES 2.79 FDR<0.001 spin 0.048 r +0.366
  Glycolysis 2.01 0.001 0.056 +0.308
  Ribosome 2.70 <0.001 0.022 +0.340
  Proteasome 2.41 <0.001 0.036 +0.385
  MHC class I 2.23 <0.001 0.020 +0.408   ← 2026-09-16 修正（原 +0.323 来自 GO:0002478）
  （spin 列为**表型层口径、同向单侧**，2026-10-09 行对齐质心更新）
"""
import csv
import io
import os

import numpy as np
from decimal import Decimal, ROUND_HALF_UP


def rnd(x, nd):
    """稿件口径：四舍五入到 nd 位小数"""
    q = Decimal(1).scaleb(-nd)
    return '%.*f' % (nd, Decimal(str(float(x))).quantize(q, rounding=ROUND_HALF_UP))

ROOT = 'F:/dcmtonii/run62_crosssectional'
AHBA = ROOT + '/PSF_ahba_100'
REF = ROOT + '/PSF_references/gene_sets'
W = ROOT + '/PSF_pls_behavior/断连稿写作_20260907'
HERE = W + '/figures_20260915'

REPS = [
    ('Oxidative phosphorylation', 'KEGG'),
    ('Glycolysis / Gluconeogenesis', 'KEGG'),
    ('Ribosome', 'KEGG'),
    ('Proteasome', 'KEGG'),
    ('antigen processing and presentation of exogenous peptide antigen '
     'via MHC class I, TAP-dependent (GO:0002479)', 'GO'),
]
EXPECT = {  # 正文 Table 3 值，用于 assert
    'Oxidative phosphorylation': (2.79, 0.048, 0.366),
    'Glycolysis / Gluconeogenesis': (2.01, 0.056, 0.308),
    'Ribosome': (2.70, 0.022, 0.340),
    'Proteasome': (2.41, 0.036, 0.385),
    # ⚠ 2026-09-16 修正：原写 0.323，该值实为 GO:0002478（103→93 基因）的表达图相关；
    # 本行声明的是 GO:0002479（73→66 基因）→ 同管线正确值 +0.408
    # （来源：run_sal_specificity_spin_v2_newmasks_20260915.py 修正后重跑，
    #   特征 Path_GO_MHC_I_antigen_TAP，r=+0.4076，spin p=0.001）
    'antigen processing and presentation of exogenous peptide antigen '
    'via MHC class I, TAP-dependent (GO:0002479)': (2.23, 0.020, 0.408),
}
SHORT = {REPS[0][0]: 'Oxidative phosphorylation', REPS[1][0]: 'Glycolysis / Gluconeogenesis',
         REPS[2][0]: 'Ribosome', REPS[3][0]: 'Proteasome', REPS[4][0]: 'MHC class I (GO:0002479)'}


def read_gmt(path):
    d = {}
    for line in io.open(path, encoding='utf-8'):
        parts = line.rstrip('\n').split('\t')
        if len(parts) < 3:
            continue
        d[parts[0]] = [g for g in parts[2:] if g]
    return d


# ---------------------------------------------------------------- 排序表
rank_tbl = list(csv.DictReader(io.open(
    AHBA + '/PSF_pls_gene_ranking_cortex100_residY_newmasks_20260915.csv',
    encoding='utf-8-sig')))
order = sorted(rank_tbl, key=lambda x: int(x['rank']))
genes = [x['gene'] for x in order]
load = np.array([float(x['loading']) for x in order])
N = len(genes)
gidx = {g: i for i, g in enumerate(genes)}
print('排序基因 %d 个' % N)

gmt = {'KEGG': read_gmt(REF + '/KEGG_2021_Human.gmt'),
       'GO': read_gmt(REF + '/GO_Biological_Process_2021.gmt')}


# ---------------------------------------------------------------- a) running ES
def running_es(gene_set, p=1.0):
    hit = np.zeros(N, bool)
    for g in gene_set:
        i = gidx.get(g)
        if i is not None:
            hit[i] = True
    nh = int(hit.sum())
    if nh == 0:
        return None, None, 0
    w = np.abs(load) ** p
    wr = np.where(hit, w, 0.0)
    tot = wr.sum()
    inc = wr / tot                       # 命中处步长
    dec = 1.0 / (N - nh)                 # 未命中处步长
    run = np.cumsum(np.where(hit, inc, -dec))
    return run, hit, nh


rows_es, rows_rep = [], []
for name, db in REPS:
    run, hit, nh = running_es(gmt[db][name])
    es_obs = float(run.max()) if abs(run.max()) >= abs(run.min()) else float(run.min())
    src = list(csv.DictReader(io.open(
        AHBA + '/gseapy_prerank_cortex100_residY_w1_newmasks_20260915_%s.csv' % db,
        encoding='utf-8-sig')))
    g_es = float([x for x in src if x['Term'] == name][0]['ES'])
    print('  %-42s 自算 ES %+.4f ｜ gseapy ES %+.4f ｜ 基因 %d'
          % (SHORT[name][:42], es_obs, g_es, nh))
    assert abs(es_obs - g_es) < 2e-3, '自算 ES 与 gseapy 不符：%s' % name
    # 抽样写盘（全曲线太长，按对数间隔取 900 点，末点补全）
    # q 为**1 基**秩；数组索引 = q-1
    q = np.unique(np.round(np.concatenate(
        [np.geomspace(1, N, 900), [N]])).astype(int))
    q = q[(q >= 1) & (q <= N)]
    for q_ in q:
        rows_es.append(dict(pathway=SHORT[name], rank=int(q_),
                            running_es=float(run[q_ - 1])))
    # 命中位置（画 tic 用，最多 500 个）
    hp = np.where(hit)[0]
    if len(hp) > 500:
        hp = hp[np.linspace(0, len(hp) - 1, 500).astype(int)]
    for i in hp:
        rows_es.append(dict(pathway=SHORT[name] + ' |HIT', rank=int(i + 1),
                            running_es=float(run[i])))
    rows_rep.append(dict(pathway=SHORT[name], db=db, n_genes=nh, es=es_obs))

with io.open(HERE + '/fig5_es_curves.csv', 'w', encoding='utf-8-sig', newline='') as f:
    w_ = csv.DictWriter(f, fieldnames=['pathway', 'rank', 'running_es'])
    w_.writeheader()
    w_.writerows(rows_es)
print('已写出 fig5_es_curves.csv（%d 行）' % len(rows_es))

# ---------------------------------------------------------------- b) 火山（全部条目）
vol = []
for db, path in [('KEGG', 'KEGG'), ('GO', 'GO')]:
    for x in csv.DictReader(io.open(
            AHBA + '/gseapy_prerank_cortex100_residY_w1_newmasks_20260915_%s.csv' % path,
            encoding='utf-8-sig')):
        vol.append(dict(db=db, term=x['Term'], nes=float(x['NES']),
                        fdr=float(x['FDR q-val']),
                        significant=float(x['FDR q-val']) < 0.05,
                        gene_set_size=int(x['Tag %'].split('/')[0].strip())
                        if '/' in str(x['Tag %']) else ''))
with io.open(HERE + '/fig5_volcano.csv', 'w', encoding='utf-8-sig', newline='') as f:
    w_ = csv.DictWriter(f, fieldnames=['db', 'term', 'nes', 'fdr', 'significant',
                                       'gene_set_size'])
    w_.writeheader()
    w_.writerows(vol)
n_sig_go = sum(1 for x in vol if x['db'] == 'GO' and x['significant'])
n_sig_kg = sum(1 for x in vol if x['db'] == 'KEGG' and x['significant'])
print('已写出 fig5_volcano.csv（GO %d 显著 / %d，KEGG %d 显著 / %d）'
      % (n_sig_go, sum(1 for x in vol if x['db'] == 'GO'),
         n_sig_kg, sum(1 for x in vol if x['db'] == 'KEGG')))
assert (n_sig_go, n_sig_kg) == (134, 50), '显著条数应与正文一致（134 / 50）'

# ---------------------------------------------------------------- c) spin 存活
spin = list(csv.DictReader(io.open(
    AHBA + '/permutation_newmasks_20260915/PSF_gsea_pathway_spin_official_result.csv',
    encoding='utf-8-sig')))
rows_spin = [dict(pathway=x['pathway'], gmt=x['gmt'], nes=float(x['NES']),
                  fdr=float(x['fdr']), p_spin=float(x['p_spin_one']),
                  p_spin_two=float(x['p_spin']),
                  survived=float(x['p_spin_one']) < 0.05) for x in spin]
with io.open(HERE + '/fig5_spin.csv', 'w', encoding='utf-8-sig', newline='') as f:
    w_ = csv.DictWriter(f, fieldnames=['pathway', 'gmt', 'nes', 'fdr', 'p_spin',
                                       'p_spin_two', 'survived'])
    w_.writeheader()
    w_.writerows(rows_spin)
alive = sum(1 for x in rows_spin if x['survived'])
print('已写出 fig5_spin.csv（存活 %d / %d）' % (alive, len(rows_spin)))
assert (len(rows_spin), alive) == (184, 127), '应为 184 条 / 存活 127（同向单侧）'

# ---------------------------------------------------------------- d) 代表通路表达图
exp = list(csv.reader(io.open(AHBA + '/PSF_ahba_100_expression.csv', encoding='utf-8-sig')))
exp_hdr = exp[0]
exp_regions = [r[0] for r in exp[1:]]
def _v(x):
    x = str(x).strip()
    return float(x) if x not in ('', 'nan', 'NA', 'NaN') else np.nan


E = np.array([[_v(v) for v in r[1:]] for r in exp[1:]])   # 空单元 → NaN（RH_SomMot_2 全 NaN）
col = {g: i for i, g in enumerate(exp_hdr[1:])}
print('表达矩阵 %d 区 × %d 基因' % E.shape)

sal = list(csv.DictReader(io.open(
    ROOT + '/PSF_pls_behavior/'
           'PSF_pls_discon_salience_residY_age_sex_cortex100_newmasks_20260915.csv',
    encoding='utf-8-sig')))
sal_map = {x['region']: float(x['salience']) for x in sal}
# 有效区 = salience 有值 ∩ 表达非全 NaN（RH_SomMot_2 在 AHBA 无采样 → 全 NaN）
valid = [i for i, r_ in enumerate(exp_regions)
         if r_ in sal_map and np.isfinite(E[i]).any()]
dropped = [exp_regions[i] for i in range(len(exp_regions))
           if exp_regions[i] in sal_map and not np.isfinite(E[i]).any()]
print('有效区 %d（排除全 NaN：%s）' % (len(valid), dropped))
assert len(valid) == 99, '应为 99 个有效区'

load_v = np.array([sal_map[exp_regions[i]] for i in valid])

# 权威 r 来源（稿件引用的输出文件；文件内为 4 位小数）
SRC_RESULT = (AHBA + '/permutation_newmasks_20260915/'
                     'PSF_sal_specificity_spin_result_newmasks_20260915.csv')
SRC_PROTEO = (AHBA + '/permutation_newmasks_20260915/'
                     'PSF_sal_specificity_spin_proteostasis_newmasks_20260915.csv')
_rsrc = {}
for _p in (SRC_RESULT, SRC_PROTEO):
    for _x in csv.DictReader(io.open(_p, encoding='utf-8-sig')):
        if _x['mode'] == 'salience':
            _rsrc[_x['feature']] = float(_x['r'])
FEAT_KEY = {REPS[0][0]: 'Path_KEGG_OXPHOS', REPS[1][0]: 'Path_KEGG_Glycolysis',
            REPS[2][0]: 'Path_KEGG_Ribosome', REPS[3][0]: 'Path_KEGG_Proteasome',
            REPS[4][0]: 'Path_GO_MHC_I_antigen_TAP'}

rows_map, rows_r = [], []
for name, db in REPS:
    gs = [g for g in gmt[db][name] if g in col]
    sub = E[np.ix_(valid, [col[g] for g in gs])]
    # 每基因跨区 Z（忽略 NaN）；再取区内跨基因**中位** Z
    mu = np.nanmean(sub, axis=0)
    sd = np.nanstd(sub, axis=0, ddof=0)
    z = (sub - mu) / sd
    m = np.nanmedian(z, axis=1)
    r_calc = float(np.corrcoef(m, load_v)[0, 1])
    r_src = _rsrc[FEAT_KEY[name]]
    assert abs(r_calc - r_src) < 1e-3, \
        '重算 r 与权威 CSV 不符：%s（%.5f vs %.4f）' % (name, r_calc, r_src)
    r = r_src                      # 输出用权威文件值（与稿件同源）
    rows_r.append(dict(pathway=SHORT[name], n_genes=len(gs), r=r,
                       n_parcels=len(valid), r_recomputed=round(r_calc, 5)))
    for reg, v in zip([exp_regions[i] for i in valid], m):
        rows_map.append(dict(pathway=SHORT[name], region=reg, median_z=float(v)))
    print('  %-42s 基因 %d ｜ r(权威) = %+.4f ｜ 重算 %+.5f ｜ 正文 %+.3f'
          % (SHORT[name][:42], len(gs), r, r_calc, EXPECT[name][2]))
    assert rnd(r, 3) == rnd(EXPECT[name][2], 3), '表达图 r 与正文不符：%s' % name

with io.open(HERE + '/fig5_pathway_maps.csv', 'w', encoding='utf-8-sig', newline='') as f:
    w_ = csv.DictWriter(f, fieldnames=['pathway', 'region', 'median_z'])
    w_.writeheader()
    w_.writerows(rows_map)
with io.open(HERE + '/fig5_pathway_map_r.csv', 'w', encoding='utf-8-sig', newline='') as f:
    w_ = csv.DictWriter(f, fieldnames=['pathway', 'n_genes', 'r', 'n_parcels',
                                       'r_recomputed'])
    w_.writeheader()
    w_.writerows(rows_r)
print('已写出 fig5_pathway_maps.csv / fig5_pathway_map_r.csv')

# ---------------------------------------------------------------- 代表通路汇总
sp_map = {x['pathway']: x for x in spin}
rep_rows = []
for name, db in REPS:
    key = SHORT[name]
    src = list(csv.DictReader(io.open(
        AHBA + '/gseapy_prerank_cortex100_residY_w1_newmasks_20260915_%s.csv' % db,
        encoding='utf-8-sig')))
    row = [x for x in src if x['Term'] == name][0]
    r_map = [x for x in rows_r if x['pathway'] == key][0]
    rep_rows.append(dict(pathway=key, db=db,
                         label='GO:0002479' if db == 'GO' else key,
                         n_genes=r_map['n_genes'],
                         nes=float(row['NES']), fdr=float(row['FDR q-val']),
                         p_spin=float(sp_map[name]['p_spin_one']),
                         r=r_map['r']))
with io.open(HERE + '/fig5_representatives.csv', 'w', encoding='utf-8-sig',
             newline='') as f:
    w_ = csv.DictWriter(f, fieldnames=['pathway', 'db', 'label', 'n_genes', 'nes',
                                       'fdr', 'p_spin', 'r'])
    w_.writeheader()
    w_.writerows(rep_rows)
print('已写出 fig5_representatives.csv')
for z in rep_rows:
    print('   %-42s NES %.2f FDR %.4f spin %.3f r %+.3f'
          % (z['label'][:42], z['nes'], z['fdr'], z['p_spin'], z['r']))

# 与 Table 3 逐项对账
for z in rep_rows:
    full = [k for k in EXPECT if SHORT[k] == z['pathway']][0]
    e_nes, e_spin, e_r = EXPECT[full]
    # 稿件用**半进位**报到 2/3 位小数；Python round() 是银行家舍入且受二进制表示影响
    # （round(0.3655,3)→0.365，稿件为 0.366）→ 统一用 Decimal ROUND_HALF_UP。
    assert rnd(z['nes'], 2) == rnd(e_nes, 2), z['pathway']
    assert rnd(z['p_spin'], 3) == rnd(e_spin, 3), z['pathway']
    assert rnd(z['r'], 3) == rnd(e_r, 3), z['pathway']
print('\n✓ 5 条代表通路的 NES / spin p / r 与 Table 3 完全一致')
