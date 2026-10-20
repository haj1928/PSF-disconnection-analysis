# -*- coding: utf-8 -*-
"""补充图数据准备（2026-09-16）—— 图 S1 与四张已移出稿件的补充图

S1  受体 21×21 空间相关矩阵（源：修订_revision_20260910/receptor_21x21_correlation.csv）
S2  21 张受体图谱皮层分区值（源：figures_20260915/unused_receptor_parcel_values.csv）
S4  LV1 loading 与临床/认知量表的相关（源：supplement/newmasks_20260915/S6a_...csv）
S5  疾病通路 × 氧化磷酸化基因集的重叠（源：KEGG GMT 现算）
S6  病灶体积守恒 QC（源：supplement/newmasks_20260915/volume_conservation_QC_N124_...csv）

注：Fig. S3（断连矩阵热图）需先由**新掩膜** LQT 输出汇总 124×100 矩阵，
    现有 CSV 为 20260903 旧口径，故本轮不出图（见报告）。
"""
import csv
import io
import os
import sys

ROOT = 'F:/dcmtonii/run62_crosssectional'
W = ROOT + '/PSF_pls_behavior/断连稿写作_20260907'
REV = W + '/修订_revision_20260910'
SUP = REV + '/supplement/newmasks_20260915'
HERE = W + '/figures_20260915'
GMT = ROOT + '/PSF_references/gene_sets'


def wr(name, rows, fields=None):
    fld = fields or list(rows[0].keys())
    with io.open(HERE + '/' + name, 'w', encoding='utf-8-sig', newline='') as f:
        w_ = csv.DictWriter(f, fieldnames=fld)
        w_.writeheader()
        w_.writerows(rows)
    print('  已写出 %s（%d 行）' % (name, len(rows)))


# ---------------- S1：21×21 相关矩阵 ----------------
print('S1 受体 21×21 相关矩阵')
raw = list(csv.reader(io.open(REV + '/receptor_21x21_correlation.csv', encoding='utf-8-sig')))
head = raw[0][1:]
mat = []
for row in raw[1:]:
    if not row or not row[0]:
        continue
    for j, v in enumerate(row[1:1 + len(head)]):
        mat.append(dict(atlas=row[0], other=head[j], r=float(v)))
wr('figS1_receptor_corr_long.csv', mat)
print('  矩阵 %d × %d' % (len(head), len(head)))
# 附图用宽表
rows = []
for row in raw[1:]:
    if not row or not row[0]:
        continue
    d = {'atlas': row[0]}
    for j, v in enumerate(row[1:1 + len(head)]):
        d[head[j]] = float(v)
    rows.append(d)
wr('figS1_receptor_corr.csv', rows, ['atlas'] + head)

# ---------------- S2：21 张图谱分区值 ----------------
print('S2 受体图谱皮层分布')
s2 = list(csv.DictReader(io.open(HERE + '/unused_receptor_parcel_values.csv',
                                 encoding='utf-8-sig')))
assert len(s2) == 21, len(s2)
wr('figS1_meta.csv', [dict(atlas=x['atlas'], system=x['neurotransmitter'],
                           target=x['target'], res=x['res']) for x in s2])

# ---------------- unused：LV1 vs 临床量表（已移出稿件） ----------------
print('S4 LV1 与临床量表相关')
s4 = list(csv.DictReader(io.open(SUP + '/S6a_LV1_vs_clinical_scales.csv',
                                 encoding='utf-8-sig')))
NICE = {'FSS-9 total': 'FSS-9 total', 'HAMD': 'HAMD', 'HAMA': 'HAMA',
        'AES-C': 'AES-C (apathy)', 'MMSE': 'MMSE', 'MoCA': 'MoCA',
        'NIHSS': 'NIHSS', 'mRS': 'mRS', 'BI': 'Barthel index',
        'Age': 'Age'}
out = [dict(scale=NICE.get(x['scale'], x['scale']), n=int(x['n']),
            rho=float(x['spearman_rho_vs_LV1loading']), p=float(x['p']))
       for x in s4]
out.sort(key=lambda z: -z['rho'])
wr('unused_clinical_scales.csv', out)

# ---------------- unused：疾病通路 × OXPHOS 重叠（已移出稿件） ----------------
print('S5 疾病通路 × OXPHOS 基因重叠')


def load_gmt(fn):
    d = {}
    for line in io.open(GMT + '/' + fn, encoding='utf-8', errors='ignore'):
        p = line.rstrip('\n').split('\t')
        if len(p) >= 3:
            d[p[0]] = [g for g in p[2:] if g]
    return d


kegg = load_gmt('KEGG_2021_Human.gmt')
oxphos = set(kegg['Oxidative phosphorylation'])
DISEASE = ['Thermogenesis', 'Huntington disease', 'Parkinson disease',
           'Prion disease', 'Alzheimer disease']
rows = []
for name in DISEASE + ['Proteasome', 'Ribosome']:
    gs = set(kegg[name])
    ov = len(gs & oxphos)
    rows.append(dict(pathway=name, n_genes=len(gs), n_overlap=ov,
                     pct_overlap=round(100.0 * ov / len(gs), 1)))
rows.sort(key=lambda z: -z['pct_overlap'])
wr('unused_disease_oxphos_overlap.csv', rows)
inter = set()
for name in ['Thermogenesis', 'Huntington disease', 'Parkinson disease', 'Prion disease']:
    inter = set(kegg[name]) if not inter else (inter & set(kegg[name]))
print('    四病通路交集 %d 个基因；全部属于 OXPHOS：%s'
      % (len(inter), inter <= oxphos))

# ---------------- unused：体积守恒 QC（已移出稿件） ----------------
print('S6 体积守恒 QC')
s6 = list(csv.DictReader(io.open(SUP + '/volume_conservation_QC_N124_newmasks_20260915.csv',
                                encoding='utf-8-sig')))
print('    列名：', list(s6[0].keys()))
wr('unused_volume_qc.csv', s6)
print('\n完成。')
