# -*- coding: utf-8 -*-
"""补充图绘图（2026-09-16）

Fig. S1  受体 21×21 空间相关矩阵热图（128 mm）—— 稿件在用
（以下四张稿件未用，EMIT_UNUSED=False 时不导出）
  21 张受体图谱的皮层分布（190 mm，7×3 网格）
  LV1 loading 与临床/认知量表相关（90 mm 棒棒糖）
  疾病相关 KEGG 通路 × 氧化磷酸化基因集重叠（90 mm 堆叠横条）
  病灶体积守恒质量控制（90 mm 散点）

数据全部读 CSV（`data_figS_20260916.py` 产出）。
"""
import csv
import io
import os
import sys

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
from nilearn import datasets, plotting

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import figstyle_20260916 as fs
import fig2_behavioral_pls_20260916 as f3

fs.apply_style(7.0)
mpl.rcParams['savefig.bbox'] = None
HERE = os.path.dirname(os.path.abspath(__file__))


def rd(fn):
    return list(csv.DictReader(io.open(os.path.join(HERE, fn), encoding='utf-8-sig')))


# 2026-09-17（方案 B）：这两张图已从稿件移出（正文零引用）。置 True 可重新导出到 _unused_20260917/。
EMIT_UNUSED = False


UNUSED_DIR = os.path.join(HERE, '_unused_20260917')     # 已移出稿件的图统一落这里


def save(fig, name, w_mm, h_mm):
    assert abs(fig.get_size_inches()[0] * 25.4 - w_mm) < 0.2
    out = UNUSED_DIR if name.startswith('unused_') else HERE
    os.makedirs(out, exist_ok=True)
    fs.save_all(fig, name, outdir=out, formats=('pdf', 'svg', 'eps', 'png'),
                bbox_inches=None)
    from PIL import Image
    w_px, h_px = Image.open(os.path.join(out, name + '.png')).size
    print('  %s → %.1f × %.1f mm（目标 %.1f × %.1f）'
          % (name, w_px / 300 * 25.4, h_px / 300 * 25.4, w_mm, h_mm))
    plt.close(fig)


SYS_CN = {'Acetylcholine': 'ACh', 'Opioid': 'Opioid', 'GABA': 'GABA',
          'Dopamine': 'DA', 'Serotonin': '5-HT', 'Noradrenaline': 'NA',
          'Glutamate': 'Glu', 'Synaptic': 'SV'}
SHORT = {'a4b2-nACh': 'α4β2', 'mu-opioid': 'μ-opioid', 'GABA-A': 'GABA-A',
         'D1': 'D1', '5-HT1A': '5-HT1A', 'NET': 'NET', 'GABA-A-a5': 'GABA-A α5',
         'DAT': 'DAT', '5-HT1A(2)': '5-HT1A (2)', 'M1 mAChR': 'M1',
         'D2/D3': 'D2/D3', '5-HT2A': '5-HT2A', '5-HT1B': '5-HT1B', '5-HT4': '5-HT4',
         'VAChT': 'VAChT', '5-HTT': '5-HT', 'mGluR5': 'mGluR5',
         'D2/D3(2)': 'D2/D3 (2)', '5-HT2A(2)': '5-HT2A (2)', 'SV2A': 'SV2A',
         'NMDA': 'NMDA'}

# ================= S1：21×21 相关矩阵 =================
print('S1 相关矩阵热图')
wide = rd('figS1_receptor_corr.csv')
meta = {x['atlas']: x for x in rd('figS1_meta.csv')}
order = [x['atlas'] for x in wide]
labels = ['%s\n%s' % (SHORT.get(meta[a]['target'], meta[a]['target']),
                      SYS_CN.get(meta[a]['system'], meta[a]['system'])) for a in order]
M = np.array([[float(x[a]) for a in order] for x in wide])
assert M.shape == (21, 21)

fig = plt.figure(figsize=(fs.W_ONEHALF * fs.MM, 128 * fs.MM))
axm = fig.add_axes([0.175, 0.145, 0.665, 0.70])
im = axm.imshow(M, cmap=fs.CMAP_DIV, vmin=-1, vmax=1)
axm.set_xticks(range(21))
axm.set_yticks(range(21))
axm.set_xticklabels(labels, fontsize=4.4, rotation=90)
axm.set_yticklabels(labels, fontsize=4.4)
axm.tick_params(length=0)
for sp in axm.spines.values():
    sp.set_linewidth(0.5)
# 系统分块界线
sysseq = [meta[a]['system'] for a in order]
for i in range(1, 21):
    if sysseq[i] != sysseq[i - 1]:
        axm.axhline(i - 0.5, color=fs.BLACK, lw=0.5)
        axm.axvline(i - 0.5, color=fs.BLACK, lw=0.5)
cax = fig.add_axes([0.855, 0.30, 0.022, 0.40])
cb = fig.colorbar(im, cax=cax)
cb.set_label('Pearson r', fontsize=6.5)
cb.ax.tick_params(labelsize=6, width=0.6, length=2)
cb.outline.set_linewidth(0.5)
fig.text(0.175, 0.885, 'a', fontsize=9, fontweight='bold')
fig.text(0.175, 0.055, 'Atlas identity (target / neurotransmitter system); '
         'black lines separate systems', fontsize=6, color=fs.GREY)
# PC1 方差占比 = 相关矩阵最大特征值 / 迹（= 高特异度算法；勿用 eigvalsh(corrcoef(M))）
_ev = np.linalg.eigvalsh(M)
print('  PC1 方差占比：%.1f%%' % (100 * _ev[-1] / np.trace(M)))
assert abs(100 * _ev[-1] / np.trace(M) - 41.7) < 0.05, 'PC1 与图注 41.7%% 不符'
save(fig, 'FigS1_receptor_corr_matrix', fs.W_ONEHALF, 128.0)

# ================= unused：21 张图谱皮层分布（原补充图，已移出稿件） =================
print('S2 受体图谱皮层分布（21 张 × 2 半球，渲染中…）')
zrows = rd('unused_receptor_parcel_values.csv')   # 原始值；绘图时**逐图内部** z 化
parcels = [k for k in zrows[0].keys()
           if k.startswith('LH_') or k.startswith('RH_')]
assert len(parcels) == 100, len(parcels)
vmap = f3.annot_vertex_map()
fsa = datasets.fetch_surf_fsaverage('fsaverage5')
BG = {}
for hemi in ('left', 'right'):
    sulc = fsa['sulc_%s' % hemi]
    aparc = f3.BASE + '/PSF_repro_PET/assets/fsa5_%s_aparc.annot' % ('lh' if hemi == 'left' else 'rh')
    BG[hemi] = f3.clean_bgmap(sulc, aparc)
import nibabel as nb
n_v = {}
for hemi, hk in (('left', 'lh'), ('right', 'rh')):
    n_v[hemi] = nb.freesurfer.io.read_annot(
        os.path.join(f3.ANNOT, '%s.Schaefer2018_100Parcels_7Networks_order.annot'
                     % hk))[0].shape[0]

fig = plt.figure(figsize=(fs.W_DOUBLE * fs.MM, 186 * fs.MM))
for k, z in enumerate(zrows):
    vals = np.array([float(z[p]) for p in parcels])
    # 逐图内部 z 化（图注口径：各图独立色标）；源 CSV 的 *_z 列是跨图谱全局 z 化，
    # 会使 GABA-A 等绝对量大的图谱整体饱和（2026-09-16 修复）
    vals = (vals - np.nanmean(vals)) / np.nanstd(vals)
    r_, c_ = divmod(k, 3)
    for j, hemi in enumerate(('left', 'right')):
        hk = 'LH' if hemi == 'left' else 'RH'
        sub_idx = [i for i, p in enumerate(parcels) if p.startswith(hk)]
        sub_names = [parcels[i] for i in sub_idx]
        sub_val = vals[sub_idx]
        vm = f3.vertex_map(sub_names, sub_val, vmap, n_v[hemi])
        ax = fig.add_axes([(0.055 + c_ * 0.308 + j * 0.138),
                           (0.935 - (r_ + 1) * 0.1285),
                           0.135, 0.125], projection='3d')
        plotting.plot_surf_stat_map(
            surf_mesh=fsa['infl_%s' % hemi], stat_map=vm, bg_map=BG[hemi],
            hemi=hemi, view='lateral', cmap=fs.CMAP_DIV, vmin=-2.5, vmax=2.5,
            threshold=0, symmetric_cbar=True, colorbar=False, axes=ax, figure=fig)
        for _c in ax.collections:
            _c.set_rasterized(True)
    axl = fig.add_axes([0.055 + c_ * 0.308, 0.935 - (r_ + 1) * 0.1285 + 0.118, 0.27, 0.012])
    axl.axis('off')
    axl.text(0.5, 0.5, '%s (%s, %s)' % (SHORT.get(z['target'], z['target']),
                                        z['neurotransmitter'], z['res']),
             fontsize=5.8, ha='center', va='center', color=fs.BLACK)
cax = fig.add_axes([0.30, 0.035, 0.40, 0.012])
cb = fig.colorbar(plt.cm.ScalarMappable(norm=mpl.colors.Normalize(-2.5, 2.5),
                                        cmap=fs.CMAP_DIV), cax=cax, orientation='horizontal')
cb.set_label('Regional mean (z)', fontsize=6.5)
cb.ax.tick_params(labelsize=6, width=0.6, length=2)
cb.outline.set_linewidth(0.5)
if EMIT_UNUSED:
    save(fig, 'unused_receptor_cortical_maps', fs.W_DOUBLE, 186.0)

# ================= unused：LV1 vs 临床量表（已移出稿件） =================
print('S4 LV1 与临床量表')
d4 = rd('unused_clinical_scales.csv')
fig = plt.figure(figsize=(fs.W_SINGLE * fs.MM, 82 * fs.MM))
ax = fig.add_axes([0.34, 0.135, 0.61, 0.83])
ys = np.arange(len(d4))[::-1]
for y_, z in zip(ys, d4):
    rho = float(z['rho'])
    p_ = float(z['p'])
    sig = p_ < 0.05
    c_ = (fs.ACCENT_RED if rho > 0 else fs.CATEGORICAL[0]) if sig else fs.GREY
    ax.plot([0, rho], [y_, y_], '-', color=c_, lw=0.8, alpha=1.0 if sig else 0.55)
    ax.plot([rho], [y_], 'o', ms=3.2, mfc=c_ if sig else 'white',
            mec=c_, mew=0.7)
    ax.text(rho + (0.022 if rho >= 0 else -0.022), y_, '%.3f' % rho,
            fontsize=5.8, va='center', ha='left' if rho >= 0 else 'right',
            color=fs.BLACK if sig else fs.GREY)
ax.axvline(0, color=fs.BLACK, lw=0.6)
ax.set_yticks(ys)
ax.set_yticklabels(['%s (n = %s)' % (z['scale'], z['n']) for z in d4], fontsize=6.2)
ax.set_xlim(-0.48, 0.42)
ax.set_ylim(-0.7, len(d4) - 0.3)
ax.set_xlabel('Spearman rho with LV1 loading', fontsize=7)
ax.tick_params(axis='x', labelsize=6.5, width=0.6, length=2)
ax.tick_params(axis='y', length=0)
for sp in ('top', 'right', 'left'):
    ax.spines[sp].set_visible(False)
fig.text(0.34, 0.955, 'filled = p < 0.05', fontsize=6, color=fs.GREY)
if EMIT_UNUSED:
    save(fig, 'unused_clinical_scales', fs.W_SINGLE, 82.0)

# ================= unused：疾病通路 × OXPHOS 重叠（已移出稿件） =================
print('S5 疾病通路重叠')
d5 = rd('unused_disease_oxphos_overlap.csv')
fig = plt.figure(figsize=(fs.W_SINGLE * fs.MM, 66 * fs.MM))
ax = fig.add_axes([0.42, 0.20, 0.54, 0.74])
ys = np.arange(len(d5))[::-1]
for y_, z in zip(ys, d5):
    ov, tot = int(z['n_overlap']), int(z['n_genes'])
    ax.barh(y_, ov, height=0.62, color=fs.ACCENT_RED, edgecolor='none')
    ax.barh(y_, tot - ov, left=ov, height=0.62, color=fs.GREY, alpha=0.35,
            edgecolor='none')
    ax.text(tot + 6, y_, '%d / %d (%.0f%%)' % (ov, tot, float(z['pct_overlap'])),
            fontsize=5.8, va='center', ha='left', color=fs.BLACK)
ax.set_yticks(ys)
ax.set_yticklabels([z['pathway'] for z in d5], fontsize=6.2)
ax.set_xlim(0, 430)
ax.set_ylim(-0.7, len(d5) - 0.3)
ax.set_xlabel('Genes shared with oxidative phosphorylation', fontsize=6.8)
ax.tick_params(axis='x', labelsize=6.5, width=0.6, length=2)
ax.tick_params(axis='y', length=0)
for sp in ('top', 'right', 'left'):
    ax.spines[sp].set_visible(False)
fig.text(0.42, 0.955, 'red = shared with OXPHOS, grey = pathway-specific',
         fontsize=6, color=fs.GREY)
if EMIT_UNUSED:
    save(fig, 'unused_disease_oxphos_overlap', fs.W_SINGLE, 66.0)

# ================= unused：体积守恒 QC（已移出稿件） =================
print('S6 体积守恒 QC')
d6 = rd('unused_volume_qc.csv')
nat = np.array([float(x['native_mL']) for x in d6])
mni = np.array([float(x['mni_mL']) for x in d6])
fig = plt.figure(figsize=(fs.W_SINGLE * fs.MM, 88 * fs.MM))
ax = fig.add_axes([0.155, 0.145, 0.80, 0.80])
lim = (0.04, 200)
ax.plot(lim, lim, '-', color=fs.BLACK, lw=0.6, zorder=2)
ax.plot(lim, [lim[0] * 2, lim[1] * 2], ':', color=fs.GREY, lw=0.6, zorder=2)
ax.plot(lim, [lim[0] / 2, lim[1] / 2], ':', color=fs.GREY, lw=0.6, zorder=2)
ax.plot(nat, mni, 'o', ms=2.2, mfc=fs.CATEGORICAL[0], mec='none', alpha=0.75, zorder=3)
ax.set_xscale('log')
ax.set_yscale('log')
ax.set_xlim(*lim)
ax.set_ylim(0.06, 300)
ax.set_xlabel('Native-space lesion volume (mL)', fontsize=7)
ax.set_ylabel('MNI152-space lesion volume (mL)', fontsize=7)
ax.tick_params(labelsize=6.5, width=0.6, length=2)
ax.text(0.06, 0.93, 'solid = identity, dotted = 0.5x / 2x tolerance\n'
        'median MNI/native ratio = 1.27 (n = 124)',
        transform=ax.transAxes, fontsize=6, color=fs.GREY, va='top')
for sp in ('top', 'right'):
    ax.spines[sp].set_visible(False)
if EMIT_UNUSED:
    save(fig, 'unused_volume_conservation_QC', fs.W_SINGLE, 88.0)
print('\n完成：稿件在用 = 图 S1（受体相关矩阵）；'
      '已移出稿件 = 21 张图谱皮层分布、临床量表、疾病通路×OXPHOS、体积守恒 QC 四张'
      '（EMIT_UNUSED=True 可重新导出）。')
