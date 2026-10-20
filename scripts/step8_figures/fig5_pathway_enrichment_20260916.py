# -*- coding: utf-8 -*-
"""Fig. 5 通路富集（Elsevier 双栏 190 mm）。

a  GSEA running-ES 曲线（5 条代表通路，含命中 rug）
b  NES–FDR 火山（全部 GO 5,301 + KEGG 312 条；显著者实心）
c  184 条受检验通路的 NES vs 表型层 spin p（存活 127 / 淘汰 57）；
   5 条代表通路以 5 色标出（2026-09-18 由 b 移入，b 的坐标系无法分辨它们）
d  5 条代表通路的皮层表达图（区内中位 Z，左右半球外侧面，共享色标）

数据全部读 CSV（`data_fig5_pathway_20260916.py` 产出）。
表面渲染复用 Fig. 3 的辅助函数（`fig2_behavioral_pls_20260916.py`）。
"""
import csv
import io
import os
import sys

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.legend_handler import HandlerTuple
from nilearn import datasets, plotting

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import figstyle_20260916 as fs
import fig2_behavioral_pls_20260916 as f3
from decimal import Decimal, ROUND_HALF_UP


def rnd3(x):
    """稿件口径：半进位到 3 位小数并带符号。

    2026-09-16 修正：原用 '%.3f'，对源值 0.3655 得 +0.365，而稿件按半进位为
    +0.366 → 图-稿末位不一致。改为 Decimal(ROUND_HALF_UP)（与
    `data_fig5_pathway_20260916.py` 的 rnd() 同口径）。
    """
    v = Decimal(str(float(x))).quantize(Decimal('0.001'), rounding=ROUND_HALF_UP)
    return ('+' if v >= 0 else '-') + str(abs(v))


fs.apply_style(7.0)
HERE = os.path.dirname(os.path.abspath(__file__))
FIG_W_MM, FIG_H_MM = 190.0, 158.0
FIG_W, FIG_H = FIG_W_MM * fs.MM, FIG_H_MM * fs.MM


def rect(x, y, w, h):
    return [x / FIG_W_MM, y / FIG_H_MM, w / FIG_W_MM, h / FIG_H_MM]


def rd(n):
    return list(csv.DictReader(io.open(os.path.join(HERE, n), encoding='utf-8-sig')))


es = rd('fig5_es_curves.csv')
vol = rd('fig5_volcano.csv')
spin = rd('fig5_spin.csv')
maps = rd('fig5_pathway_maps.csv')
reps = rd('fig5_representatives.csv')

ORDER = [r['pathway'] for r in reps]
# 2026-10-09：a 图例的 'Glycolysis' → 全名「Glycolysis / Gluconeogenesis」（与全稿 glcnm 统一）；
# 'OXPHOS' 为公认缩写、语义完整，保留。
DLEG = [{'Oxidative phosphorylation': 'OXPHOS',
         'Glycolysis / Gluconeogenesis': 'Glycolysis / Gluconeogenesis',
         'Ribosome': 'Ribosome', 'Proteasome': 'Proteasome',
         'MHC class I (GO:0002479)': 'MHC class I (GO:0002479)'}[p] for p in ORDER]
# d 面板行标签：MHC 用短名（全名会越过 c 面板右缘）；
# 「Glycolysis / Gluconeogenesis」折为两行（单行会压 a 面板 "10,000" 刻度——护栏⑦实测）。
# 2026-10-10 实测：折行块宽 16.6 mm（由 'Gluconeogenesis' 行决定），左置仍嫌长 ⇒ 全部标签
# 右缘由 125.0 右移至 127.9 mm（= MAP_X + 1.9）：左缘 111.3（距 '10,000' 右缘 110.1 约 1.2 mm），
# 右缘距脑图左缘（3D 轴内缩 2.9 mm，实测 128.9–129.3）约 1.0–1.4 mm。护栏⑦实装把关。
DLAB = [{'MHC class I (GO:0002479)': 'MHC class I',
         'Glycolysis / Gluconeogenesis': 'Glycolysis /\nGluconeogenesis'}.get(x, x) for x in DLEG]
# 5 条曲线用 figstyle 的 5 槽类别色板（2026-09-18：绿位 #1B7837 → #009E73，
# 原组合在蓝色盲下蓝/绿 ΔE 仅 8.3；新组合四种色觉最小 ΔE 全部 ≥15）
COL = {p: c for p, c in zip(ORDER, fs.CATEGORICAL5)}

fig = plt.figure(figsize=(FIG_W, FIG_H))

# =============================================================== a ES 曲线
ax = fig.add_axes(rect(13.0, 112.0, 99.0, 40.0))
rug_y = -0.10
for i, p in enumerate(ORDER):
    cur = [x for x in es if x['pathway'] == p]
    x = np.array([float(z['rank']) for z in cur])
    y = np.array([float(z['running_es']) for z in cur])
    o = np.argsort(x)
    ax.plot(x[o], y[o], color=COL[p], lw=1.0)
    hit = [float(z['rank']) for z in es if z['pathway'] == p + ' |HIT']
    ax.vlines(hit, rug_y - 0.022 * (4 - i) - 0.022, rug_y - 0.022 * (4 - i),
              color=COL[p], lw=0.35, alpha=0.8)
ax.axhline(0, color=fs.GREY, lw=0.6, ls=(0, (2.5, 1.6)))
ax.set_xscale('log')
ax.set_xlim(1, 17000)
ax.set_xticks([1, 10, 100, 1000, 10000])
ax.set_xticklabels(['1', '10', '100', '1,000', '10,000'])
ax.set_ylim(-0.34, 0.72)
ax.set_yticks([0.0, 0.2, 0.4, 0.6])
ax.set_xlabel('Gene rank in w1')
ax.set_ylabel('Running enrichment score')
ax.text(0.985, 0.03, 'line = running ES   ticks = gene-set hits',
        transform=ax.transAxes, fontsize=6.0, color=fs.GREY, ha='right', va='bottom')
ax.text(-0.005, 1.06, 'a', transform=ax.transAxes, fontsize=8, fontweight='bold',
        va='bottom', ha='right')
handles = [mpl.lines.Line2D([], [], color=COL[p], lw=1.2, label=DLEG[i])
           for i, p in enumerate(ORDER)]
ax.legend(handles=handles, loc='upper left', frameon=False, fontsize=6.0,
          handletextpad=0.5, borderpad=0.1, labelspacing=0.32)

# =============================================================== b 火山
axb = fig.add_axes(rect(14.0, 62.0, 98.0, 40.0))
nes = np.array([float(x['nes']) for x in vol])
fdr = np.array([float(x['fdr']) for x in vol])
sig = np.array([x['significant'] == 'True' for x in vol])
db = [x['db'] for x in vol]
lg = -np.log10(np.clip(fdr, 1e-4, 1))
for lab, mask, c, mk, ms in [('GO, q ≥ 0.05', np.array(db) == 'GO', fs.GREY, 'o', 1.0),
                             ('KEGG, q ≥ 0.05', np.array(db) == 'KEGG', fs.GREY, '^', 1.4)]:
    m = mask & (~sig)
    axb.plot(nes[m], lg[m], mk, ms=ms, color=c, alpha=0.45, mec='none', label=None)
for lab, mask, c, mk, ms in [('GO, q < 0.05', np.array(db) == 'GO', fs.CATEGORICAL[0], 'o', 1.6),
                             ('KEGG, q < 0.05', np.array(db) == 'KEGG', fs.ACCENT_RED, '^', 2.2)]:
    m = mask & sig
    axb.plot(nes[m], lg[m], mk, ms=ms, color=c, mec='none', label=lab)
axb.axhline(-np.log10(0.05), color=fs.BLACK, lw=0.6, ls=(0, (2.2, 1.6)))
axb.text(3.65, -np.log10(0.05) + 0.06, 'q = 0.05', fontsize=5.9, color=fs.BLACK,
         ha='right', va='bottom')
# 2026-09-18：**删去** 5 条代表通路在 b 面板的引线标注，改由 c 面板承担指认。
# 实测依据（`_trial_fig6bc_options_20260918.py`）：
#   ① OXPHOS / Ribosome / Proteasome 的 raw FDR q 浮点下溢为 0（全库 28/5,613 条），
#      经 clip(fdr, 1e-4, 1) 后三者 −log10 q 全 = 4.0000 → **y 轴完全重合**，
#      在 b 的坐标系里几何上不可分辨，换任何引线画法都无解；
#   ② 屏幕上 OXPHOS–Ribosome 仅相距 1.11 mm、Proteasome–MHC-I 仅 2.12 mm；
#   ③ 5 条引线总长 375.4 mm = b 轴宽（98 mm）的 3.8 倍，横穿数据区造成交叉。
#   b 的职责回归"全景 + 显著集中；代表通路的指认见 c。
# 2026-09-24（郝医生裁定 C15）：b 面板用「**加黑边的彩色大点**」标出 5 条代表通路。
# ⚠️ 不用引线 —— 2026-09-18 记录：OXPHOS/Ribosome/Proteasome 的 raw FDR q 下溢为 0，
#    经 clip(1e-4) 后 −log10 q 全 = 4.0000（y 轴重合），5 条引线总长 375.4 mm（轴宽 3.8 倍）
#    会横穿数据区 ⇒ 只画点、不画线，颜色与 a、c 一致，由 c 面板负责文字指认。
# ⚠️ 不能用 term 精确匹配：reps 的 pathway 是人工短名（如 'MHC class I (GO:0002479)'），
#    而 volcano 的 term 是 GO 全名（'antigen processing ... via MHC class I, TAP-dependent (GO:0002479)'）
#    ⇒ 实测只命中 4/5，MHC-I 那条**会静默丢掉**。改用 (db, NES) 匹配，并断言唯一命中。
for _i, _r in enumerate(reps):
    _c = [z for z in vol
          if z['db'] == _r['db'] and abs(float(z['nes']) - float(_r['nes'])) < 1e-6]
    assert len(_c) == 1, '代表通路 %r 在 volcano 里命中 %d 条' % (_r['pathway'], len(_c))
    axb.plot([float(_c[0]['nes'])],
             [-np.log10(np.clip(float(_c[0]['fdr']), 1e-4, 1))],
             'o', ms=3.4, mfc=COL[ORDER[_i]], mec=fs.BLACK, mew=0.5, zorder=6)
axb.set_xlabel('NES')
axb.set_ylabel('$-\\log_{10}$ FDR $q$')
axb.set_xlim(-3.9, 3.9)
axb.text(-0.005, 1.06, 'b', transform=axb.transAxes, fontsize=8, fontweight='bold',
         va='bottom', ha='right')
# 2026-09-19（郝医生指出）：原 `loc='lower left'` 的图例**压住了左下的灰点**
# —— 不显著点的 NES 范围 −2.235 ~ +1.810、lg 0 ~ 1.295，图例占的
# x ∈ [−3.9, −1.55]、y ∈ [−0.2, 0.9] 正是点最密的区域之一。
# 移到左上：该区 x ∈ [−3.9, −2.4]、y ≥ 3.0 完全空白（NES 最小 −2.235），
# 且与 a、c 两个面板的图例位置统一为 'upper left'。
# 2026-09-24（郝医生）：b 面板图例同样加第 3 条「representative (n = 5)」（5 色点组），
# 与 c 面板同法 ⇒ 图例条目数 == 图注类别数（3 类）。
_hb, _lb = axb.get_legend_handles_labels()
_hb.append(tuple(mpl.lines.Line2D([], [], marker='o', ms=2.0, color=c,
                                  mec=fs.BLACK, mew=0.4, linestyle='none')
                 for c in fs.CATEGORICAL5))
_lb.append('representative (n = 5)')
axb.legend(_hb, _lb, loc='upper left', frameon=False, fontsize=6.0, handletextpad=0.4,
           borderpad=0.1, labelspacing=0.3, markerscale=1.6,
           handler_map={tuple: HandlerTuple(ndivide=None)})

# =============================================================== c spin 存活
axc = fig.add_axes(rect(14.0, 12.0, 98.0, 40.0))
s_nes = np.array([float(x['nes']) for x in spin])
s_p = np.array([float(x['p_spin']) for x in spin])
s_ok = np.array([x['survived'] == 'True' for x in spin])
s_lg = -np.log10(np.clip(s_p, 1e-4, 1))
# 2026-09-19（郝医生裁定 S2）：c 面板背景改**灰阶**，把颜色全部让给 5 条代表通路。
# 原因：CATEGORICAL5 里 OXPHOS 的 #2166AC 与「存活」蓝**完全相同**、MHC-I 的 #B2182B
# 与「淘汰」红**完全相同** → 同面板内「蓝」「红」各两义（靠尺寸+黑边勉强区分）。
# 淘汰与否其实已由「是否落在 spin p = 0.05 虚线下方」表达，颜色本属冗余。
GREY_IN, GREY_OUT = '#6E6E6E', '#CBCBCB'      # 深灰 = 存活，浅灰 = 淘汰
axc.plot(s_nes[~s_ok], s_lg[~s_ok], 'o', ms=1.9, color=GREY_OUT,
         mec='none', label='eliminated (%d)' % int((~s_ok).sum()))
axc.plot(s_nes[s_ok], s_lg[s_ok], 'o', ms=1.9, color=GREY_IN,
         mec='none', alpha=0.9, label='retained (%d)' % int(s_ok.sum()))
axc.axhline(-np.log10(0.05), color=fs.BLACK, lw=0.6, ls=(0, (2.2, 1.6)))
LAB_X, LAB_HA = 2.93, 'right'
axc.text(LAB_X, -np.log10(0.05) + 0.06, '$p_{\\mathrm{spin}}$ = 0.05', fontsize=5.9,
         color=fs.BLACK, ha=LAB_HA, va='bottom')
# 自检：判读标签不得压住数据点。⚠️ 文本框必须**由 LAB_X/LAB_HA 推导**，不能写死坐标——
# 第一版把框写死在 2.93 处，反向验证（把标注挪回左端 1.52）时断言照样通过 = 假护栏。
# 近似：12 字符 ≈ 0.20 个 NES 单位、0.10 个 −log10 p 单位。
_w, _h = 0.20, 0.10
_lx0, _lx1 = (LAB_X - _w, LAB_X) if LAB_HA == 'right' else (LAB_X, LAB_X + _w)
_ly0 = -np.log10(0.05) + 0.06
_hit = int(((s_nes >= _lx0) & (s_nes <= _lx1) &
            (s_lg >= _ly0) & (s_lg <= _ly0 + _h)).sum())
assert _hit == 0, 'c 面板判读标签「spin p = 0.05」在 x∈[%.2f,%.2f] 处压住 %d 个数据点' % (_lx0, _lx1, _hit)
axc.set_xlabel('NES')
axc.set_ylabel('$-\\log_{10}$ $p_{\\mathrm{spin}}$')
# 2026-09-18：184 条显著通路的 NES 全在 +1.589 ~ +2.786（无一条 < 1），
# 原 xlim(-2.6, 3.6) 有 81% 轴宽为空、其中 x<1 段占 58% 轴宽且 0 个点 → 收紧到数据范围。
axc.set_xlim(1.45, 2.95)
axc.text(-0.005, 1.06, 'c', transform=axc.transAxes, fontsize=8, fontweight='bold',
         va='bottom', ha='right')
# 2026-09-18：5 条代表通路的指认**由 b 移入 c**。c 的横轴即 NES，
# 5 条在此的最小屏幕间距为 6.57 mm（b 面板仅 1.11 mm）→ 几何可分；
# 且 4/5 条 spin p < 0.05（糖酵解/糖异生 0.056 边缘），标出即「代表通路 4/5 通过表型层 spin 检验」的完整报告。
REP_C = [int(np.argmin(np.abs(s_nes - float(r_['nes'])))) for r_ in reps]
# 自检：标注的 5 个点必须真的对应 5 条代表通路（既不能用模糊匹配，也不能落到别的点上）
for k, i in enumerate(REP_C):
    assert abs(s_nes[i] - float(reps[k]['nes'])) < 1e-6, \
        'c 面板第 %d 条代表通路 NES 匹配失败（%.6f vs %.6f）' \
        % (k, s_nes[i], float(reps[k]['nes']))
assert len(set(REP_C)) == 5, 'c 面板 5 条代表通路落到了重复的点上：%s' % REP_C
# 自检：代表通路 4/5 存活（糖酵解/糖异生 0.056 边缘为预期；图注为「4 of the 5」完整报告）
_n_ok5 = sum(1 for i in REP_C if s_p[i] < 0.05)
assert _n_ok5 == 4, \
    'c 面板代表通路存活数 %d（应为 4/5；糖酵解 0.056 边缘）' % _n_ok5
_excl5 = [ORDER[k] for k, i in enumerate(REP_C) if not (s_p[i] < 0.05)]
assert _excl5 == ['Glycolysis / Gluconeogenesis'], \
    '唯一未存活者应为糖酵解/糖异生：%s' % _excl5
for k, i in enumerate(REP_C):
    axc.plot(s_nes[i], s_lg[i], 'o', ms=5.2, color=COL[ORDER[k]],
             mec=fs.BLACK, mew=0.65, zorder=6)
# 图例第 3 条用 5 色点组（HandlerTuple），使**图例条目数 == 图注类别数**（3 类）；
# 描边必须与图上实际点一致（黑边），否则图例与图元不符
_h, _l = axc.get_legend_handles_labels()
_h.append(tuple(mpl.lines.Line2D([], [], marker='o', ms=3.0, color=c,
                                 mec=fs.BLACK, mew=0.4, linestyle='none')
                for c in fs.CATEGORICAL5))
_l.append('representative (n = 5)')
# 2026-10-10：从 upper right **复原为 upper left**（与 a/b 面板统一布局——郝医生 09-19 裁定；
# 10-09 外迁系误判）。复核实测：最高点（secondary alcohol metabolic process, 1.861/3.000）
# 仅落在透明外框的**空白角**（框角内 0.4 px），距最近图例文字 6.4 mm、全 184 点与图例内容
# 最小净距 0.78 mm（ryanodine 点 vs 'representative' 行）⇒ 无任何视觉遮挡。
# 护栏⑤已同步精化为「内容级」判据（点标记 vs 图例文字/手柄，净距 >= 0.5 mm）。
axc.legend(_h, _l, loc='upper left', frameon=False, fontsize=6.0,
           handletextpad=0.4, borderpad=0.1, labelspacing=0.3,
           handler_map={tuple: HandlerTuple(ndivide=None)})

# =============================================================== d 表达图
vmap = f3.annot_vertex_map()
fsa = datasets.fetch_surf_fsaverage('fsaverage5')
BG = {}
vmax = max(abs(float(x['median_z'])) for x in maps)
vmax = float(np.ceil(vmax * 10) / 10)
# 2026-10-09：无数据区（RH_SomMot_2）用黑色标出——与 Fig.4c 及图注「黑色区」一致
# （机制照抄 Fig.4 真生成器：cmap 副本 + set_under('black') + 该区顶点赋 -vmax*10）。
CMAP_F5D = fs.CMAP_DIV.copy()
CMAP_F5D.set_under('black')
import nibabel as nb_
N_VERT = {'left': nb_.load(fsa['sulc_left']).darrays[0].data.shape[0],
          'right': nb_.load(fsa['sulc_right']).darrays[0].data.shape[0]}
for hemi, fn in (('left', 'lh'), ('right', 'rh')):
    BG[hemi] = f3.clean_bgmap(fsa['sulc_%s' % hemi],
                              f3.ANNOT + '/%s.Schaefer2018_100Parcels_7Networks_order.annot' % fn)

names = sorted({x['region'] for x in maps})
MAP_X, MAP_W, ROW_P, ROW_H = 126.0, 26.0, 28.0, 25.0   # 右移，给行标签让位
D_Y0 = 12.0                      # 最底一行地图的下沿
_DLAB_OBJS = []                  # 护栏⑥用：收集 d 行标签对象
for i, p in enumerate(ORDER):
    y0 = D_Y0 + (len(ORDER) - 1 - i) * ROW_P
    for j, hemi in enumerate(('left', 'right')):
        sub = {x['region']: float(x['median_z']) for x in maps if x['pathway'] == p}
        hemi_names = [n for n in names if (n.startswith('LH') if hemi == 'left'
                                           else n.startswith('RH'))]
        vals = [sub.get(n, np.nan) for n in hemi_names]
        vmap_h = {k: v for k, v in vmap.items() if
                  (k.startswith('LH') if hemi == 'left' else k.startswith('RH'))}
        nv = N_VERT[hemi]
        v = f3.vertex_map(hemi_names, vals, vmap_h, nv)
        # 只把**无表达数据的那个区**（RH_SomMot_2）标黑；其余非 parcel 顶点保持背景，不得涂黑
        if 'RH_SomMot_2' in vmap_h:
            v[vmap_h['RH_SomMot_2']] = -vmax * 10.0     # 配合 set_under('black')
        v[~np.isfinite(v)] = 0.0                        # 皮层外/内侧壁等非 parcel 顶点 → 中性
        axm = fig.add_axes(rect(MAP_X + j * (MAP_W + 1.2), y0, MAP_W, ROW_H),
                            projection='3d')     # nilearn 表面图需要 3D 轴
        plotting.plot_surf_stat_map(
            fsa['infl_%s' % hemi], v, hemi=hemi, view='lateral',
            bg_map=BG[hemi], cmap=CMAP_F5D, vmin=-vmax, vmax=vmax,
            threshold=0.0, symmetric_cbar=True, colorbar=False,
            axes=axm, figure=fig)
        for _coll in axm.collections:      # 脑表面栅格化（否则矢量文件 ~38 MB）
            _coll.set_rasterized(True)
        axm.text2D(0.03, 0.05, 'L' if hemi == 'left' else 'R',
                   transform=axm.transAxes, fontsize=6.0, color=fs.GREY,
                   ha='left', va='bottom')
    _label_y = y0 + ROW_H / 2 - (0.8 if '\n' in DLAB[i] else 0.0)   # 折行块（3 行）下移 0.8 mm 避 a 轴底线
    _DLAB_OBJS.append(fig.text(
        (MAP_X + 1.9) / FIG_W_MM, _label_y / FIG_H_MM,       # 右缘 127.9 mm（见 DLAB 注释）
        '%s\nr = %s' % (DLAB[i], rnd3([r_ for r_ in reps
                                       if r_['pathway'] == p][0]['r'])),
        fontsize=6.1, color=fs.BLACK, ha='right', va='center', linespacing=1.35))

D_TOP = D_Y0 + (len(ORDER) - 1) * ROW_P + ROW_H
D_BOT = D_Y0
fig.text((MAP_X - 2.0) / FIG_W_MM, (D_TOP + 2.0) / FIG_H_MM, 'd',
         fontsize=8, fontweight='bold', ha='left', va='bottom')
cax = fig.add_axes([(MAP_X + 2 * (MAP_W + 1.2) - 2.0) / FIG_W_MM,
                    (D_BOT + 4.0) / FIG_H_MM, 2.6 / FIG_W_MM,
                    (D_TOP - D_BOT - 8.0) / FIG_H_MM])
sm = mpl.cm.ScalarMappable(norm=mpl.colors.Normalize(-vmax, vmax), cmap=fs.CMAP_DIV)
cb = fig.colorbar(sm, cax=cax, orientation='vertical')
cb.set_label('Median $Z$', fontsize=6.4, labelpad=2.0)
cb.ax.tick_params(labelsize=6.0, width=0.5, length=1.6)
cb.outline.set_linewidth(0.5)

# =============================================================== 自检
fig.canvas.draw()
rend = fig.canvas.get_renderer()
fb = fig.get_window_extent()
qa = []
for a in fig.axes:
    if not a.get_visible():
        continue
    for t in list(a.texts) + [a.xaxis.label, a.yaxis.label]:
        if not t.get_visible() or not t.get_text().strip():
            continue
        bb = t.get_tightbbox(rend)
        if bb is None:
            continue
        out = max(fb.x0 - bb.x0, fb.y0 - bb.y0, bb.x1 - fb.x1, bb.y1 - fb.y1)
        if out > -0.5:
            qa.append((t.get_text()[:32].replace('\n', '|'), round(out, 1)))
over = [] if not qa else qa
print('  越界文字：%d 个' % len(over))
for o in over:
    print('     %-34s 越界 %.1f px' % o)

print('  自检：')
print('    a 曲线 %d 条' % len(ORDER))
print('    b 火山 %d 条（显著 %d）' % (len(vol), int(sig.sum())))
print('    c spin %d 条（存活 %d / 淘汰 %d）'
      % (len(spin), int(s_ok.sum()), int((~s_ok).sum())))
print('    d 表达图 %d 条通路 × 2 半球（共享 vmax ±%.1f）' % (len(ORDER), vmax))
assert len(vol) == 5613 and int(sig.sum()) == 184
assert len(spin) == 184 and int(s_ok.sum()) == 127
print('    ✓ 条数与正文一致（GO+KEGG 5,613；显著 184；spin 存活 127 / 淘汰 57）')

# ---- 2026-09-18 新增护栏：代表通路的指认必须在 c，且 5 点几何可分 ----
# ① 反向护栏：b 面板若被加回引线，会出现 Annotation 对象（annotate 的产物）→ 立即失败。
#    （注意：不能只数 axb.texts——面板字母 'b' 与 'q = 0.05' 都是普通 Text）
_nann = sum(isinstance(t, mpl.text.Annotation) for t in axb.texts)
assert _nann == 0, 'b 面板仍有 %d 条引线标注（annotate）——应已全部删除' % _nann
# ② 正向护栏：c 面板 5 点最小屏幕间距必须 > 4 mm，否则标注相互重叠
_px = [axc.transData.transform((s_nes[i], s_lg[i])) for i in REP_C]
_d_mm = min(np.hypot(_px[a][0] - _px[b][0], _px[a][1] - _px[b][1])
            for a in range(5) for b in range(a + 1, 5)) / fig.dpi * 25.4
print('    b 面板引线 %d 条（应为 0）；c 面板 5 条代表通路最小间距 %.2f mm'
      % (_nann, _d_mm))
assert _d_mm > 4.0, 'c 面板 5 条代表通路过近（最小 %.2f mm），标注会相互重叠' % _d_mm
print('    ✓ 代表通路 4/5 存活（未存活者：糖酵解/糖异生 0.056 边缘）；p = %s'
      % ', '.join('%.3f' % float(s_p[i]) for i in REP_C))
# ③ 图例条目数必须等于图注类别数（c 面板 3 类：淘汰 / 存活 / 代表通路）
_lg_c = [t.get_text() for t in axc.get_legend().get_texts()]
assert len(_lg_c) == 3, 'c 面板图例 %d 条（应为 3）：%s' % (len(_lg_c), _lg_c)
print('    ✓ c 面板图例 3 条：%s' % ' / '.join(_lg_c))
# ④ S2 护栏：c 面板前景（5 条代表通路）与背景（存活/淘汰灰阶）不得撞色
def _hx(c):
    return mpl.colors.to_hex(c).lower()


def _rgb(h):
    h = h.lstrip('#')
    return tuple(int(h[i:i + 2], 16) for i in (0, 2, 4))


_fg5 = {_hx(c) for c in fs.CATEGORICAL5}
_bg = {_hx(GREY_IN), _hx(GREY_OUT)}
_clash = _fg5 & _bg
assert not _clash, 'c 面板背景色与代表通路色撞色（同色两义）：%s' % _clash
_dc = sum((a - b) ** 2 for a, b in zip(_rgb(_hx(GREY_IN)), _rgb(_hx(GREY_OUT)))) ** .5
_draw = {_hx(ln.get_color()) for ln in axc.lines}
assert _bg <= _draw, 'c 面板图上未用到声明的灰阶色：%s' % (_bg - _draw)
print('    ✓ c 面板 S2 生效：背景 %s(存活)/%s(淘汰)，RGB 距离 %.1f；与 5 色前景无撞色'
      % (_hx(GREY_IN), _hx(GREY_OUT), _dc))
assert _dc > 60, 'c 面板灰阶两色过近（RGB 距离 %.1f），1.9 pt 小点上不可分辨' % _dc

# ⑤ 图例不得压住数据点（2026-09-19 建；2026-10-10 精化为「内容级」判据）
#    缘由：c 面板最高点（secondary alcohol…, 1.861/3.000）仅落在左上图例**透明外框的空白角**
#    （框角内 0.4 px，距最近文字 6.4 mm、全 184 点与内容最小净距 0.78 mm）——无视觉遮挡；
#    原「点中心在 bbox 内」判据据此误报，曾把 c 图例迁往 upper right（破坏 a/b/c 统一左上布局）。
#    新判据：数据点（扣各自标记半径）与**图例文字/手柄内容框**的净距 >= 0.5 mm；bbox 内点数仅打印。
def _leg_net_clear_mm(_leg, _xy, _rad_px, _rend):
    """点到图例内容（文字 bbox 左扩 35 px 含手柄）的最小净距（mm，扣标记半径）。"""
    _cb = [(t.get_window_extent(_rend).x0 - 35.0, t.get_window_extent(_rend).x1,
            t.get_window_extent(_rend).y0, t.get_window_extent(_rend).y1)
           for t in _leg.get_texts()]
    _best = 1e9
    for (_px_, _py_), _r_ in zip(_xy, _rad_px):
        for (_x0, _x1, _y0, _y1) in _cb:
            _dx = max(_x0 - _px_, _px_ - _x1, 0.0)
            _dy = max(_y0 - _py_, _py_ - _y1, 0.0)
            _best = min(_best, float(np.hypot(_dx, _dy)) - _r_)
    return _best * 25.4 / fig.dpi


_a_pts = np.array([(float(z['rank']), float(z['running_es']))
                   for z in es if '|HIT' not in z['pathway']])
_pt2px = fig.dpi / 72.0                       # pt → 显示像素（标记 ms 为直径 pt）
_rad_a = np.zeros(len(_a_pts))                # a 为曲线（无标记）
_rad_b = np.array([((1.6 if r['db'] == 'GO' else 2.2) if r['significant'] == 'True'
                    else (1.0 if r['db'] == 'GO' else 1.4)) / 2.0 * _pt2px for r in vol])
_rad_c = np.full(len(s_nes), 1.9 / 2.0 * _pt2px)
_rad_c[REP_C] = 5.2 / 2.0 * _pt2px
for _ax_, _nm, _pts_, _rad in [(ax, 'a', _a_pts, _rad_a),
                               (axb, 'b', np.column_stack([nes, lg]), _rad_b),
                               (axc, 'c', np.column_stack([s_nes, s_lg]), _rad_c)]:
    _leg = _ax_.get_legend()
    assert _leg is not None, '%s 面板无图例' % _nm
    _xy = _ax_.transData.transform(_pts_)
    _bb = _leg.get_window_extent(rend)
    _in = int(((_xy[:, 0] >= _bb.x0) & (_xy[:, 0] <= _bb.x1) &
               (_xy[:, 1] >= _bb.y0) & (_xy[:, 1] <= _bb.y1)).sum())
    _clr = _leg_net_clear_mm(_leg, _xy, _rad, rend)
    print('    %s 图例：bbox 内 %d 点（参考）；数据点—图例内容最小净距 %.2f mm'
          % (_nm, _in, _clr))
    assert _clr >= 0.5, \
        '%s 面板数据点距图例内容净距 %.2f mm < 0.5（图例压点！）' % (_nm, _clr)
print('    ✓ 三个面板：数据点与图例内容（文字/手柄）净距均 >= 0.5 mm')

# ⑥ d 行标签不得压住 b/c 面板（2026-10-09 新增；Gly 全名扩宽后实测）
_panel_boxes = {'b': axb.get_window_extent(rend), 'c': axc.get_window_extent(rend)}
for _t in _DLAB_OBJS:
    _bb = _t.get_window_extent(rend)
    for _pn, _pb in _panel_boxes.items():
        _ov = (_bb.x0 < _pb.x1) and (_bb.x1 > _pb.x0) and (_bb.y0 < _pb.y1) and (_bb.y1 > _pb.y0)
        assert not _ov, 'd 行标签 %r 压住 %s 面板' % (_t.get_text().split('\n')[0], _pn)
print('    ✓ d 行标签 5 条均未压住 b/c 面板')

# ⑦ d 行标签不得与 a 面板的刻度/轴文字相碰（2026-10-09 建；2026-10-10 实装 assert）
#    背景：Gly 全名折行后块宽实测 16.6 mm，左置仍压 '10,000'（右缘 110.1）；
#    修复：标签右缘 125.0 → 127.9 mm（= MAP_X + 1.9，见 DLAB 注释）。本护栏把实测几何固化。
_a_texts = ([t for t in ax.get_xticklabels() if t.get_text().strip()]
            + [ax.xaxis.label] + list(ax.texts))
_hit7 = []
for _t in _DLAB_OBJS:
    _bb = _t.get_window_extent(rend)
    for _ta in _a_texts:
        if not _ta.get_visible() or not _ta.get_text().strip():
            continue
        _ab = _ta.get_window_extent(rend)
        if (_bb.x0 < _ab.x1) and (_bb.x1 > _ab.x0) and (_bb.y0 < _ab.y1) and (_bb.y1 > _ab.y0):
            _hit7.append((_t.get_text().split('\n')[0][:24], _ta.get_text()[:16]))
assert not _hit7, 'd 行标签压住 a 面板文字：%s' % _hit7
# 同时固化「标签右缘 — 脑图左缘」间距（脑图左缘 = 轴左缘 + 2.9 mm 内缩，2026-10-10 实测）
_pxmm = fig.get_window_extent(rend).width / FIG_W_MM
_gap_mm = (MAP_X + 2.9) - max(_t.get_window_extent(rend).x1 / _pxmm for _t in _DLAB_OBJS)
assert _gap_mm >= 0.5, 'd 行标签右缘距脑图左缘仅 %.2f mm（应 >= 0.5）' % _gap_mm
print('    ✓ d 行标签 5 条：未碰 a 面板文字，右缘距脑图左缘 %.2f mm' % _gap_mm)

mpl.rcParams['savefig.bbox'] = None
fs.save_all(fig, 'Fig5_pathway_enrichment', outdir=HERE,
            formats=('pdf', 'svg', 'eps', 'png'), bbox_inches=None)
from PIL import Image
w_px, h_px = Image.open(os.path.join(HERE, 'Fig5_pathway_enrichment.png')).size
print('  成品像素 %d × %d → %.1f × %.1f mm @300dpi（目标 %.1f × %.1f）'
      % (w_px, h_px, w_px / 300 * 25.4, h_px / 300 * 25.4, FIG_W_MM, FIG_H_MM))
