# -*- coding: utf-8 -*-
"""
fig4_celltype_20260916.py —— Fig. 4 疲劳相关基因的细胞类型富集

面板（2026-09-18 改版：删除原 c 面板「20 种子热图」，原 d 升为 c）
  a 7 类细胞基因集在 w1 排序中的中位秩（点 + 10,000 次基因置换空分布 95% CI；
    竖线为全基因组中位秩基线 = 7,817）
  b 7 类表达图与 LV1 loading 的空间相关 r（棒棒糖；实心 = spin p < 0.05；
    数字前缀 p = 以区别于横轴的 r）
  c HPA 脑表达过滤前后（空心 = 过滤前，实心 = 过滤后；颜色 = 过滤后是否 BH 显著）

版面（2026-09-18 第二轮）：a / b / c **三张并列等宽**（50 / 45 / 45 mm，190 × 62 mm），
  三图共用同一 ylim (-0.9, 6.9) 故行完全对齐；**行标签只画在 a**（b / c 行序相同），
  三图的字母位置统一在 (0.0, 1.04)。a 面板的 ← front-enriched / back-enriched → 已删除
  （图注已说明基线两侧的含义）。

为什么删原 c 面板（2026-09-18 实测）
  7 类中 4 类（小胶质/少突/兴奋/抑制）在 20 个种子下的 p 全部钉在置换下限 1e-4
  （−log10 p = 4.00，恰为色标上限）→ 行内零变化、20 列不携带信息；且「仅 OPC 0/20、
  其余 20/20」正文 §3.4 已完整写出 → 双重冗余。若日后需要重新展示种子稳健性，
  宜改为「中位 [范围]」范围条，而不是 7×20 饱和热图。

图型依据（scipilot chart_selection）
  * 分类 × 点估计 + 区间 → 点 + 误差棒（**不用均值柱**：柱状看不出是否偏离置换零分布）
  * 分类 × 有正负的连续量 → 棒棒糖，0 参考线
  * 同一批对象的前后对比 → 哑铃图（不是分组柱状）

数据（脚本读取，无硬编码）
  PSF_celltype_enrich_seidlitz_Dataset5_newmasks_20260915.csv
  PSF_celltype_enrich_seidlitz_Dataset5_hpafilter_newmasks_20260915.csv
  permutation_newmasks_20260915/PSF_sal_specificity_spin_Dataset5_result_newmasks_20260915.csv
"""
import csv
import io
import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import matplotlib as mpl
import matplotlib.pyplot as plt

import figstyle_20260916 as fs

HERE = os.path.dirname(os.path.abspath(__file__))
BASE = 'F:/dcmtonii/run62_crosssectional'
AHBA = BASE + '/PSF_ahba_100'
REV = BASE + '/PSF_pls_behavior/断连稿写作_20260907/修订_revision_20260910'

fs.apply_style(base_size=7.0)
FIG_W, FIG_H = 190.0, 68.0
GENOME_MEDIAN_RANK = 7817.0        # (15,633 + 1) / 2，全基因组中位秩基线

NAME = {
    'Oligodendrocytes': 'Oligodendrocytes',
    'Micro': 'Microglia',
    'Astro': 'Astrocytes',
    'Endo': 'Endothelial',
    'OPC': 'OPC',
    'Excitatory_Neurons': 'Excitatory neurons',
    'Inhibitory_Neurons': 'Inhibitory neurons',
}


def rect(x, y, w, h):
    return [x / FIG_W, y / FIG_H, w / FIG_W, h / FIG_H]


def rd(p):
    return list(csv.DictReader(io.open(p, encoding='utf-8-sig')))


def ci95(s):
    return [float(v) for v in s.strip('[]').split(',')]


def load():
    enr = rd(AHBA + '/PSF_celltype_enrich_seidlitz_Dataset5_newmasks_20260915.csv')
    hpa = rd(AHBA + '/PSF_celltype_enrich_seidlitz_Dataset5_hpafilter_newmasks_20260915.csv')
    spa = rd(AHBA + '/permutation_newmasks_20260915/'
             'PSF_sal_specificity_spin_Dataset5_result_newmasks_20260915.csv')

    m = {}
    for r in enr:
        m[r['cell_type']] = dict(
            n=int(r['n_marker']), rank=float(r['median_rank']),
            ci=ci95(r['null_95CI']), p=float(r['p_prand_twotail']),
            bh=float(r['fdr_bh']))
    for r in hpa:
        m[r['cell_type']]['bh_hpa'] = float(r['fdr_bh'])
    # 2026-09-19：**负荷侧与疲劳侧同源** —— 同一个 spin 结果文件里就有
    # mode='discon_mean'（断连负荷图）的 7 行，同一套球面旋转零模型
    # （1,000 次、seed = 0、双侧）。所以画进同一面板是同质可比，
    # 不需要新数据、不需要新零模型。
    for r in spa:
        if not r['feature'].startswith('Cell_'):
            continue
        key = r['feature'][len('Cell_'):]
        if key not in m:
            continue
        if r['mode'] == 'salience':
            m[key]['r'] = float(r['r'])
            m[key]['spin_p'] = float(r['P_SPIN'])
        elif r['mode'] == 'discon_mean':
            m[key]['r_bur'] = float(r['r'])
            m[key]['spin_p_bur'] = float(r['P_SPIN'])
    order = sorted(m, key=lambda k: m[k]['rank'])
    # 护栏①：7 类必须两套值齐全（防止上游文件改了 mode 取值而静默丢一半数据）
    _miss = [k for k in order if not {'r', 'spin_p', 'r_bur', 'spin_p_bur'} <= set(m[k])]
    assert not _miss, '以下细胞类型缺疲劳侧或负荷侧数据：%s' % _miss
    return m, order


def panel_a(fig, m, order):
    ax = fig.add_axes(rect(33.0, 18.0, 50.0, 42.0))
    ys = np.arange(len(order))[::-1]
    for y, k in zip(ys, order):
        d = m[k]
        col = fs.CATEGORICAL[0] if d['bh'] < 0.05 else fs.GREY
        ax.plot(d['ci'], [y, y], color=col, lw=0.8, alpha=0.55, zorder=2)
        ax.plot([d['ci'][0], d['ci'][1]], [y, y], color=col, lw=0, marker='|',
                ms=3.0, alpha=0.9, zorder=2)
        ax.plot([d['rank']], [y], 'o', ms=3.8, color=col, zorder=3)
    ax.axvline(GENOME_MEDIAN_RANK, color=fs.BLACK, lw=0.6, ls=(0, (2.4, 1.6)), zorder=1)
    # 2026-09-18：标注改到坐标区上方单行（原在区内 y=6.45 向下两行，压住最上行少突胶质细胞
    # 的误差棒，其空分布 95% CI = [7276, 8359]），并补上基线数值。
    ax.text(GENOME_MEDIAN_RANK + 70, len(order) - 0.02 + 0.30,
            'genome-wide median rank\n(7,817)', fontsize=5.9, color=fs.BLACK,
            va='bottom', ha='left', linespacing=1.25, clip_on=False)
    ax.set_yticks(ys)
    # 三张图并列：仅 a 显示行标签，b / c 共用（行序完全相同，见脚本头与图注）
    # 2026-09-19：**恢复基因集大小**（郝医生裁定 Table 2 删除 → 数值必须有新家）。
    # 原（09-18）为「不与 Table 2 重复」而删掉，结果 Fig. 5 图注反而写成
    # 「各类基因集大小见 Table 2」——图反向依赖表。加回后依赖解除。
    ax.set_yticklabels(['%s (%d)' % (NAME[k], m[k]['n']) for k in order], fontsize=6.2)
    ax.set_ylim(-0.9, len(order) - 0.1)
    ax.set_xlim(4200, 10600)
    ax.set_xticks([5000, 7000, 9000])          # 面板宽 50 mm，6 个刻度会挤在一起
    ax.set_xlabel('Median rank in the w1 gene ranking', labelpad=1.5)
    ax.xaxis.set_label_coords(0.5, -0.135)
    ax.tick_params(length=1.8, pad=1.2)
    # 2026-09-18：方向词（← front-enriched / back-enriched →）已按郝医生意见删除
    # ——  图注已说明「低于基线为排序前端富集、高于基线为排序后端」，图上不必重复。
    ax.text(0.0, 1.04, 'a', transform=ax.transAxes, ha='left', va='bottom',
            fontsize=8, fontweight='bold', color=fs.BLACK)
    print('  a 中位秩（按秩升序）：%s'
          % ', '.join('%s %.0f' % (k, m[k]['rank']) for k in order))


def panel_b(fig, m, order):
    """空间层：7 类表达图与 **LV1 loading（疲劳模式）** 及 **断连负荷图** 的相关。

    2026-09-19 郝医生裁定：Table 2 与 Fig. 5 数据高度重复（细胞类型、方向、中位秩、
    BH p、空间 r/p 在图上都有）→ Table 2 降为补充表，**断连负荷对照并入本面板**。
    画法：同轴双系列。圆点 = 与 LV1 loading（颜色仍只编码**疲劳侧** spin p < 0.05，
    不引入第二套颜色语义）；空心方块 = 与断连负荷图；点线连接同一细胞类型的两个值，
    使「解离」（小胶质 +0.290 vs −0.212）一眼可见。
    ⚠️ 连线用**点线**而非实线：a 面板的误差棒也是水平线，实线会造成误读。
    """
    ax = fig.add_axes(rect(90.0, 18.0, 45.0, 42.0))
    ys = np.arange(len(order))[::-1]
    _r_labels = []                             # 收集 p 标签，供 main() 做实测越界检查
    for y, k in zip(ys, order):
        d = m[k]
        # ⚠️ 2026-09-19：**每个符号用颜色编码"它自己的"显著性**，不要让方块跟随行颜色。
        # 若方块跟随行色：少突负荷侧 p = 0.629（不显著）却会画成蓝方块 → 误导。
        # 全图统一语义仍是「蓝 = 该值达 p < 0.05、灰 = 未达」，与 a–c 一致。
        sig_f = d['spin_p'] < 0.05
        sig_b = d['spin_p_bur'] < 0.05
        col_f = fs.CATEGORICAL[0] if sig_f else fs.GREY
        col_b = fs.CATEGORICAL[0] if sig_b else fs.GREY
        # 配对点线（下层）
        ax.plot([d['r_bur'], d['r']], [y, y], color=fs.GREY, lw=0.6,
                ls=(0, (1.6, 1.6)), alpha=0.55, zorder=1)
        # 负荷侧：空心方块（空心与疲劳侧的实心/空心解耦，靠形状区分系列）
        ax.plot([d['r_bur']], [y], marker='s', ms=2.7, mfc='white',
                mec=col_b, mew=0.7, zorder=3)
        # 疲劳侧：圆点（实心 = spin p < 0.05），保持 09-18 以来的编码
        ax.plot([d['r']], [y], marker='o', ms=3.8, color=col_f,
                mfc=col_f if sig_f else 'white', mec=col_f, zorder=4)
        _lt = ax.text(d['r'] + (0.024 if d['r'] >= 0 else -0.024), y,
                      'p = %.3f' % d['spin_p'], fontsize=5.8, color=col_f,
                      ha='left' if d['r'] >= 0 else 'right', va='center')
        _r_labels.append(_lt)
    # 参考线一律虚线：a 的 7,817 与 c 的 0.05 都是 (0,(2.4,1.6)) 黑虚线；
    # 本条(2026-09-18 前)是全稿仅有的两条黑实线参考线之一（另一条在已随 2026-09-22 修订删除的旧 Fig. 7a），已统一。
    ax.axvline(0, color=fs.BLACK, lw=0.6, ls=(0, (2.4, 1.6)), zorder=1)
    ax.set_yticks(ys)
    ax.set_yticklabels([''] * len(order))      # 行标签见 a（顺序相同）
    ax.set_ylim(-0.9, len(order) - 0.1)        # 与 a / c 同 ylim → 三图行对齐
    ax.set_xlim(-0.40, 0.58)                   # 容得下 "p = 0.005"；负荷最小值 −0.212 亦在内
    # 2026-09-19：轴标签须同时涵盖两个系列（符号分工由全图图注说明，b 不另设图例
    # —— 本图惯例是「符号/颜色由 a–c 共用图注统一说明」，且 45 mm 宽放不下
    # 2 条图例而不压住数据：左下角压抑制性/兴奋性行，左上角压少突行）。
    ax.set_xlabel('$r$ (expression map vs LV1 loading / burden)', labelpad=1.5)
    ax.xaxis.set_label_coords(0.5, -0.135)
    ax.tick_params(length=1.8, pad=1.2, labelsize=6.2)
    ax.text(0.0, 1.04, 'b', transform=ax.transAxes, ha='left', va='bottom',
            fontsize=8, fontweight='bold', color=fs.BLACK)
    print('  b 疲劳侧 r / spin p：%s'
          % ', '.join('%s %+.4f/%.4f' % (k, m[k]['r'], m[k]['spin_p']) for k in order))
    print('  b 负荷侧 r / spin p：%s'
          % ', '.join('%s %+.4f/%.4f' % (k, m[k]['r_bur'], m[k]['spin_p_bur']) for k in order))
    # ---- 2026-09-19 新增护栏 ----
    # ② 两个点都必须落在 xlim 内（负荷点超出轴会被静默裁掉，图上只剩疲劳点）
    _xl = ax.get_xlim()
    _out = [k for k in order
            if not (_xl[0] + 0.02 < m[k]['r_bur'] < _xl[1] - 0.02
                    and _xl[0] + 0.02 < m[k]['r'] < _xl[1] - 0.02)]
    assert not _out, 'b 面板有值落在 xlim 外（会被裁掉）：%s' % _out
    # ③ 标签越界：**交给 main() 用实测渲染宽度检查**（本函数内 fig 尚未 draw，
    #    测不了文字 bbox）。2026-09-19：原实现用硬编码字符宽度估算，且容限写成
    #    max(|xlim[0]|, |xlim[1]|) —— 反向验证证明它**永远通过**（假护栏），已改。
    ax._r_labels = _r_labels
    # ④ 两侧值必须真的来自不同 mode（防上游改列名后两套值变成同一套）
    _same = [k for k in order if abs(m[k]['r'] - m[k]['r_bur']) < 1e-9]
    assert not _same, 'b 面板疲劳侧与负荷侧数值完全相同（疑似读错 mode）：%s' % _same
    print('     ✓ 负荷点全部在轴内；两套值不同源')


def panel_c(fig, m, order):
    """HPA 脑表达过滤前后（原 panel_d，2026-09-18 升为 c）

    编码（原图与图注不符，已统一）：
      空心圆 = 过滤前（全部基因）；实心圆 = 过滤后
      颜色   = **过滤后**是否 BH 显著（蓝 = p < 0.05，灰 = 否）
    注意：小胶质/少突/兴奋/抑制四类过滤前后 p 完全相同（均处于置换下限 1e-4），
    两点重合 → 图上呈现为「实心点嵌在空心环内」，即"未变化"。
    """
    ax = fig.add_axes(rect(142.0, 18.0, 45.0, 42.0))
    ys = np.arange(len(order))[::-1]
    for y, k in zip(ys, order):
        a, b = m[k]['bh'], m[k]['bh_hpa']
        col = fs.CATEGORICAL[0] if b < 0.05 else fs.GREY
        ax.plot([a, b], [y, y], color=col, lw=0.7, zorder=2, alpha=0.85)
        # 「过滤前」的空心环须**略大于**「过滤后」的实心点：四类过滤前后 p 完全相同、
        # 两点坐标重合，同尺寸会被实心点完全盖住（实测 4× 放大后仍看不见环）。
        ax.plot([a], [y], marker='o', ms=4.6, mfc='white', mec=fs.GREY,
                mew=0.8, zorder=3)
        ax.plot([b], [y], marker='o', ms=3.0, color=col, zorder=4)
    ax.axvline(0.05, color=fs.BLACK, lw=0.6, ls=(0, (2.4, 1.6)), zorder=1)
    # 面板仅 45 mm 宽：标注贴线左侧（行数据最左在 2e-4，仍在标注左界之外），避开右上图例
    ax.text(0.042, len(order) - 0.35, 'BH $p$ = 0.05', fontsize=5.7, color=fs.BLACK,
            va='top', ha='right')
    ax.set_xscale('log')
    ax.set_xlim(5e-5, 4.0)
    ax.set_xticks([1e-4, 1e-3, 1e-2, 1e-1, 1])
    ax.set_xticklabels(['$10^{-4}$', '$10^{-3}$', '$10^{-2}$', '$10^{-1}$', '1'])
    ax.xaxis.set_minor_locator(mpl.ticker.NullLocator())
    ax.set_yticks(ys)
    ax.set_yticklabels([''] * len(order))      # 行标签见 a（顺序相同）
    ax.set_ylim(-0.9, len(order) - 0.1)        # 与 a / b 同 ylim → 三图行对齐
    ax.set_xlabel('BH-adjusted $p$ (gene level)', labelpad=1.5)
    ax.xaxis.set_label_coords(0.5, -0.135)
    ax.tick_params(length=1.8, pad=1.2, labelsize=6.2)
    # 图例（2026-09-18 定稿：**3 类**，与图注「空心=过滤前、实心=过滤后、
    # 颜色=过滤后是否 BH 显著」一一对应；图注是三类，图例必须也是三类）
    #   位置：轴外、x 轴标签下方、单列竖排。
    #   为什么不能在轴内：BH p = 0.05 竖线落在轴内 61% 处，而三条图例右对齐至少占
    #   面板宽的 49% —— 轴内任何位置都会横跨该虚线（左上/左下又被 2e-4 与
    #   0.027–0.21 的数据占住）。45 mm 宽也放不下单行三条，故竖排。
    #   标记尺寸与轴内一致（空心环 4.6 > 实心点 3.0），避免图例出现第三种尺寸。
    ax.plot([], [], marker='o', ms=4.6, mfc='white', mec=fs.GREY, ls='none',
            label='Before HPA filter')
    ax.plot([], [], marker='o', ms=3.0, color=fs.CATEGORICAL[0], ls='none',
            label='After filter, BH $p$ < 0.05')
    ax.plot([], [], marker='o', ms=3.0, color=fs.GREY, ls='none',
            label='After filter, BH $p$ > 0.05')
    ax.legend(loc='upper left', bbox_to_anchor=(0.0, -0.215), ncol=1, frameon=False,
              fontsize=5.3, handletextpad=0.35, borderpad=0.0, labelspacing=0.28)
    ax.text(0.0, 1.04, 'c', transform=ax.transAxes, ha='left', va='bottom',
            fontsize=8, fontweight='bold', color=fs.BLACK)
    print('  c HPA 前后 BH p：%s'
          % ', '.join('%s %.4f→%.4f' % (k, m[k]['bh'], m[k]['bh_hpa']) for k in order))


def main():
    m, order = load()
    print('  读入 7 类')
    assert len(order) == 7
    # 差值检验 p（H0: r_疲劳 = r_负荷，配对置换，上游 code_celltype_specificity_contrast_20260915.py）
    # Table 2 降为补充表后，这组数出现在补充表与 Fig. 5 图注 → 此处打印以供货件取值核对。
    _pth = REV + '/downstream_newmasks_20260915/celltype_specificity_contrast_FULLPREC_20260919.csv'
    if os.path.exists(_pth):
        _ct = {r['cell_type']: r for r in rd(_pth)}
        print('  差值检验 p（疲劳 − 负荷）：%s'
              % ', '.join('%s %.3f' % (k, float(_ct[k]['p_diff'])) for k in order))
    else:
        print('  ⚠️ 未找到 %s，跳过差值检验打印' % os.path.basename(_pth))
    fig = plt.figure(figsize=(FIG_W * fs.MM, FIG_H * fs.MM))
    panel_a(fig, m, order)
    panel_b(fig, m, order)
    panel_c(fig, m, order)

    fig.canvas.draw()
    fb = fig.get_window_extent()
    rend = fig.canvas.get_renderer()
    over = []
    for ai, ax in enumerate(fig.axes):
        x0, x1 = ax.get_xlim()
        y0, y1 = ax.get_ylim()
        items = list(ax.texts) + [ax.xaxis.label, ax.yaxis.label]
        items += [t for t in ax.get_xticklabels()
                  if t.get_visible() and x0 - 1e-9 <= t.get_position()[0] <= x1 + 1e-9]
        items += [t for t in ax.get_yticklabels()
                  if t.get_visible() and y0 - 1e-9 <= t.get_position()[1] <= y1 + 1e-9]
        leg = ax.get_legend()
        if leg is not None:
            items += list(leg.get_texts())
        for t in items:
            if not t.get_visible() or not t.get_text().strip():
                continue
            bb = t.get_tightbbox(rend)
            if bb is None:
                continue
            out = max(fb.x0 - bb.x0, fb.y0 - bb.y0, bb.x1 - fb.x1, bb.y1 - fb.y1)
            if out > -1.0:
                over.append(('ax%d' % ai, t.get_text()[:32].replace('\n', '|'), round(out, 1)))
    print('  越界文字：%d 个' % len(over))
    for o in over:
        print('     %-6s %-34s 越界 %.1f px' % o)

    # ---- 护栏③（实测）：b 面板的 p 标签不得越出轴框，否则会侵入相邻的 c 面板 ----
    # 用真实渲染 bbox，不靠字符宽度估算（估算会随字号/字体/语言漂移）。
    _axb = fig.axes[1]
    _x0d = _axb.transData.transform((_axb.get_xlim()[0], 0))[0]
    _x1d = _axb.transData.transform((_axb.get_xlim()[1], 0))[0]
    _bad = []
    for _t in getattr(_axb, '_r_labels', []):
        _bb = _t.get_window_extent(rend)
        if _bb.x1 > _x1d + 0.5 or _bb.x0 < _x0d - 0.5:
            _bad.append('%s (%.1f px 越界)' % (_t.get_text(), max(_bb.x1 - _x1d, _x0d - _bb.x0)))
    assert not _bad, 'b 面板 p 标签越出轴框（会侵入相邻面板）：%s' % _bad
    print('     ✓ b 面板 %d 个 p 标签实测全部在轴框内' % len(getattr(_axb, '_r_labels', [])))

    mpl.rcParams['savefig.bbox'] = None
    fs.save_all(fig, 'Fig5_celltype', outdir=HERE,
                formats=('pdf', 'svg', 'eps', 'png'), bbox_inches=None)
    from PIL import Image
    w, h = Image.open(os.path.join(HERE, 'Fig5_celltype.png')).size
    print('  成品 %.1f × %.1f mm（目标 %.1f × %.1f）'
          % (w / 300 * 25.4, h / 300 * 25.4, FIG_W, FIG_H))
    print('  完成')


if __name__ == '__main__':
    main()
