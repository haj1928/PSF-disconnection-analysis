# -*- coding: utf-8 -*-
"""
fig1_lesion_distribution_20260916.py —— Fig. 1 病灶空间分布（三分位叠加概率图）

图型决策（scipilot 选图逻辑）
  * 数据形态：124 例二值掩码 → 任务为「空间分布 + 组间对比」，非数值分布
  * 候选：① 代表例单例图 → 信息量低、易被质疑选择性展示；② 轮廓叠加 → 124 层不可读；
          ③ 组内重叠概率图（overlay probability map）→ 卒中病灶定位研究的标准做法
  * 选定：③。按病灶体积**等例数三分位**（41/42/41）分 3 层，各出一张组内重叠概率图，
          三层**共享色标**（绝对概率，可跨层比较）

版式
  3 行（小/中/大）× 5 轴位层面（z = 0/+10/+20/+30/+40 mm）
  行标签在左、共享垂直色标在右；尺寸 Elsevier 双栏 190 mm
  配色 academic-figure-skill SEQUENTIAL（COPY VERBATIM），透明度随概率递增

坐标约定
  掩码 affine：x = 90 − i，y = −126 + j，z = −72 + k
  显示为**神经学方位**：左半球在左（x 向右增大），前部在上（y 向上增大）
"""
import csv
import io
import os
import sys

import nibabel as nib
import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.offsetbox import AnnotationBbox, HPacker, TextArea, VPacker
from matplotlib.ticker import FuncFormatter

import figstyle_20260916 as fs

HERE = os.path.dirname(os.path.abspath(__file__))
fs.apply_style(base_size=7.0)

# ------------------------------------------------------------------ 参数
SLICE_Z = [0, 10, 20, 30, 40]            # MNI z (mm)
CROP_X = (-74.0, 74.0)                   # MNI x 显示范围
CROP_Y = (-108.0, 68.0)                  # MNI y 显示范围
TEMPLATE = ('C:/Users/Administrator/.workbuddy/binaries/python/envs/default/'
            'Lib/site-packages/nilearn/datasets/data/'
            'mni_icbm152_t1_tal_nlin_sym_09a_converted.nii.gz')

FIG_W_MM = 190.0
LABEL_COL_MM = 29.0
SLICE_W_MM = 30.4
CBAR_W_MM = 70.0
CBAR_H_MM = 3.2
CBAR_X_MM = 29.0
CBAR_Y_MM = 5.0
ROW_PITCH_MM = 37.6
GRID_BOTTOM_MM = 11.0


def mm_to_index(aff, x=0.0, y=0.0, z=0.0):
    return (int(round((x - aff[0, 3]) / aff[0, 0])),
            int(round((y - aff[1, 3]) / aff[1, 1])),
            int(round((z - aff[2, 3]) / aff[2, 2])))


def crop_slice(arr, aff, z):
    """取轴位层面 → 显示数组 D[y, x]（神经学方位：左半球在左、前部在上）"""
    i_a, i_b = mm_to_index(aff, CROP_X[0])[0], mm_to_index(aff, CROP_X[1])[0]
    j_a, j_b = mm_to_index(aff, y=CROP_Y[0])[1], mm_to_index(aff, y=CROP_Y[1])[1]
    k = mm_to_index(aff, z=z)[2]
    lo_i, hi_i = min(i_a, i_b), max(i_a, i_b)
    sl = arr[lo_i:hi_i + 1, j_a:j_b + 1, k]      # (x 递减, y 递增)
    return sl[::-1, :].T                          # → (y 递增, x 递增)


def row_label(name, rng, n):
    """行标签：n 用斜体（统计符号惯例），其余正体"""
    top = TextArea(name + '  (' + rng + ')',
                   textprops=dict(fontsize=7, color=fs.BLACK))
    it = TextArea('n', textprops=dict(fontsize=7, fontstyle='italic', color=fs.BLACK))
    eq = TextArea(' = %d' % n, textprops=dict(fontsize=7, color=fs.BLACK))
    return VPacker(children=[top, HPacker(children=[it, eq], align='baseline', pad=0, sep=0)],
                   align='right', pad=0, sep=1.0)


def main():
    summ = list(csv.DictReader(io.open(os.path.join(HERE, 'fig1_tertile_summary.csv'),
                                       encoding='utf-8-sig')))
    probs, aff = [], None
    for s in summ:
        img = nib.load(os.path.join(HERE, s['out_nii']))
        probs.append(np.asarray(img.dataobj, dtype=np.float32))
        aff = np.asarray(img.affine)
        print('  读入 %s  n=%s  峰值 %.4f' % (s['out_nii'], s['n'], float(s['max_prob'])))

    from nibabel.processing import resample_from_to
    tpl = np.asarray(resample_from_to(nib.load(TEMPLATE), (probs[0].shape, aff),
                                      order=1).dataobj, dtype=np.float32)
    tpl = np.clip(tpl, 0, 255) / 255.0
    print('  模板重采样 shape=%s' % (tpl.shape,))

    vmax = max(float(np.max(p)) for p in probs)
    print('  共享色标上限 vmax = %.4f' % vmax)

    d0 = crop_slice(probs[0], aff, 0)
    print('  显示数组 shape = %s（应为 (y, x) = (177, 149)）' % (d0.shape,))
    asp = d0.shape[0] / float(d0.shape[1])           # h / w
    print('  切片高宽比 = %.4f' % asp)

    # -------------------------------------------------------------- 版面（全部按 mm 精确定位）
    slice_h = SLICE_W_MM * asp
    fig_h = GRID_BOTTOM_MM + 2 * ROW_PITCH_MM + slice_h + 6.0
    fig = plt.figure(figsize=(FIG_W_MM * fs.MM, fig_h * fs.MM))
    print('  版面 %.1f × %.1f mm；单切片 %.1f × %.1f mm' % (FIG_W_MM, fig_h, SLICE_W_MM, slice_h))

    def rect(x_mm, y_mm, w_mm, h_mm):
        return [x_mm / FIG_W_MM, y_mm / fig_h, w_mm / FIG_W_MM, h_mm / fig_h]

    cmap = fs.CMAP_SEQ.copy()
    cmap.set_under(alpha=0.0)
    # 覆盖用色标：`fs.CMAP_OVERLAY` = SEQUENTIAL 抬高冷端的版本（#F7FBFF 近白，
    # 叠在浅灰模板上不可辨）→ 见 QA 报告「有意偏离」。
    # 2026-09-17：原在此处写死色值，现收进统一色板（颜色一字未改，仅集中管理）。
    cmap_ov = fs.CMAP_OVERLAY
    norm = mpl.colors.Normalize(vmin=0.0, vmax=vmax)

    labels = [('Small', '< 1.2 mL'), ('Medium', '1.2–6.5 mL'), ('Large', '> 6.5 mL')]

    row_bboxes = []
    for r in range(3):
        y0 = GRID_BOTTOM_MM + (2 - r) * ROW_PITCH_MM
        ax_first = None
        for c, z in enumerate(SLICE_Z):
            x0 = LABEL_COL_MM + c * SLICE_W_MM
            ax = fig.add_axes(rect(x0, y0, SLICE_W_MM, slice_h))
            if c == 0:
                ax_first = ax
            t = crop_slice(tpl, aff, z)
            p = crop_slice(probs[r], aff, z)
            ax.imshow(np.ma.masked_where(t <= 0.11, t), cmap='gray',
                      vmin=0.06, vmax=0.80, alpha=0.60,
                      origin='lower', interpolation='nearest', aspect='auto',
                      rasterized=True)
            frac = np.clip(p / vmax, 0, 1)
            rgba = cmap_ov(frac)
            # 固定不透明度：概率**只由颜色**编码，与色标一一对应（避免 alpha+颜色双重编码）
            rgba[..., 3] = np.where(p > 0, 0.75, 0.0)
            ax.imshow(rgba, origin='lower', interpolation='nearest', aspect='auto',
                      rasterized=True)
            ax.set_axis_off()
            if r == 0:
                ax.text(0.5, 1.01, 'z = %+d' % z if z else 'z = 0',
                        transform=ax.transAxes, ha='center', va='bottom',
                        fontsize=6.8, color=fs.BLACK)
        # 行标签 + 面板字母（锚定该行**第一**个子图）
        ab = AnnotationBbox(row_label(labels[r][0], labels[r][1], int(summ[r]['n'])),
                            (-0.05, 0.5), xycoords='axes fraction',
                            box_alignment=(1.0, 0.5), frameon=False,
                            annotation_clip=False)
        ax_first.add_artist(ab)
        row_bboxes.append(ab)
        ax_first.text(-0.05, 1.045, 'abc'[r], transform=ax_first.transAxes,
                      ha='right', va='bottom', fontsize=8, fontweight='bold',
                      color=fs.BLACK)

    # -------------------------------------------------------------- 色标（横向，置于网格下方）
    cax = fig.add_axes(rect(CBAR_X_MM, CBAR_Y_MM, CBAR_W_MM, CBAR_H_MM))
    sm = mpl.cm.ScalarMappable(norm=norm, cmap=cmap_ov)
    cb = fig.colorbar(sm, cax=cax, orientation='horizontal')
    cb.ax.tick_params(labelsize=6.3, width=0.5, length=1.8, pad=1.2)
    cb.outline.set_linewidth(0.5)
    cb.set_ticks([0.0, 0.05, 0.10, 0.15, 0.20])
    cb.ax.xaxis.set_major_formatter(
        FuncFormatter(lambda v, _: ('%.2f' % v).rstrip('0').rstrip('.')))
    fig.text(CBAR_X_MM / FIG_W_MM, (CBAR_Y_MM + CBAR_H_MM + 1.0) / fig_h,
             'Lesion overlap probability', ha='left', va='bottom',
             fontsize=6.8, color=fs.BLACK)

    # -------------------------------------------------------------- 数据自检
    bad = 0
    for r in range(3):
        for c, z in enumerate(SLICE_Z):
            p = crop_slice(probs[r], aff, z)
            t = crop_slice(tpl, aff, z)
            if p.shape != t.shape:
                print('  [!] 行%d 列%d 形状不一致' % (r, c)); bad += 1
            if not np.isfinite(p).all():
                print('  [!] 行%d 列%d 含非有限值' % (r, c)); bad += 1
    print('  数值/形状自检：%s' % ('PASS' if bad == 0 else 'FAIL(%d)' % bad))
    print('  三分位峰值：%s（共享 vmax=%.4f）' % ([round(float(p.max()), 4) for p in probs], vmax))

    # -------------------------------------------------------------- 版式自检（文字不出界）
    fig.canvas.draw()
    fb = fig.get_window_extent()
    qa_artists = []
    for r in range(3):                       # 行标签盒（用 AnnotationBbox 自身 extent）
        qa_artists.append(('row-label-%s' % 'abc'[r], row_bboxes[r]))
    for t in list(fig.texts) + [t for ax in fig.axes for t in ax.texts] \
            + list(cb.ax.get_xticklabels()):
        if t.get_visible() and t.get_text().strip():
            qa_artists.append((t.get_text()[:34], t))
    over = []
    rend = fig.canvas.get_renderer()
    for tag, art in qa_artists:
        try:
            bb = (art.get_tightbbox(rend) if hasattr(art, 'get_tightbbox')
                  else art.get_window_extent(rend))
        except Exception:
            try:
                bb = art.get_window_extent()
            except Exception:
                continue
        if bb is None:
            continue
        out = max(fb.x0 - bb.x0, fb.y0 - bb.y0, bb.x1 - fb.x1, bb.y1 - fb.y1)
        if out > -1.0:
            over.append((tag, round(bb.x0 - fb.x0, 1), round(bb.y0 - fb.y0, 1),
                         round(bb.x1 - fb.x0, 1), round(bb.y1 - fb.y0, 1), round(out, 1)))
    print('  越界文字：%d 个（阈值 −1.0 px @ 画布 dpi）' % len(over))
    for o in over:
        print('     %-36s x %.1f–%.1f  y %.1f–%.1f  越界 %.1f px' % o)
    print('  画布 %.0f × %.0f px' % (fb.x1, fb.y1))
    # 行标签实测宽度（mm）
    for r in range(3):
        bb = row_bboxes[r].get_tightbbox(rend)
        print('  行标签 %s：x %.1f–%.1f mm（左边界应 ≥ 0）'
              % ('abc'[r], bb.x0 / fb.x1 * FIG_W_MM, bb.x1 / fb.x1 * FIG_W_MM))

    # 投稿用精确版面尺寸：关闭 tight 裁剪（版式自检已确认 0 越界，无需裁剪）
    mpl.rcParams['savefig.bbox'] = None
    fs.save_all(fig, 'Fig1_lesion_distribution', outdir=HERE,
                formats=('pdf', 'svg', 'eps', 'png'), bbox_inches=None)
    png = os.path.join(HERE, 'Fig1_lesion_distribution.png')
    from PIL import Image
    w_px, h_px = Image.open(png).size
    print('  成品像素 %d × %d → 折算 %.1f × %.1f mm @300dpi（目标 %.1f × %.1f）'
          % (w_px, h_px, w_px / 300 * 25.4, h_px / 300 * 25.4, FIG_W_MM, fig_h))
    print('  完成')


if __name__ == '__main__':
    main()
