# -*- coding: utf-8 -*-
"""
fig3_robustness_20260916.py —— Fig. 3 断连-疲劳共变模式的稳健性

面板（2026-09-24 郝医生裁定：**删 b 面板，全图仅保留单面板散点**）
  a 脑-行为散点（Lx vs Ly，n = 124）——双层带（均数 CI + 个体 PI），右下留白放统计量
  （原 b 面板「多口径置换 p 与 corr 双栏点图」：其 4 个口径的 p 与 corr 已在正文 §3.2
    逐条写出，图面冗余，2026-09-24 删除。口径明细仍保留在提交包
    derived_data/fig3_sensitivity_specs.csv，供复现核对脚本使用。
    更早删除：原 a「碎石图」→「置换零分布」、原 d「逐例剔除曲线」。
    单面板后**不再有面板字母**。）

图型依据（scipilot chart_selection / academic-figure-skill）
  * 2 连续 → 散点 + 回归线（n = 124 ≥ 30）
  * ⚠️ 已实测否掉「分箱均值 ± CI」聚合方案：4 个四分位 CI 全重叠（Q1 vs Q4 Welch
    p = 0.176、Cohen d = 0.35），n 124 → 4 会摧毁这个中小效应

数据（全部由脚本读取，无硬编码）
  fig3_subject_scores.csv      —— 个体 Lx / Ly / 体积
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
fs.apply_style(base_size=7.0)

# ---------------------------------------------------------------- 版面（mm，左下原点）
FIG_W, FIG_H = 190.0, 97.0
# 单面板（2026-09-24 删 b 后）：下留 12 mm 放 x 刻度与轴标签，上留 5 mm。
# 散点高取 80 mm 使纵横比 ≈ 1.8（70 mm 时过扁，两个带会被拉成大楔形、视觉上压过数据）。
A_RECT = (34.0, 12.0, 144.0, 80.0)      # a：散点


def rect(x, y, w, h):
    return [x / FIG_W, y / FIG_H, w / FIG_W, h / FIG_H]


def read_csv(p):
    return list(csv.DictReader(io.open(p, encoding='utf-8-sig')))


# ================================================================ a：脑-行为散点
def panel_a(fig, subj):
    ax = fig.add_axes(rect(*A_RECT))
    lx = np.array([float(r['Lx']) for r in subj])
    ly = np.array([float(r['Ly']) for r in subj])
    r = float(np.corrcoef(lx, ly)[0, 1])
    n = len(lx)
    ax.scatter(lx, ly, s=7.0, facecolor=fs.CATEGORICAL[0], edgecolor='none', alpha=0.55,
               zorder=3)
    # 回归线 + 两类「带」（2026-09-24 郝医生裁定增补）
    #   内层 = 均数 95% 置信带：只描述【回归线本身】的不确定度
    #   外层 = 个体 95% 预测带：描述【单个被试观测值】的散布范围（半宽约为内层的 8.7 倍）
    #   ⚠️ 原图只画内层 ⇒ 实测 96.8% 的散点落在带外，读图者会误判为"数据与趋势不一致"。
    #      根因：se 公式缺了个体残差方差项（+1）。两层都画，读者才能同时读到
    #      "线有多准"与"点该散多开"，散点几乎全落在预测带内。
    b1, b0 = np.polyfit(lx, ly, 1)
    xs = np.linspace(lx.min(), lx.max(), 120)
    yhat = b1 * xs + b0
    dof = n - 2
    resid = ly - (b1 * lx + b0)
    s2 = float(resid @ resid) / dof
    xbar = lx.mean()
    sxx = float(((lx - xbar) ** 2).sum())
    tcrit = 1.9799                        # t(122, 0.975)，n = 124
    se_ci = np.sqrt(s2 * (1.0 / n + (xs - xbar) ** 2 / sxx))          # 均数 CI
    se_pi = np.sqrt(s2 * (1.0 + 1.0 / n + (xs - xbar) ** 2 / sxx))    # 个体 PI（多 +1）
    ax.fill_between(xs, yhat - tcrit * se_pi, yhat + tcrit * se_pi, color=fs.ACCENT_RED,
                    alpha=0.08, lw=0, zorder=1)                        # 外层：个体预测带
    ax.fill_between(xs, yhat - tcrit * se_ci, yhat + tcrit * se_ci, color=fs.ACCENT_RED,
                    alpha=0.20, lw=0, zorder=2)                        # 内层：均数置信带
    ax.plot(xs, yhat, color=fs.ACCENT_RED, lw=0.9, zorder=4)
    # 顶部预留标注空间：高杠杆个体位于图幅右上角，其上方原本无数据
    ylo, yhi = ly.min() - 0.45, ly.max() + 0.60   # 2026-09-22：不再为最大病灶标注留白
    ax.set_ylim(ylo, yhi)
    ax.set_xlim(lx.min() - 0.06, lx.max() + 0.06)
    # 2026-09-22：最大病灶（高杠杆）个体的空心圈与标注已随该检验一并移除
    ax.set_xlabel('Brain score, $L_x$ (LV1)', labelpad=1.5)
    ax.set_ylabel('Fatigue score, $L_y$ (LV1)', labelpad=1.5)
    ax.tick_params(length=1.8, pad=1.2)
    # 统计量置于右下留白区（该区无散点、无置信带）
    ax.text(0.975, 0.045, 'corr($L_x$, $L_y$) = %+.3f\n$n$ = %d' % (r, n),
            transform=ax.transAxes, ha='right', va='bottom', fontsize=6.8,
            color=fs.BLACK, linespacing=1.45, zorder=6)
    print('  a 散点：corr(Lx,Ly) = %+.4f（应 ≈ +0.3174）' % r)


def main():
    subj = read_csv(os.path.join(HERE, 'fig3_subject_scores.csv'))

    fig = plt.figure(figsize=(FIG_W * fs.MM, FIG_H * fs.MM))
    panel_a(fig, subj)

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
                over.append(('ax%d' % ai, t.get_text()[:34].replace('\n', '|'), round(out, 1)))
    print('  越界文字：%d 个' % len(over))
    for o in over:
        print('     %-6s %-36s 越界 %.1f px' % o)

    mpl.rcParams['savefig.bbox'] = None
    fs.save_all(fig, 'Fig3_robustness', outdir=HERE,
                formats=('pdf', 'svg', 'eps', 'png'), bbox_inches=None)
    from PIL import Image
    w, h = Image.open(os.path.join(HERE, 'Fig3_robustness.png')).size
    print('  成品 %.1f × %.1f mm（目标 %.1f × %.1f）' % (w / 300 * 25.4, h / 300 * 25.4,
                                                         FIG_W, FIG_H))
    print('  完成')


if __name__ == '__main__':
    main()
