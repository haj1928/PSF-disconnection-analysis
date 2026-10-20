# -*- coding: utf-8 -*-
"""
Fig. 2  行为 PLS 结果（皮层 100 区 · 表面渲染）
2026-09-16 | 数据：PSF_pls_discon_salience_residY_age_sex_cortex100_newmasks_20260915.csv

面板
  a  LV1 载荷在 100 个皮层区的分布（Schaefer100 fsa5 表面渲染）
  b  |BSR| > 2 的稳健区域分布（BSR 由 5,000 次 bootstrap, seed = 0 估计）

图型决策（技能 Step 1）
  * 数据形态：1 个连续变量定义在 100 个皮层区上 → 空间分布图（脑表面）；**不用**折线/柱状，
    因为 x 轴是空间而非有序量（scipilot《chart_selection》：分类不同量纲不得连线）
  * **色标必须与数据类型匹配**，两个面板不同：
      - a 面板载荷有正有负（−0.062 ~ +0.201，20/100 区为负）→ **双向色标**（蓝−白−红，以 0 为中心）
      - b 面板只显示 24 个稳健区，其 BSR **全为正**且都 > 2 → **单极强度色标**（浅红→深红，
        范围 [2.0, 3.2]），只表达"超过阈值的强度"，色标上不出现 |BSR| ≤ 2 的任何区间
  * ⚠️ 2026-09-17 修正：原先两个面板的**脑图**都用双向色标，但**色条**却画成了顺序色标
    （`cmap=fs.CMAP_SEQ`），导致"色条只有白→蓝、脑图却大片红"的自相矛盾；且 b 面板色条
    跨 −3.12~+3.12，中段白色正对应被排除的 |BSR| < 2 区间。现两者已统一。
  * 稳健区域以外的区（|BSR| ≤ 2）在 b 面板中置 NaN、**不渲染**（露出灰底），
    图注亦注明"仅显示 |BSR| > 2 的 24 个稳健区"

输出：PDF / SVG / EPS / 可编辑 PPTX / 300 dpi PNG（见 figstyle_20260916.py）
"""
import csv
import io
import os
import sys

import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
import nibabel as nb
from nilearn import datasets, plotting

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import figstyle_20260916 as fs

HERE = os.path.dirname(os.path.abspath(__file__))
BASE = 'F:/dcmtonii/run62_crosssectional'
SAL = BASE + '/PSF_pls_behavior/PSF_pls_discon_salience_residY_age_sex_cortex100_newmasks_20260915.csv'
ANNOT = BASE + '/PSF_repro_PET/assets/schaefer100_fsa5'
OUT = HERE

TITLE = ('LV1 loading' if False else None)


def load_parcels():
    rows = list(csv.DictReader(io.open(SAL, encoding='utf-8-sig')))
    assert len(rows) == 100, '载荷文件应为 100 行，实为 %d' % len(rows)
    names = [r['region'] for r in rows]
    sal = np.array([float(r['salience']) for r in rows])
    bsr = np.array([float(r['BSR']) for r in rows])
    sig = np.array([r['sig'].strip().lower() == 'true' for r in rows])
    return names, sal, bsr, sig


def annot_vertex_map():
    """返回 {parcel_name: 顶点索引数组}（左右半球分开）"""
    out = {}
    for hemi, fn in (('LH', 'lh.Schaefer2018_100Parcels_7Networks_order.annot'),
                     ('RH', 'rh.Schaefer2018_100Parcels_7Networks_order.annot')):
        labels, ctab, names = nb.freesurfer.io.read_annot(os.path.join(ANNOT, fn))
        for lab, nm in zip(np.unique(labels), [None] * len(np.unique(labels))):
            pass
        for i, nm in enumerate(names):
            key = nm.decode() if isinstance(nm, bytes) else nm
            key = key.replace('7Networks_', '')
            idx = np.where(labels == i)[0]
            if idx.size:
                out[key] = idx
    return out


def vertex_map(names, values, vmap, n_vert):
    out = np.full(n_vert, np.nan)
    for nm, v in zip(names, values):
        idx = vmap.get(nm)
        assert idx is not None, '分区名未匹配：%s' % nm
        out[idx] = v
    return out



def clean_bgmap(sulc_path, aparc_path):
    """把内侧壁（aparc 标签 0 = Unknown/medial wall）的 sulc 置为中性值，
    再做「背景板 L2」亮度映射：sulc（2–98 百分位）线性映射进 [0.10, 0.38]。
    背景：nilearn `Greys` 为 0=白 / 1=黑，映射后背景最深只到浅灰，避免
    inflated 表面内侧视图出现黑色团块；nilearn 默认不做内侧壁处理。"""
    import nibabel as _nb
    bg = np.asarray(_nb.load(sulc_path).darrays[0].data, dtype=float)
    labels = _nb.freesurfer.io.read_annot(aparc_path)[0]
    neutral = float(np.median(bg[labels != 0]))
    bg2 = bg.copy()
    bg2[labels == 0] = neutral
    v = bg2[np.isfinite(bg2)]
    p2, p98 = np.percentile(v, 2.0), np.percentile(v, 98.0)
    t = np.clip((bg2 - p2) / (p98 - p2), 0.0, 1.0)
    return 0.10 + t * (0.38 - 0.10)

# ---------------------------------------------------------------- 版面常量
# 供 PPTX 版（fig2_make_editable_pptx_20260917.py）复用，避免两处几何各自漂移。
FIG_H_MM = 136.0
GRID = dict(height_ratios=[1.0, 1.0, 0.62], left=0.035, right=0.995,
            top=0.975, bottom=0.068, wspace=0.0, hspace=0.02)
VIEWS = [('left', 'lateral'), ('right', 'lateral'), ('left', 'medial'), ('right', 'medial')]
# 2026-09-17 修正：原 a/b 用 text2D(-0.12) → 字形落在距左边缘 0.5–1.2 mm（几乎压边），
# 而 c 用 text(-0.32) → 23.1 mm，两者相差 21.9 mm，三个字母既不齐又贴边。
# 现统一为**同一图坐标 x + 各自对齐本行内容顶边**。
LETTER_X, LETTER_PAD = 0.016, 0.007     # 距左边缘约 3.0 mm；字母顶边高出内容顶边约 0.9 mm
CONTENT_DY = {0: 0.040, 1: 0.0385, 2: 0.0043}   # 轴顶→内容顶的实测偏移（3D 轴留白大、2D 轴几乎无）
CB_X, CB_W, CB_H = 0.42, 0.16, 0.014
BSR_LO, BSR_HI = 2.0, 3.2               # b 面板色标只覆盖稳健区实际取值区间


def prepare():
    """数据 → 顶点图 + 表面几何。图版与 PPTX 版共用，避免各写一份。"""
    names, sal, bsr, sig = load_parcels()
    print('=== Fig. 3 自检 ===')
    print('  载荷：n=%d，范围 %.3f–%.3f（负值 %d 个）'
          % (len(sal), sal.min(), sal.max(), int((sal < 0).sum())))
    print('  稳健区（|BSR|>2）：%d 个（正文口径 24）' % int(sig.sum()))
    bsr_sig = np.array([b for b, s in zip(bsr, sig) if s])
    print('  a 面板：载荷含负值 %d 个 → **双向色标**（蓝−白−红，以 0 为中心，'
          '显示范围 ±%.3f）' % (int((sal < 0).sum()), np.abs(sal).max()))
    print('  b 面板：稳健区 BSR 全为正 %s（%s）→ **单极强度色标**（浅红→深红，范围 2.0–3.2）'
          % (bool((bsr_sig > 0).all()), '%.3f–%.3f' % (bsr_sig.min(), bsr_sig.max())))
    assert (bsr_sig > 0).all(), 'b 面板采用单极红色标，要求稳健区 BSR 全为正；若出现负值需改回双向'
    assert int(sig.sum()) == 24, '稳健区数 %d ≠ 声明值 24' % int(sig.sum())

    vmap = annot_vertex_map()
    fsa = datasets.fetch_surf_fsaverage('fsaverage5')
    BG = {}
    for hemi in ('left', 'right'):
        sulc = fsa['sulc_%s' % hemi]
        aparc = BASE + '/PSF_repro_PET/assets/fsa5_%s_aparc.annot' % ('lh' if hemi == 'left' else 'rh')
        BG[hemi] = clean_bgmap(sulc, aparc)
    out = {}
    for hemi, hk in (('left', 'LH'), ('right', 'RH')):
        n_v = nb.freesurfer.io.read_annot(
            os.path.join(ANNOT, ('lh' if hemi == 'left' else 'rh') +
                         '.Schaefer2018_100Parcels_7Networks_order.annot'))[0].shape[0]
        sub_names = [n for n in names if n.startswith(hk)]
        idx = [i for i, n in enumerate(names) if n.startswith(hk)]
        out[hemi] = dict(
            n_v=n_v,
            sal=vertex_map(sub_names, sal[idx], vmap, n_v),
            bsr_sig=vertex_map(sub_names, np.where(sig[idx], bsr[idx], np.nan), vmap, n_v),
        )
    return dict(names=names, sal=sal, bsr=bsr, sig=sig, bsr_sig=bsr_sig,
                fsa=fsa, BG=BG, out=out,
                vmax_sal=float(np.nanmax(np.abs(sal))),
                vmax_bsr=float(np.nanmax(np.abs(bsr_sig))))


def panel_specs(vmax_sal):
    """两个面板的「色标 + 显示范围」——色条与脑图共用同一份，结构上不会再分叉。"""
    return (
        dict(key='sal', tag='a', cmap=fs.CMAP_DIV, vmin=-vmax_sal, vmax=vmax_sal,
             thr=0.0, sym=True, cbtitle='LV1 loading', cb_dy=0.0241),
        dict(key='bsr_sig', tag='b', cmap=fs.CMAP_SEQ_RED, vmin=BSR_LO, vmax=BSR_HI,
             thr=BSR_LO, sym=False, cbtitle='BSR (|BSR| > 2)', cb_dy=0.0273),
    )


def main():
    fs.apply_style()
    D = prepare()
    names, sal, bsr, sig = D['names'], D['sal'], D['bsr'], D['sig']
    fsa, BG, out = D['fsa'], D['BG'], D['out']

    # 高度 130 → 136 mm、bottom 0.035 → 0.068（2026-09-17）：原 bottom 只给轴下方留 4.6 mm，
    # 而「刻度数值(6.5pt) + x 轴标签(7pt) + 各处 padding」实测需约 6.4 mm → x 轴标签最下缘
    # 距图底仅 0.25 mm（几乎压边）。加高 6 mm 并抬高轴，刚好腾出约 2.9 mm 下边距，
    # 且**三行面板的绝对高度不变**（脑图不缩水）。
    fig = plt.figure(figsize=(fs.W_DOUBLE * fs.MM, FIG_H_MM * fs.MM))
    gs = fig.add_gridspec(3, 4, **GRID)
    views = VIEWS
    # b 面板只画稳健强度 → 色标范围就是稳健区实际取值区间 [2, 3.2]，不出现 |BSR| ≤ 2 的任何区间

    # 色条与脑图**共用同一套 cmap 与同一组 vmin/vmax**（见文件头 2026-09-17 修正说明）
    PANELS = panel_specs(D['vmax_sal'])

    row_pos = {}
    for row, P in enumerate(PANELS):
        for col, (hemi, view) in enumerate(views):
            ax = fig.add_subplot(gs[row, col], projection='3d')
            plotting.plot_surf_stat_map(
                surf_mesh=fsa['infl_%s' % hemi], stat_map=out[hemi][P['key']],
                bg_map=BG[hemi], hemi=hemi, view=view,
                cmap=P['cmap'], vmin=P['vmin'], vmax=P['vmax'], threshold=P['thr'],
                symmetric_cbar=P['sym'], colorbar=False, axes=ax, figure=fig,
            )
            ax.set_title('')  # 标题留白；视野由 3D 自动
            if col == 0:
                row_pos[row] = ax.get_position()
    # 面板字母 a / b：同一 x，各自贴本行内容顶边
    for row, P in enumerate(PANELS):
        pos = row_pos[row]
        fig.text(LETTER_X, pos.y1 - CONTENT_DY[row] + LETTER_PAD, P['tag'],
                 ha='left', va='top', fontsize=9, fontweight='bold', color=fs.BLACK)
    # 每行一个色条（位置跟随本行轴区下缘）
    for row, P in enumerate(PANELS):
        cby = row_pos[row].y0 + P['cb_dy']
        cax = fig.add_axes([CB_X, cby, CB_W, CB_H])
        norm = mpl.colors.Normalize(vmin=P['vmin'], vmax=P['vmax'])
        cb = fig.colorbar(mpl.cm.ScalarMappable(norm=norm, cmap=P['cmap']), cax=cax,
                          orientation='horizontal')
        cb.ax.tick_params(labelsize=6, width=0.6, length=2)
        cb.outline.set_linewidth(0.5)
        cax.set_title(P['cbtitle'], fontsize=6.5, pad=2)

    # ---------------- c 面板：24 个稳健区的网络归属（2026-09-16 新增）----------------
    cnt = {}
    for n_, s_ in zip(names, sig):
        if not s_:
            continue
        cnt[n_.split('_')[1]] = cnt.get(n_.split('_')[1], 0) + 1
    # 排序：先按区数降序，同数时按正文列举顺序（默认模式 → 背侧注意 → 额顶控制 → 躯体运动 → 腹侧注意/显著）
    ORD = {'Default': 0, 'DorsAttn': 1, 'Cont': 2, 'SomMot': 3,
           'SalVentAttn': 4, 'Limbic': 5, 'Vis': 6}
    items = sorted(cnt.items(), key=lambda kv: (-kv[1], ORD.get(kv[0], 99)))
    print('  c 网络归属：' + '、'.join('%s %d' % (k_, v_) for k_, v_ in items))
    assert sum(v for _, v in items) == 24, 'c 面板计数 %d ≠ 24' % sum(v for _, v in items)
    NLAB = {'Default': 'Default mode', 'DorsAttn': 'Dorsal attention',
            'Cont': 'Frontoparietal control', 'SomMot': 'Somatomotor',
            'SalVentAttn': 'Ventral attention / salience',
            'Limbic': 'Limbic', 'Vis': 'Visual'}
    axc = fig.add_subplot(gs[2, 1:3])
    y = np.arange(len(items))[::-1]
    axc.barh(y, [v for _, v in items], height=0.6,
             color=fs.CATEGORICAL[0], edgecolor='none')
    axc.set_yticks(y)
    axc.set_yticklabels([NLAB.get(k_, k_) for k_, _ in items], fontsize=6.5)
    for yy_, (k_, v_) in zip(y, items):
        axc.text(v_ + 0.16, yy_, str(v_), va='center', ha='left',
                 fontsize=6.2, color=fs.BLACK)
    # labelpad 1.5（原为默认 4.0）：底部留白已由 gridspec bottom 补足，
    # 标签同时上提一点，避免整块"刻度数值 + 标签"压在图底。
    axc.set_xlabel('Parcels with |BSR| > 2', fontsize=7, labelpad=1.5)
    axc.set_xlim(0, max(v for _, v in items) + 1.3)
    axc.tick_params(axis='x', labelsize=6.5, width=0.6, length=2)
    axc.tick_params(axis='y', length=0)
    for sp in ('top', 'right'):
        axc.spines[sp].set_visible(False)
    # 面板字母 c：与 a/b 同一图坐标 x，对齐 c 面板内容顶边（原为 -0.32/1.04 → 23.1 mm 处，与 a/b 差 21.9 mm）
    posc = axc.get_position()
    fig.text(LETTER_X, posc.y1 - CONTENT_DY[2] + LETTER_PAD, 'c',
             ha='left', va='top', fontsize=9, fontweight='bold', color=fs.BLACK)

    # bbox_inches=None（2026-09-17）：原用默认的 'tight' 会**按内容裁剪**，于是最靠边的元素
    # （面板字母）永远只剩约 0.35 mm 边距——把字母往右移 3 mm 也没用，因为裁掉的正是那段余白。
    # ⚠️ 光传 bbox_inches=None 不够：figstyle 里 `savefig.bbox='tight'` 会在参数为 None 时
    #    接管（matplotlib 的回落规则）→ 必须同时关掉 rcParam。其余 8 张图早就这么做了，
    #    只有本脚本漏掉，这正是 Fig. 3 成品被裁边的原因。
    mpl.rcParams['savefig.bbox'] = None
    fs.save_all(fig, 'Fig2_behavioral_PLS', outdir=OUT, bbox_inches=None)
    print('  布局：a/b 两行各 4 视野（左外、右外、左内、右内）+ c 网络归属；'
          '与图注一致（Fig. 3. a 载荷 · b 稳健区 · c 网络归属）')
    return 0


if __name__ == '__main__':
    sys.exit(main())
