# ====================================================================
# v5 色彩契约（2026-10-08，全稿统一到 RdBu_r 取样）
#   本文件由 _palette_v5_build_20261008.py 生成，勿手改。
#   DIVERGING=['#053061', '#3480B9', '#A5CEE3', '#F7F6F6', '#F7B596', '#C6413E', '#67001F']
#   ACCENT_RED=#BA2832  CATEGORICAL[0]=#2A71B2
#   V5_KEEP_SINGLE=False
# ====================================================================
# -*- coding: utf-8 -*-
"""
figstyle_20260916.py —— PSF 断连稿统一出图样式模块
按三个绘图技能的规范合并：

* 期刊尺寸/字体/线宽：Elsevier（NeuroImage）单栏 90 mm / 1.5 栏 140 mm / 双栏 190 mm；
  Arial；正文标签 7–9 pt（下限 6 pt）；线宽 0.6 pt（academic-figure-skill `journal-specs.md`
  只列 Nature 89/183 mm，本稿按 Elsevier 取值——已在出图报告中注明该偏离）
* 配色：`academic-figure-skill/references/color-palettes.md` 的 **COPY VERBATIM** 色板
  （Nature/Cell/Science 语义色），另按 scipilot 的色盲安全要求补充 Okabe-Ito 备选
* 排版：无网格、只留左/下轴脊、图例不压数据（academic-figure-skill 反模式清单 AP-*）
* 导出：矢量优先（PDF/EPS/SVG），位图 300 dpi；`pdf.fonttype=42` + `svg.fonttype='none'`
  保证文字可编辑（郝医生要求「图都要能编辑」）
"""
import os

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap

# ------------------------------------------------------------------ 期刊尺寸（mm）
MM = 1.0 / 25.4
W_SINGLE = 90.0
W_ONEHALF = 140.0
W_DOUBLE = 190.0

# ------------------------------------- 配色（academic-figure-skill COPY VERBATIM）
# 2026-09-21 换色板（郝医生批准 6 色方案）
#   旧：["#2166AC", "#B2182B", "#1B7837", "#F1A340", "#762A83", "#666666"]
#       四种色觉模拟最小 ΔE：正常 43.6 ／ 红 15.5 ／ **绿 12.4** ／ **蓝 14.2**
#       （绿盲与蓝盲低于本稿 ΔE ≥ 15 的判定标准；此前只校验过 CATEGORICAL5，6 色这块漏了）
#   新：H4+I2（左海方案 8 与 9 的组合）→ 正常 32.2 ／ 红 28.3 ／ **绿 23.4** ／ **蓝 23.0**
#   ⚠️ **索引语义必须保住**（本稿按索引取单色，不是 zip 分配）：
#       [0] 蓝 = 主强调色（10+ 处：Fig.1/3/4/5/6/7/8/S1/S2 的显著标记与主数据色）
#       [2] 绿 = 次强调色（Fig. 7c 的 μ-阿片折线、Fig. S1 的区间带与参考线）
#       故 6 色按色相重排为 蓝 → 浅蓝 → 橄榄绿 → 浅黄 → 橙红 → 深青，使 [0] 仍是蓝、[2] 仍是绿。
#   2026-09-21 第二轮（郝医生裁定「全部统一」）：CATEGORICAL5 / DIVERGING / SEQUENTIAL 已一并统一。
#   模拟方法：Machado 2009 线性矩阵 + CIELab ΔE76。
CATEGORICAL = ["#2A71B2", "#72BCD5", "#C4CC99", "#FFD06F", "#E76254", "#01303B"]
# 5 槽类别色板：用于"同一张图里有 5 个分类"的面板（目前只有 Fig. 5a 的 5 条代表通路）。
# 2026-09-21 第二轮：换成 **I3+H2**（左海方案 9 与 8 的组合）。
#   旧 #2166AC#F1A340#009E73#762A83#B2182B → 红 15.5（贴边）／绿 16.3／蓝 21.3
#   新                            → 正常 34.1 ／ 红 **28.3** ／ 绿 30.8 ／ 蓝 23.0，四种色觉全 ≥15
CATEGORICAL5 = ["#E76254", "#C4CC99", "#01303B", "#72BCD5", "#FFD06F"]
#   2026-09-21 第三轮（郝医生裁定「包 B」）：第 0 槽由 #A2221F 换成 **#E76254**，
#   把深红 #A2221F 让给 ACCENT_RED —— 原因：原来 [0]=#A2221F 与强调红 #B2182B 在
#   Fig.6 同图内 ΔE 仅 3.7（近乎同色，C2 判据 FAIL）。换后 C2 3.7 → **22.5**。
CATEGORICAL_EXTENDED = [
    "#2166AC", "#B2182B", "#1B7837", "#F1A340", "#762A83", "#666666",
    "#4393C3", "#D6604D", "#5AAE61", "#B35806", "#9970AB", "#999999",
]
# 2026-09-21 第二轮：左端 #2166AC → 新 CATEGORICAL[0] **#376795**，消除全稿"两套蓝"；
# 红端 = ACCENT_RED（第三轮改为 #A2221F），故发散板红端与强调红**同色**。
#   全稿一种蓝（#376795）、一种红（#A2221F），两者都在师兄 R 文件里。
# 校验：正常 62.6 ／ 红 59.8 ／ 绿 65.0 ／ 蓝 60.3；L* 42→97→38 = 倒 V ✓
DIVERGING = ["#053061", "#3480B9", "#A5CEE3", "#F7F6F6", "#F7B596", "#C6413E", "#67001F"]
# 2026-09-21 第二轮：由 Blues 改为**暖臂**（浅米黄 → 深红）。
# 校验：正常 19.9 ／ 红 18.8 ／ 绿 20.1 ／ 蓝 22.5；L* 86→44→27 单调递减 ✓
SEQUENTIAL = ["#FDDBC7", "#DC6E57", "#A21328"]
# ColorBrewer Blues 7 级（**历史遗留：当前无脚本引用**，保留以便回退对照）
BLUES7 = ["#F7FBFF", "#E4EEF7", "#9ECAE1", "#6BAED6", "#4292C6", "#08519C", "#08306B"]
# 覆盖式概率图（叠在浅灰模板上）：低端不能近白，否则与浅灰模板糊在一起。
# 2026-09-21 第二轮：随 SEQUENTIAL 改为暖臂，低端取 #F1D4AA（浅米黄，仍比模板深）。
# 校验：相邻级 正常 10.7 ／ 红 8.8 ／ 绿 10.1 ／ 蓝 6.3；L* 86→81→44→34→27 单调 ✓
SEQUENTIAL_OVERLAY = ["#FDDBC7", "#F19E7D", "#CF5246", "#AB162A", "#67001F"]
ACCENT_RED = "#BA2832"    # v5：RdBu_r t=0.875；WCAG 白底对比 7.54（旧 #B2182B 为 6.87）
GREY = "#999999"
BLACK = "#222222"

# scipilot 技能推荐的色盲安全备选（Okabe-Ito）
OKABE_ITO = ["#0072B2", "#D55E00", "#009E73", "#E69F00", "#CC79A7", "#56B4E9", "#000000"]

CMAP_SEQ = LinearSegmentedColormap.from_list('afs_seq', SEQUENTIAL)
CMAP_BLUES7 = LinearSegmentedColormap.from_list('afs_blues7', BLUES7)
CMAP_OVERLAY = LinearSegmentedColormap.from_list('afs_overlay', SEQUENTIAL_OVERLAY)
CMAP_DIV = LinearSegmentedColormap.from_list('afs_div', DIVERGING)
CMAP_QUAL = LinearSegmentedColormap.from_list('afs_qual', CATEGORICAL_EXTENDED)

# ---------------------------------------------------------------- 单极强度色标
# 用于「只有一个方向、但要表达强弱」的量（如 |BSR|、富集强度）。
#
# **取 DIVERGING 正向臂的一段**，而不是另起一套红色 —— 这样同一张图里
# 「面板 a 的正向红」与「面板 b 的强度红」是**严格同一批颜色**，不会看着"不齐"。
# （2026-09-17 修正：原用 ColorBrewer RdBu 的 #FDDBC7→#EF8A62→#B2182B，
#  其中 #EF8A62 是**橙调**、并不在 DIVERGING 这条臂上，导致 b 面板偏橙、
#  与 a 面板的粉红明显不同调。）
#
# 低端取 t=0.58（#ECD3D6）而非臂的中点 t=0.50（#F7F6F6）：脑表面 inflated 渲染的
# 高光接近纯白（实测灰底像素里出现 255），近白色斑块会与背景糊在一起；
# t=0.58 仍带 25 级彩度、肉眼可辨，且本来就属于 a 面板的用色范围。
SEQ_RED_T = (0.50, 0.625, 0.75, 0.875, 1.00)
CMAP_SEQ_RED = LinearSegmentedColormap.from_list(
    'afs_seq_red', [CMAP_DIV(t) for t in SEQ_RED_T])


def f3(v, nd=3, sign=False):
    """图上数值标签的**统一格式化**：Decimal 半进位（ROUND_HALF_UP）。

    为什么必须用它而不是 '%.3f'：稿件采用半进位，而 '%.3f' 作用于二进制浮点，
    在 X.XXX5 这类十进制边界上会掉到相邻偶数或受浮点表示偏移影响，
    导致**图上标签与稿件差 1 个末位**。本稿已三次栽在此处：
      · Fig. 7 d 面板 OXPHOS：源值 0.3655 → '%.3f' 得 +0.365、稿件 +0.366
      · α4β2 / μ-阿片偏相关（现由正文 Table 3 承载）：全精度源值 0.3014746 / 0.0094504 →
        '%.3f' 得 0.301 / 0.009（与现稿一致；旧稿 0.302 / 0.010 系中间量双重舍入，2026-10-10 修正）
    凡在图上打印数值标签，一律走本函数。
    """
    from decimal import Decimal, ROUND_HALF_UP
    d = Decimal(str(float(v))).quantize(Decimal(1).scaleb(-nd), rounding=ROUND_HALF_UP)
    s = ('%%.%df' % nd) % d
    if sign and d >= 0:
        s = '+' + s
    return s


def apply_style(base_size=7.0):
    """全局样式（每个出图脚本第一行调用）"""
    plt.rcParams.update({
        'font.family': 'sans-serif',
        'font.sans-serif': ['Arial', 'Helvetica', 'Liberation Sans', 'DejaVu Sans'],
        'font.size': base_size,
        'axes.titlesize': base_size + 1,
        'axes.labelsize': base_size,
        'xtick.labelsize': base_size - 0.5,
        'ytick.labelsize': base_size - 0.5,
        'legend.fontsize': base_size - 0.5,
        'axes.linewidth': 0.6,
        'axes.spines.top': False,
        'axes.spines.right': False,
        'xtick.major.width': 0.6,
        'ytick.major.width': 0.6,
        'xtick.major.size': 2.0,
        'ytick.major.size': 2.0,
        'grid.linewidth': 0.4,
        'lines.linewidth': 0.8,
        'axes.grid': False,
        'figure.dpi': 300,
        'savefig.dpi': 300,
        'savefig.bbox': 'tight',
        'savefig.pad_inches': 0.01,
        'pdf.fonttype': 42,        # 文字不转曲 → Illustrator 可编辑
        'ps.fonttype': 42,
        'svg.fonttype': 'none',
    })


def panel_labels(axes, labels=None, style='lower', dx=-0.06, dy=1.06):
    """子图标签。本稿图注沿用 a/b/c（小写）以与既有图注一致
    —— academic-figure-skill 对 Elsevier 给的是 (A)(B)(C)，此处**有意偏离**并已记录。"""
    axes = list(axes)
    if labels is None:
        labels = [chr(ord('a') + i) for i in range(len(axes))]
    for ax, lab in zip(axes, labels):
        txt = lab if style == 'lower' else '(%s)' % lab
        ax.text(dx, dy, txt, transform=ax.transAxes, fontsize=8, fontweight='bold',
                va='top', ha='right', color=BLACK)
    return labels


def save_all(fig, basename, outdir='.', dpi=300, formats=('pdf', 'svg', 'eps', 'png'),
             vec_dpi=600, eps_dpi=300, bbox_inches='tight'):
    """五格式导出（矢量优先）；文件被占用时改存 _new 并告警。

    vec_dpi    只控制**嵌入矢量文件中的位图**（rasterized=True 的脑图/热图）分辨率，
               文字与线条仍是矢量、可编辑。
    eps_dpi    EPS 内嵌位图分辨率（EPS 的位图为无压缩编码，60 张脑切片时 600 dpi
               可致文件逾 50 MB → 默认降到 300 dpi）。
    bbox_inches='tight' 会按内容裁剪 → 成品物理尺寸小于版面标称值；投稿要求精确
               版面尺寸时传 None（此时版面 = figsize）。
    """
    written = []
    for ext in formats:
        path = os.path.join(outdir, '%s.%s' % (basename, ext))
        d = dpi if ext == 'png' else (eps_dpi if ext == 'eps' else vec_dpi)
        try:
            fig.savefig(path, dpi=d, bbox_inches=bbox_inches)
            written.append(path)
            print('  已写出：%s' % path)
        except PermissionError:
            alt = path.replace('.%s' % ext, '_new.%s' % ext)
            fig.savefig(alt, dpi=d, bbox_inches=bbox_inches)
            written.append(alt)
            print('  [!] %s 被占用 → 已改存 %s' % (os.path.basename(path), os.path.basename(alt)))
    return written
