# -*- coding: utf-8 -*-
"""verify_manuscript_claims.py —— 逐条核对稿件结论

给独立验证者用。**只读本仓库 `derived_data/` 里的 CSV，不依赖作者本机任何路径**，
也不联网。每一条"稿件声称"都对应一段可复算的计算，结果以 PASS / FAIL / 需重跑 呈现。

用法
    python verify_manuscript_claims.py
    （脚本自动定位 derived_data/：支持「提交包根目录」与「审阅包
      08_验证与审计/01_数值自检/」两种布局，任意工作目录均可运行。）

依赖：pandas, numpy, scipy（见 requirements.txt）

关于小数比对
    稿件按**四舍五入到三位**（半进位）报告。本脚本把重算值用同样的规则取整后再比，
    避免 0.3826 与 0.383 被误判为不一致。容差为最后一位的半个单位。

剔除的条目（2026-09-19 / 2026-09-22 两次稿件瘦身）
    只核对**当前稿件仍然声称的结论**。稿件正文与补充材料已删除的分析不再列入，
    否则会出现"脚本在核对一个稿件已不再声称的结论"：
      · corr(Lx, Ly) 的 patient-resampling bootstrap 95% CI（+0.287 ~ +0.485）——
        正文与 Fig. 4 图注均已删除；
      · 控制 21 张图谱第一主成分（PC1）后的偏相关（+0.384 / +0.346 / −0.256）——
        正文与补充材料均已删除；
      · 最大病灶个体的杠杆点检验——2026-09-22 随 §3.3(iii) 与 Fig. 4 的
        leverage 口径一并删除（相关台账行已于 2026-10-10 一并移除）；
      · top-300 基因列表的 Jaccard 重叠（双侧 0.19–0.28、左半球 0.44–0.65）——
        随补充图 Fig. S1 一并删除（2026-09-22），其数据文件已于 2026-10-10
        随 `legacy_removed_20260923/` 目录归档至仓库之外（移动、未删除）。
    另有少数条目的标签里注明了"正文未报"或"补充材料表 S3"——它们是正文未直接给出数值、
    但仍由图表或补充材料支撑的校验（如无协变量 corr、Fig. 5c 的点位、表 S3 的家族比较）。
"""
import io
import os
import sys
from decimal import Decimal, ROUND_HALF_UP

import numpy as np
import pandas as pd
from scipy import stats

try:
    sys.stdout.reconfigure(encoding='utf-8')
except Exception:
    pass

_HERE = os.path.dirname(os.path.abspath(__file__))


def _resolve_derived():
    """定位 derived_data/ 目录（兼容三种布局）：
      ① 提交包根（脚本旁）；② 审阅包布局（自 08_验证与审计/01_数值自检/ 向上查找
      `06_数据/01_结果数据_derived_data/derived_data/`）；③ 任意工作目录（两者兼搜）。
    """
    cands = []
    for base in (_HERE, os.getcwd()):
        cur = base
        for _ in range(4):
            cands.append(os.path.join(cur, 'derived_data'))
            cands.append(os.path.join(cur, '06_数据', '01_结果数据_derived_data', 'derived_data'))
            cur = os.path.dirname(cur)
    for c in cands:
        if os.path.isdir(c):
            return c + os.sep
    raise FileNotFoundError(
        '未能定位 derived_data/：请在提交包根目录运行，或确认审阅包内 '
        '06_数据/01_结果数据_derived_data/derived_data/ 存在。')


D = _resolve_derived()
results = []          # (章节, 声称, 重算值, 结论, 备注)


def r3(x):
    """半进位到三位小数（与稿件一致）"""
    return float(Decimal(str(float(x))).quantize(Decimal('0.001'), rounding=ROUND_HALF_UP))


def r2(x):
    return float(Decimal(str(float(x))).quantize(Decimal('0.01'), rounding=ROUND_HALF_UP))


def r1(x):
    """百分比类数值稿件只报一位小数（41.7% / 62.7%）"""
    return float(Decimal(str(float(x))).quantize(Decimal('0.1'), rounding=ROUND_HALF_UP))


def check(section, claim, got, expected, note=''):
    ok = (got is not None) and abs(float(got) - float(expected)) < 5e-4
    results.append((section, claim, got, expected, 'PASS' if ok else 'FAIL', note))
    return ok


def load(name):
    return pd.read_csv(D + name)


def load_ix(name):
    """首列为行索引（相关矩阵、多种子矩阵都是这种格式）"""
    return pd.read_csv(D + name, index_col=0)


def main():
    # ================= §3.2 行为 PLS =================
    # ⚠️ 2026-09-19 稿件瘦身后（正文 −20%）：corr(Lx, Ly) 的 patient-resampling bootstrap
    #    95% CI（+0.287 ~ +0.485）已从稿件正文与 Fig. 4 图注中删除，故不再作为"稿件声称"核对。
    #    本 CSV 仍随包提供——下面 corr 观测值（+0.317）就取自同一文件的 obs_corr 列。
    boot = load('pls_LV1_bootstrap_CI.csv').iloc[0]
    check('3.2', 'corr(Lx, Ly) = +0.317', r3(boot['obs_corr']), 0.317)

    sal = load('pls_LV1_brain_loading_cortex100.csv')
    check('3.2', '载荷为正的区域 = 80/100', int((sal['salience'] > 0).sum()), 80)
    check('3.2', '最小载荷 = −0.062', r3(sal['salience'].min()), -0.062)
    check('3.2', '|BSR| > 2 的稳健区域 = 24', int((sal['BSR'].abs() > 2).sum()), 24)

    # 病灶体积 × FSS-9 的 Spearman（正文 §3.2(ii) 报告 ρ = +0.024, p = 0.791）
    lv = load('lesion_volume_descriptives.csv').set_index('variable')
    check('3.2', '病灶体积 × FSS-9 Spearman rho = +0.024', r3(lv.loc['FSS9总分', 'rho']), 0.024)
    check('3.2', '病灶体积 × FSS-9 Spearman p = 0.791', r3(lv.loc['FSS9总分', 'p']), 0.791)
    # 网络分布
    net = sal.loc[sal['BSR'].abs() > 2, 'region'].str.extract(r'^[LR]H_([A-Za-z]+?)_')[0]
    net = net.replace({'SomMot': 'SomMot', 'Default': 'Default', 'Cont': 'Cont',
                       'DorsAttn': 'DorsAttn', 'SalVentAttn': 'SalVentAttn'})
    for name, exp in (('Default', 8), ('DorsAttn', 5), ('Cont', 5), ('SomMot', 5)):
        check('3.2', '稳健区域中 %s 网络 = %d' % (name, exp), int((net == name).sum()), exp)

    sens = load('sensitivity_specs.csv')
    for spec_key, exp, label in (
            ('无协变量', 0.017, '无协变量 置换 p = 0.017'),
            ('+log10 体积（仅 Y 残差）', 0.035, '+log10 体积 置换 p = 0.035'),
            ('FSS-9 总分单调秩变换', 0.044, 'FSS-9 单调秩变换 置换 p = 0.044')):
        row = sens.loc[sens['spec'] == spec_key, 'perm_p']
        check('3.2', label, r3(row.iloc[0]) if len(row) else None, exp)

    for spec_key, exp, label in (
            # 稿件正文只报该口径的置换 p（= 0.017）；corr 值原以旧 Fig. 4b 面板呈现（该面板已删，值保留供复算）
            ('无协变量', 0.322, '无协变量 corr = +0.322（正文未报；旧 Fig. 4b 已删，值供复算）'),
            ('+log10 体积（仅 Y 残差）', 0.314, '+log10 体积 corr = +0.314'),
            ('FSS-9 总分单调秩变换', 0.326, 'FSS-9 单调秩变换 corr = +0.326')):
        row = sens.loc[sens['spec'] == spec_key, 'corr']
        check('3.2', label, r3(row.iloc[0]) if len(row) else None, exp)

    # ⚠️ 2026-09-22 稿件修订：**最大病灶个体的杠杆点检验已整段删除**——原 §3.3 的 (iii) 段
    #    与 Fig. 4 的两行 leverage 口径均不再报告，故这 5 条一并移出核对。
    #    sensitivity_specs.csv 里 group = 'leverage' 的两行已于 2026-10-10 一并移除
    #    （该表其余口径仍在使用）。
    #    注：|BSR|>2 区域数 = 24 本身仍由上面"|BSR| > 2 的稳健区域 = 24"一条核对。


    # ---- LV1 置换零分布（§3.2 的置换 p = 0.023 与 94.0% 方差解释）----
    #      注：该面板已于 2026-09-18 从 Fig. 4 移除（数值仍在正文与 Fig. 3 图注中报告），数据保留供核对
    sp = load('pls_singval_spectrum.csv')
    check('3.2', 'LV1 解释交叉协方差 = 94.0%',
          r1(sp.loc[sp['LV'] == 1, 'varexp_pct'].iloc[0]), 94.0)
    check('3.2', 'LV2-LV9 合计 = 6.0%',
          r1(sp.loc[sp['LV'] > 1, 'varexp_pct'].sum()), 6.0)
    nsum = load('pls_perm_null_summary.csv').set_index('LV')
    obs_sv = float(nsum.loc[1, 'observed_singval'])
    check('3.2', '观测 LV1 奇异值 = 3.82（内部量，正文未报）', r2(obs_sv), 3.82)
    check('3.2', '空分布中不小于观测值的置换数 = 116',
          int(nsum.loc[1, 'n_exceed']), 116)
    null1 = load('pls_perm_null_singvals.csv')['LV1'].values
    nperm = len(null1)
    p_null = (int((null1 > obs_sv).sum()) + 1) / (nperm + 1)
    check('3.2', '由零分布复算的主分析置换 p = 0.023', r3(p_null), 0.023,
          note='%d 次置换' % nperm)
    check('3.2', '摘要表 perm_p 与零分布矩阵复算值一致',
          r3(p_null), r3(float(nsum.loc[1, 'perm_p'])))

    # ================= §3.3 供体稳健性（LODO）=================
    # 供体编号（9861 等）必须按字符串读，否则 pandas 会推断成整数、索引取不到
    lc = pd.read_csv(D + 'lodo_consistency_by_donor.csv',
                     dtype={'donor': str}).set_index('donor')
    bil = lc[lc['hemisphere'] == 'bilateral']['pearson_r'].astype(float)
    lef = lc[lc['hemisphere'] == 'left-only']['pearson_r'].astype(float)
    check('3.3', '其余 4 例供体 LODO r 下限 = 0.959', r3(lef.min()), 0.959)
    check('3.3', '其余 4 例供体 LODO r 上限 = 0.966', r3(lef.max()), 0.966)
    check('3.3', '供体 9861（双侧）LODO r = +0.783',
          r3(lc.loc['9861', 'pearson_r']), 0.783)
    check('3.3', '供体 10021（双侧）LODO r = +0.829',
          r3(lc.loc['10021', 'pearson_r']), 0.829)

    # 公共有效区（90 个）约束下重算
    cr = pd.read_csv(D + 'lodo_common90_r_by_donor.csv', dtype={'donor': str})
    check('3.3', '约束区数（7 组矩阵公共有效区）= 90',
          int(cr['n_common_regions'].iloc[0]), 90)
    cri = cr.set_index('donor')['r_common_regions'].astype(float)
    check('3.3', '公共 90 区重算：9861 r = +0.744', r3(cri.loc['9861']), 0.744)
    check('3.3', '公共 90 区重算：10021 r = +0.839', r3(cri.loc['10021']), 0.839)

    # ⚠️ 2026-09-22 稿件修订：**补充图 Fig. S1（基因排序 LODO 图）整图删除**，
    #    top-300 基因列表的 Jaccard 重叠随之不再是稿件的声称，原 3 条核对已移出。
    #    数据文件 lodo_top300_jaccard_by_donor.csv 已于 2026-10-10 随
    #    `legacy_removed_20260923/` 目录归档至仓库之外（移动、未删除）。
    #    仍在核对的 §3.3 数值：各供体 LODO r（lodo_consistency_by_donor.csv，上）
    #    与公共 90 区重算 r（lodo_common90_r_by_donor.csv，下）——二者均为正文句内报告。

    # ================= §3.4 细胞类型层 =================
    ct = load('celltype_enrichment_Dataset5.csv').set_index('cell_type')
    for ct_name, rank, label in (('Micro', 6242, '小胶质 中位秩 = 6,242'),
                                 ('Oligodendrocytes', 5232, '少突胶质 中位秩 = 5,232')):
        check('3.4', label, int(ct.loc[ct_name, 'median_rank']), rank)
    check('3.4', 'OPC BH-adjusted p = 0.993', r3(ct.loc['OPC', 'fdr_bh']), 0.993)

    sp = load('celltype_spatial_spin_Dataset5.csv')
    sp_sal = sp[sp['mode'] == 'salience'].set_index('feature')
    sp_dis = sp[sp['mode'] == 'discon_mean'].set_index('feature')
    check('3.4', '小胶质 空间 r = +0.290', r3(sp_sal.loc['Cell_Micro', 'r']), 0.290)
    check('3.4', '小胶质 spin p = 0.044', r3(sp_sal.loc['Cell_Micro', 'P_SPIN']), 0.044)
    check('3.4', '少突胶质 空间 r = +0.212', r3(sp_sal.loc['Cell_Oligodendrocytes', 'r']), 0.212)
    check('3.4', '少突胶质 spin p = 0.020', r3(sp_sal.loc['Cell_Oligodendrocytes', 'P_SPIN']), 0.020)
    check('3.4', '小胶质 vs 断连负荷 r = −0.212',
          r3(sp_dis.loc['Cell_Micro', 'r']), -0.212)
    check('3.4', '小胶质 vs 断连负荷 p = 0.176',
          r3(sp_dis.loc['Cell_Micro', 'P_SPIN']), 0.176)

    # 2026-10-10 追加：偏相关列与抑制性神经元行（原缺口见 celltype_specificity_contrast_说明.md）
    ccp = load('celltype_specificity_contrast.csv').set_index('cell_type')
    for ct_name, rv, pv in (('Micro', 0.258, 0.051),
                            ('Oligodendrocytes', 0.207, 0.019)):
        check('3.4', '%s 控负荷后偏 r = %+.3f（表 S1）' % (ct_name, rv),
              r3(ccp.loc[ct_name, 'r_partial']), rv)
        check('3.4', '%s 控负荷后偏 r spin p = %.3f（表 S1）' % (ct_name, pv),
              r3(ccp.loc[ct_name, 'p_partial']), pv)
    check('3.4', '抑制性神经元 空间 r = −0.120（表 S1）',
          r3(ccp.loc['Inhibitory_Neurons', 'r_salience']), -0.120)
    check('3.4', '抑制性神经元 spin p = 0.034（表 S1）',
          r3(ccp.loc['Inhibitory_Neurons', 'p_salience']), 0.034)
    check('3.4', '抑制性神经元 控负荷后偏 r = −0.124（表 S1）',
          r3(ccp.loc['Inhibitory_Neurons', 'r_partial']), -0.124)
    check('3.4', '抑制性神经元 控负荷后偏 r spin p = 0.028（表 S1）',
          r3(ccp.loc['Inhibitory_Neurons', 'p_partial']), 0.028)

    # §3.6 受体层由 ms 行开始
    ms = load_ix('multiseed_celltype_Dataset5_p.csv')      # 20 seeds × 7 cell types
    check('3.4', '小胶质在 20/20 个种子中 p<0.05', int((ms['Micro'] < 0.05).sum()), 20)
    check('3.4', '少突胶质在 20/20 个种子中 p<0.05',
          int((ms['Oligodendrocytes'] < 0.05).sum()), 20)
    check('3.4', 'OPC 在 0/20 个种子中 p<0.05', int((ms['OPC'] < 0.05).sum()), 0)

    # ================= §3.5 通路层 =================
    go = load('gsea_GO_w1.csv')
    kegg = load('gsea_KEGG_w1.csv')
    qcol_go = [c for c in go.columns if c.lower() in ('fdr q-val', 'qvalue')][0]
    qcol_kg = [c for c in kegg.columns if c.lower() in ('fdr q-val', 'qvalue')][0]
    check('3.5', 'GO 检验条目数 = 5,301', len(go), 5301)
    check('3.5', 'GO 显著条目（q<0.05）= 134', int((go[qcol_go] < 0.05).sum()), 134)
    check('3.5', 'KEGG 检验条目数 = 312', len(kegg), 312)
    check('3.5', 'KEGG 显著条目（q<0.05）= 50', int((kegg[qcol_kg] < 0.05).sum()), 50)

    core = load('core_pathway_results.csv')
    op = core.loc[core['pathway'].str.contains('Oxidative phosphorylation', case=False, na=False), 'NES']
    # 稿件此处报 2 位小数（2.79），用 r2 比较；全精度值为 2.7860799
    check('3.5', '氧化磷酸化 NES = 2.79（稿件报 2 位）', r2(op.iloc[0]) if len(op) else None, 2.79)

    pspin = load('gsea_pathway_spin_survival.csv')
    check('3.5', '通路层 spin 检验条数 = 184', len(pspin), 184)
    # 2026-10-10 同步：10-09 通路层口径更新后，表型层（同向单侧）存活列更名为 p_spin_one
    # （127 存活 / 57 淘汰，与正文 §3.5(ii) 及 Fig. 5c 图注一致；旧期望 143/41 系旧版数据）。
    check('3.5', '通路层 spin 存活（p<0.05）= 127', int((pspin['p_spin_one'] < 0.05).sum()), 127)
    check('3.5', '通路层 spin 未存活 = 57', int((pspin['p_spin_one'] >= 0.05).sum()), 57)

    # 正文不再点名具体的"未存活"条目（只报存活 143 / 淘汰 41），
    # 但这两个点仍在 Fig. 5c 的 184 点之内 → 保留为图表级校验
    for key, exp, label in (
            ('SRP-dependent cotranslational', 0.065,
             'SRP 协同翻译 spin p = 0.065（正文未点名；Fig. 5c 含该点）'),
            ('cytochrome c to oxygen', 0.105,
             '细胞色素 c 至氧的电子传递 spin p = 0.105（正文未点名；Fig. 5c 含该点）')):
        row = core.loc[core['pathway'].str.contains(key, case=False, na=False), 'p_spin']
        check('3.5', label, r3(row.iloc[0]) if len(row) else None, exp)

    # ================= §3.6 受体层 =================
    rec = load('receptor_mapping_21atlases.csv').set_index('target')
    for tgt, name, rv, pv, fdr in (
            ('a4b2-nACh', 'α4β2', 0.383, 0.060, 0.489),
            ('mu-opioid', 'μ-阿片', 0.274, 0.168, 0.489),
            ('GABA-A', 'GABA-A', -0.255, 0.140, 0.489)):
        row = rec.loc[tgt]
        check('3.6', '%s r = %+0.3f（Table 3）' % (name, rv), r3(row['r']), rv)
        check('3.6', '%s spin p = %.3f（Table 3）' % (name, pv), r3(row['spin_p']), pv)
        check('3.6', '%s BH 校正 p = %.3f（Table 3）' % (name, fdr), r3(row['fdr']), fdr)
    check('3.6', '21 次检验 BH 校正后显著图谱数 = 0', int((rec['fdr'] < 0.05).sum()), 0)
    check('3.6', '21 次检验 BH 校正后不显著图谱数 = 21', int((rec['fdr'] >= 0.05).sum()), 21)

    # 图谱两两相关：**必须从 100 区原始值全精度重算**。
    # 陷阱：receptor_21x21_correlation.csv 里的值是已舍入到 4 位的中间量；
    # 若在它上面再做 3 位半进位，−0.2484959… 会先被存成 −0.2485、再进位成 −0.249，
    # 与稿件的 −0.248 造成**假 FAIL**。舍入只能在全精度值上做一次。
    pv = load('receptor_parcel_values_21x100.csv')
    X = pv.iloc[:, 5:].values.astype(float)      # 前 5 列为元数据，其后 100 个皮层区
    C = np.corrcoef(X)
    tgt_list = list(pv['target'])

    def pcorr(a, b):
        return C[tgt_list.index(a), tgt_list.index(b)]

    check('3.6', 'α4β2 与 μ-阿片 r = +0.515（由 100 区重算）',
          r3(pcorr('a4b2-nACh', 'mu-opioid')), 0.515)
    check('3.6', 'μ-阿片 与 GABA-A r = −0.540（由 100 区重算）',
          r3(pcorr('mu-opioid', 'GABA-A')), -0.540)
    check('3.6', 'α4β2 与 GABA-A r = −0.248（由 100 区重算，临界值）',
          r3(pcorr('a4b2-nACh', 'GABA-A')), -0.248)

    corr21 = load_ix('receptor_21x21_correlation.csv')
    maxdiff = float(np.abs(C - corr21.values).max())
    results.append(('3.6', '包内 21×21 矩阵可由 100 区数值复算（最大逐格差 < 1e-4）',
                    '%.2e' % maxdiff, '<1e-4', 'PASS' if maxdiff < 1e-4 else 'FAIL', ''))

    # PC1 方差解释比（可由相关矩阵的特征值直接得到）
    def pc1_var(rows_idx):
        w = np.linalg.eigvalsh(np.corrcoef(X[rows_idx]))[::-1]
        return w[0] / w.sum() * 100.0

    idx3 = [tgt_list.index(t) for t in ('a4b2-nACh', 'mu-opioid', 'GABA-A')]
    check("3.6", "21 张图谱整体 PC1 方差解释 = 41.7%", r1(pc1_var(range(len(tgt_list)))), 41.7)
    check('3.6', '三张显著图谱 PC1 方差解释 = 62.7%', r1(pc1_var(idx3)), 62.7)

    # ⚠️ 2026-09-19 稿件瘦身后：控制 21 张图谱第一主成分（PC1）后的偏相关已**整条删除**
    #    ——正文与补充材料均已不再报告，
    #    故这 3 条声称与其数据文件 receptor_partial_control_PC1.csv 一并移出本包。
    #    注意：上面两行"PC1 方差解释比"仍保留——它们描述的是**图谱间共线程度**，稿件仍在报告。

    # 2026-09-22 稿件修订（+ 2026-10-10 受体层口径定稿）：原 Fig. 7（PET 受体图谱图）整图删除，
    # 其数值由**正文 Table 3**承载——偏相关三条即 Table 3 的 α4β2 +0.301 / μ-阿片 +0.009 /
    # GABA-A −0.151（2026-10-10 由全精度 +0.3014746 / +0.0094504 / −0.1511337 一次性半进位修正，
    # 取代旧值 +0.302 / +0.010 / −0.151——旧系 4 位中间量双重舍入所致）。
    par = load('receptor_partial_correlation.csv').set_index('target')
    for tgt, name, rv, pv in (
            ('a4b2-nACh', 'α4β2', 0.301, 0.074),
            ('mu-opioid', 'μ-阿片', 0.009, 0.469),
            ('GABA-A', 'GABA-A', -0.151, 0.268)):
        check('3.6', '控制彼此后 %s 偏 r = %+0.3f（Table 3）' % (name, rv), r3(par.loc[tgt, 'partial_r']), rv)
        check('3.6', '控制彼此后 %s 偏 r spin p = %.3f（Table 3）' % (name, pv), r3(par.loc[tgt, 'partial_p_spin']), pv)

    ch = load('receptor_cholinergic_family.csv').set_index('tracer_target')
    check('3.6', 'VAChT r = +0.058（Table 3 其余胆碱能图谱）', r3(ch.loc['VAChT', 'r']), 0.058)
    check('3.6', 'M1 r = −0.048（Table 3 其余胆碱能图谱）', r3(ch.loc['M1 mAChR', 'r']), -0.048)
    check('3.6', '控另两张胆碱能后 α4β2 偏 r = +0.403（Table 3）',
          r3(ch.loc['a4b2-nACh', 'partial_r_ctrl_other_two_cholinergic']), 0.403)
    check('3.6', '控另两张胆碱能后 α4β2 偏 r spin p = 0.041（Table 3；表内唯一 p<0.05）',
          r3(ch.loc['a4b2-nACh', 'partial_p_spin']), 0.041)
    check('3.6', '控另两张胆碱能后 VAChT 偏 r spin p = 0.210（Table 3）',
          r3(ch.loc['VAChT', 'partial_p_spin']), 0.210)
    check('3.6', '控另两张胆碱能后 M1 偏 r spin p = 0.302（Table 3）',
          r3(ch.loc['M1 mAChR', 'partial_p_spin']), 0.302)

    chp = load('receptor_cholinergic_pairwise.csv')
    check('3.6', '胆碱能两两相关文件行数 = 3', len(chp), 3)

    msp = load_ix('multiseed_receptor_spin_p_matrix.csv')   # 21 atlases × 20 seeds
    rn = lambda t: rec.loc[t, 'receptor']                   # target → 图谱文件名
    check('3.6', 'α4β2 在部分种子中 p<0.05（2/20）',
          int((msp.loc[rn('a4b2-nACh')].values < 0.05).sum()), 2)
    for tgt in ('mu-opioid', 'GABA-A'):
        check('3.6', '%s 在 20 个种子中均不显著（0/20）' % tgt,
              int((msp.loc[rn(tgt)].values < 0.05).sum()), 0)
    others = [t for t in rec.index if t != 'a4b2-nACh']
    n_any = sum(1 for t in others if (msp.loc[rn(t)].values < 0.05).sum() > 0)
    check('3.6', '其余 20 张图谱均无任何种子显著（均 0/20）', n_any, 0)
    check('3.6', 'α4β2 种子间 p 中位数 = 0.057',
          r3(np.median(msp.loc[rn('a4b2-nACh')].values)), 0.057)
    check('3.6', 'α4β2 种子间 p 最大值 = 0.065',
          r3(np.max(msp.loc[rn('a4b2-nACh')].values)), 0.065)

    # ================= 讨论：效能量化 =================
    # 稿件：以 corr(Lx, Ly) = +0.317 与双侧 α = 0.05 推算，n = 124 的检验效能约 0.95。
    # 这里按 Fisher z 近似解析式复算，只需 n 与 r，不需要任何患者级数据。
    from scipy.stats import norm
    r_obs = float(sens.loc[sens['spec'] == '主分析（Y~年龄+性别）', 'corr'].iloc[0])
    n_obs = 124
    zc = np.arctanh(r_obs)
    se = 1.0 / np.sqrt(n_obs - 3)
    crit = norm.ppf(1 - 0.05 / 2)
    power = norm.cdf(zc / se - crit) + norm.cdf(-zc / se - crit)
    check('讨论', '效能量化：n = 124、corr = +0.317 时 power ≈ 0.95',
          r2(power), 0.95, note='Fisher z 解析式，仅需 n 与 r')

    # ================= 输出 =================
    w = max(len(x[1]) for x in results) + 2
    print('=' * (w + 42))
    print('稿件结论逐条核对（重算值仅来自 derived_data/）')
    print('=' * (w + 42))
    print('%-6s %-*s %12s %12s  %s' % ('章节', w, '稿件声称', '重算值', '期望值', '结论'))
    print('-' * (w + 42))
    for sec, claim, got, exp, verdict, note in results:
        if got is None:
            gs = '—'
        elif isinstance(got, (int, float)):
            gs = '%g' % got
        else:
            gs = str(got)
        print('%-6s %-*s %12s %12s  %s%s' % (sec, w, claim[:w], gs, exp, verdict,
                                             ('  ← ' + note) if note else ''))
    npass = sum(1 for x in results if x[4] == 'PASS')
    nfail = sum(1 for x in results if x[4] == 'FAIL')
    print('-' * (w + 42))
    print('合计 %d 条：PASS %d ｜ FAIL %d' % (len(results), npass, nfail))
    print()
    print('说明：FAIL 请先确认是否为本仓库数据版本与稿件版本不一致；')
    print('      临界值务必回到全精度数据重算，不要在已舍入的中间表上再舍一次')
    print('      （receptor_21x21_correlation.csv 已舍入到 4 位，直接再用会造成假 FAIL）。')
    print('      需要重跑上游才能得到的项目（如 LV1 方差解释、HPA 过滤 GSEA 计数）')
    print('      见 README_中文.md 第 7 步。')
    return nfail


if __name__ == '__main__':
    sys.exit(0 if main() == 0 else 1)
