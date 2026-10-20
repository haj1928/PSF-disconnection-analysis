# -*- coding: utf-8 -*-
"""repro_smoke_test.py —— 包内快速健康检查（约 5 秒）

用途：换环境、或对包内容有任何改动后，先跑本脚本快速确认「关键数据未串、未损」。
本脚本只读 `derived_data/`，不联网、不依赖作者本地路径。
完整核对（稿件全部结论，96 条）请运行 `verify_manuscript_claims.py`。

检查项：
  1  关键文件齐全
  2  数据表形状（21×100 区域值 / 184 通路 / 7 细胞类型 / 15,633 基因）
  3  PLS 锚点：corr = +0.317、perm p = 0.023
  4  受体锚点：α4β2 spin p = 0.060；BH 最小 = 0.4888（0/21 显著）；偏相关全精度 0.3014746
  5  通路锚点：表型层存活 127/184
  6  细胞层：PRAND p 全部落在 (1+k)/(N+1) 网格上
  7  21×21 相关矩阵可由 100 区全精度值复算（max|diff| < 1e-4；抽查 α4β2–μ = +0.515）

用法：python repro_smoke_test.py
"""
import os
import sys

import numpy as np
import pandas as pd

try:
    sys.stdout.reconfigure(encoding="utf-8")
except Exception:
    pass

HERE = os.path.dirname(os.path.abspath(__file__))
D = os.path.join(HERE, "derived_data") + os.sep

results = []


def chk(label, ok, detail=""):
    results.append((label, bool(ok), detail))
    print("[%s] %s%s" % ("PASS" if ok else "FAIL", label, ("  ← " + detail) if detail else ""))


def load(name):
    return pd.read_csv(D + name)


print("=== repro_smoke_test.py — 快速健康检查 ===")
print()

# ---------------- 1 关键文件 ----------------
must = [
    "pls_LV1_bootstrap_CI.csv", "sensitivity_specs.csv",
    "receptor_mapping_21atlases.csv", "receptor_partial_correlation.csv",
    "receptor_parcel_values_21x100.csv", "receptor_21x21_correlation.csv",
    "gsea_pathway_spin_survival.csv", "celltype_enrichment_Dataset5.csv",
    "gene_ranking_w1_cortex100.csv", "gsea_KEGG_w1.csv",
]
missing = [f for f in must if not os.path.exists(D + f)]
chk("关键文件齐全（%d/%d）" % (len(must) - len(missing), len(must)), not missing,
    ("缺: " + ", ".join(missing)) if missing else "")

# ---------------- 2 形状 ----------------
pv = load("receptor_parcel_values_21x100.csv")
chk("21 张受体图谱 × 100 区（区域值表形状）",
    pv.shape == (21, 105), "shape=%s" % (pv.shape,))
chk("100 区区域值无 NaN/Inf",
    np.isfinite(pv.iloc[:, 5:].values.astype(float)).all())

g = load("gsea_pathway_spin_survival.csv")
chk("通路层 184 条（GO+KEGG 显著集）", len(g) == 184, "n=%d" % len(g))

ct = load("celltype_enrichment_Dataset5.csv")
chk("细胞层 7 类", len(ct) == 7, "n=%d" % len(ct))

rk = load("gene_ranking_w1_cortex100.csv")
chk("基因排序表 15,633 行", len(rk) == 15633, "n=%d" % len(rk))

# ---------------- 3 PLS 锚点 ----------------
boot = load("pls_LV1_bootstrap_CI.csv")
obs = float(boot["obs_corr"].iloc[0])
chk("PLS：obs_corr = 0.3171（圆整 0.317）", round(obs, 4) == 0.3171, "got=%.4f" % obs)

sens = load("sensitivity_specs.csv")
row = sens[sens["spec"] == "主分析（Y~年龄+性别）"].iloc[0]
chk("主分析 corr = +0.317 / perm p = 0.023（敏感性表）",
    round(float(row["corr"]), 3) == 0.317 and round(float(row["perm_p"]), 3) == 0.023,
    "corr=%.6f p=%.6f" % (float(row["corr"]), float(row["perm_p"])))

# ---------------- 4 受体锚点 ----------------
rec = load("receptor_mapping_21atlases.csv").set_index("target")
chk("受体：α4β2 spin p = 0.060", round(float(rec.loc["a4b2-nACh", "spin_p"]), 3) == 0.060)
chk("受体：21 次检验 BH 最小校正 p = 0.4888（0/21 显著）",
    round(float(rec["fdr"].min()), 4) == 0.4888 and int((rec["fdr"] < 0.05).sum()) == 0,
    "min=%.4f n_sig=%d" % (float(rec["fdr"].min()), int((rec["fdr"] < 0.05).sum())))

par = load("receptor_partial_correlation.csv").set_index("target")
pr = float(par.loc["a4b2-nACh", "partial_r"])
chk("受体：偏相关 α4β2 r = +0.3014746（全精度）",
    abs(pr - 0.3014745623749493) < 5e-6, "got=%.10f" % pr)
pp = float(par.loc["a4b2-nACh", "partial_p_spin"])
chk("受体：偏相关 α4β2 spin p = 0.0744", abs(pp - 0.07442557442557443) < 5e-6, "got=%.10f" % pp)

# ---------------- 5 通路锚点 ----------------
n_alive = int((g["p_spin_one"] < 0.05).sum())
chk("通路：表型层存活 127/184", n_alive == 127, "alive=%d" % n_alive)

# ---------------- 6 细胞层网格 ----------------
p_grid = ct["p_prand_twotail"].values.astype(float)
resid = np.abs(p_grid * 10001 - np.round(p_grid * 10001))
chk("细胞：PRAND p 全部落在 (1+k)/10001 网格上", bool((resid < 1e-3).all()),
    "max|resid|=%.2e" % float(resid.max()))

# ---------------- 7 21×21 复算 ----------------
X = pv.iloc[:, 5:].values.astype(float)
C = np.corrcoef(X)
corr21 = load("receptor_21x21_correlation.csv")
corr21 = corr21.set_index(corr21.columns[0])
maxdiff = float(np.abs(C - corr21.values).max())
chk("21×21 相关矩阵可由 100 区全精度值复算（max|diff| < 1e-4）",
    maxdiff < 1e-4, "max|diff|=%.2e" % maxdiff)
tl = list(pv["target"])
r_pair = C[tl.index("a4b2-nACh"), tl.index("mu-opioid")]
chk("抽查：α4β2 与 μ-阿片 r = +0.515", round(r_pair, 3) == 0.515, "got=%.4f" % r_pair)

# ---------------- 汇总 ----------------
npass = sum(1 for _, ok, _ in results if ok)
print()
print("-" * 52)
print("合计 %d 条：PASS %d ｜ FAIL %d" % (len(results), npass, len(results) - npass))
if npass < len(results):
    print("重要：出现 FAIL——请检查数据完整性，并运行 verify_manuscript_claims.py 定位。")
else:
    print("全部通过。完整核对（96 条）请运行：python verify_manuscript_claims.py")
sys.exit(0 if npass == len(results) else 1)
