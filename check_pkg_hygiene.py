# -*- coding: utf-8 -*-
"""check_pkg_hygiene.py —— 提交包卫生检查（打包/提交前运行）

检查：
  1  无 .bak / __pycache__ / *.tmp* 残留（递归，.git 除外）
  2  无 Office 锁文件（~$*）
  3  关键文件齐全（README×2 / CODE_PROVENANCE / DATA_AVAILABILITY / verify / smoke / 本脚本）
  4  关键目录齐全（scripts/step1–8 / PSF_repro_PET / resources / derived_data）
  5  各步骤脚本数快照 + derived_data 文件数（信息项）

用法：python check_pkg_hygiene.py
"""
import os
import sys

try:
    sys.stdout.reconfigure(encoding="utf-8")
except Exception:
    pass

HERE = os.path.dirname(os.path.abspath(__file__))
results = []


def chk(label, ok, detail=""):
    results.append((label, bool(ok)))
    print("[%s] %s%s" % ("PASS" if ok else "FAIL", label, ("  ← " + detail) if detail else ""))


print("=== check_pkg_hygiene.py — 提交包卫生检查 ===")
print()

# ---------------- 1 残留扫描 ----------------
junk = []
for dp, dn, fn in os.walk(HERE):
    parts = dp.split(os.sep)
    if ".git" in parts:
        continue
    for d in dn:
        if d == "__pycache__" or d.startswith(".bak"):
            junk.append(os.path.relpath(os.path.join(dp, d), HERE))
    for f in fn:
        if ".bak" in f or f.endswith((".tmp", ".tmp_pathway", ".tmp_receptor")):
            junk.append(os.path.relpath(os.path.join(dp, f), HERE))
chk("无 .bak / __pycache__ / *.tmp* 残留", not junk,
    ("发现 %d 项: " % len(junk)) + "; ".join(junk[:6]) + (" …" if len(junk) > 6 else "") if junk else "")

# ---------------- 2 锁文件 ----------------
locks = []
for dp, dn, fn in os.walk(HERE):
    if ".git" in dp.split(os.sep):
        continue
    for f in fn:
        if f.startswith("~$"):
            locks.append(os.path.relpath(os.path.join(dp, f), HERE))
chk("无 Office 锁文件（~$*）", not locks, "; ".join(locks[:5]))

# ---------------- 3 关键文件 ----------------
must_files = [
    "README.md", "README_中文.md", "CODE_PROVENANCE.md", "DATA_AVAILABILITY.md",
    "verify_manuscript_claims.py", "repro_smoke_test.py", "check_pkg_hygiene.py",
]
missing = [f for f in must_files if not os.path.exists(os.path.join(HERE, f))]
chk("关键文件齐全（%d/%d）" % (len(must_files) - len(missing), len(must_files)),
    not missing, ("缺: " + ", ".join(missing)) if missing else "")

# ---------------- 4 关键目录 ----------------
must_dirs = [
    "scripts/step1_lqt_disconnection", "scripts/step2_behavioral_pls",
    "scripts/step3_ahba_expression", "scripts/step4_gene_ranking_celltype",
    "scripts/step5_gsea", "scripts/step6_pet_receptor",
    "scripts/step7_specificity_controls", "scripts/step8_figures",
    "PSF_repro_PET", "resources", "derived_data",
]
missing_d = [d for d in must_dirs if not os.path.isdir(os.path.join(HERE, d))]
chk("关键目录齐全（%d/%d）" % (len(must_dirs) - len(missing_d), len(must_dirs)),
    not missing_d, ("缺: " + ", ".join(missing_d)) if missing_d else "")

# ---------------- 5 计数快照（信息项） ----------------
snap = {}
for d in ("scripts/step1_lqt_disconnection", "scripts/step2_behavioral_pls",
          "scripts/step3_ahba_expression", "scripts/step4_gene_ranking_celltype",
          "scripts/step5_gsea", "scripts/step6_pet_receptor",
          "scripts/step7_specificity_controls", "scripts/step8_figures",
          "PSF_repro_PET", "resources"):
    base = os.path.join(HERE, d)
    n = 0
    for dp, dn, fn in os.walk(base):
        if "__pycache__" in dp or ".bak" in dp:
            continue
        n += sum(1 for f in fn if f.endswith(".py"))
    snap[d] = n
print("\n[info] 脚本计数快照: " + " | ".join("%s=%d" % (k.split("/")[-1], v) for k, v in snap.items()))
n_dd = len([f for f in os.listdir(os.path.join(HERE, "derived_data"))
            if os.path.isfile(os.path.join(HERE, "derived_data", f))])
print("[info] derived_data 文件数: %d" % n_dd)

# ---------------- 汇总 ----------------
npass = sum(1 for _, ok in results if ok)
print()
print("-" * 52)
print("合计 %d 条：PASS %d ｜ FAIL %d" % (len(results), npass, len(results) - npass))
if npass == len(results):
    print("卫生检查通过。完整校验：python verify_manuscript_claims.py（96 条）；")
    print("快速健康检查：python repro_smoke_test.py。")
else:
    print("FAIL：请先清理残留（.bak/__pycache__ 移至包外，勿直接删除），再复跑。")
sys.exit(0 if npass == len(results) else 1)
