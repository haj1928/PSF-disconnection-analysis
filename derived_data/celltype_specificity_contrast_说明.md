# `celltype_specificity_contrast.csv` 列用途说明（2026-09-23 核校；2026-10-10 数值同步）

本说明回答一个问题：**这个文件里哪几列仍被正文使用、哪几列已弃用**。
结论先行：

- **现文件（`celltype_specificity_contrast.csv`）的 8 列全部仍在用** —— 它是细胞类型层
  「特异性判据」的唯一载体，正文 §3.4 的相关结论由它支撑。
- **已弃用的是 Δr 配对检验的两列 `diff` / `p_diff`** —— 它们**已经不在现文件里**，
  只保留在同目录的备份 `celltype_specificity_contrast.csv.bak_before_nodeltaR_20260922`。
- 因此本文件**不要删除**（用户要求），备份也不要删除（它是 Δr 口径的唯一留档）。

## 1. 现文件的 8 列（仍在使用）

| 列 | 含义 | 用途 |
|---|---|---|
| `cell_type` | 细胞类型（7 类，Seidlitz 2020 Dataset5 方案） | 行标签 |
| `n_genes` | 该类基因集的基因数 | 与 `celltype_enrichment_Dataset5.csv` 一致 |
| `r_salience` | 该类表达图与 **LV1 载荷图**的空间皮尔逊 *r* | §3.4 的"是否与疲劳相关模式共变" |
| `p_salience` | 上者的 **spin 置换 p**（1,000 次球面旋转） | 显著性判据 |
| `r_burden` | 该类表达图与**断连负荷图**的空间 *r* | 对照：排除"只是跟着病灶负荷走" |
| `p_burden` | 上者的 spin 置换 p | 对照的显著性 |
| `r_partial` | **控制断连负荷图后**的偏相关 *r* | ★ **特异性判据**：控掉负荷后仍与 LV1 共变才算特异 |
| `p_partial` | 上者的 spin 置换 p | 特异性判据的显著性 |

三层递进关系即特异性论证：`r_salience`（共变）→ `r_burden`（是否与负荷混淆）
→ `r_partial`（控负荷后是否还在）。**三层的 6 个数值列都被正文使用**，缺一不可。

现文件取值（3 位小数；与 `celltype_spatial_spin_Dataset5.csv` 的全精度值一致，两表均在核对脚本读取范围内）：

| cell_type | n_genes | r_salience | p_salience | r_burden | p_burden | r_partial | p_partial |
|---|---|---|---|---|---|---|---|
| Micro | 615 | +0.290 | 0.044 | −0.212 | 0.176 | +0.258 | 0.051 |
| Oligodendrocytes | 783 | +0.212 | 0.020 | −0.046 | 0.389 | +0.207 | 0.019 |
| Astro | 777 | +0.122 | 0.309 | −0.057 | 0.387 | +0.113 | 0.311 |
| Endo | 839 | +0.206 | 0.070 | −0.202 | 0.178 | +0.172 | 0.084 |
| Excitatory_Neurons | 1132 | −0.026 | 0.088 | +0.047 | 0.623 | −0.017 | 0.103 |
| Inhibitory_Neurons | 829 | −0.120 | 0.034 | −0.005 | 0.282 | −0.124 | 0.028 |
| OPC | 141 | +0.122 | 0.337 | −0.160 | 0.194 | +0.092 | 0.354 |

## 2. 已弃用的 2 列（Δr 配对检验）

| 列 | 含义 | 状态 |
|---|---|---|
| `diff` | Δr = `r_partial` − `r_salience`（控制负荷图后 *r* 的变化量） | **已弃用**（下表值为坐标修正前的历史留档） |
| `p_diff` | 该 Δr 的配对／置换检验 p | **已弃用** |

- 这两列在 **2026-09-22** 的修订中随"Δr 配对检验"一并弃用：
  判定特异性只需看 `r_partial` 本身是否显著（控负荷后是否还在），
  再对"`r` 掉了多少"做一次配对检验属于**二次推断**，且 Δr 的零分布难以用
  spin 置换恰当地构造，故稿件不再报告。
- 备份文件名 `celltype_specificity_contrast.csv.bak_before_nodeltaR_20260922`
  中的 `before_nodeltaR` 即指"删除 Δr 列之前"的版本；把现文件与备份对比
  可以看到差异**恰好**是这两列被去掉，其余 8 列的数值逐格相同（未做任何数值改动）。
- 备份中 7 行的 `diff` / `p_diff` 取值留档如下：

| cell_type | diff | p_diff |
|---|---|---|
| Micro | +0.5021 | 0.0020 |
| Oligodendrocytes | +0.2585 | 0.1159 |
| Astro | +0.1792 | 0.2697 |
| Endo | +0.4080 | 0.0130 |
| Excitatory_Neurons | −0.0727 | 0.6264 |
| Inhibitory_Neurons | −0.1148 | 0.4585 |
| OPC | +0.2820 | 0.0819 |

## 3. 与核对脚本的关系

- `verify_manuscript_claims.py` 的 §3.4 **同时读取 `celltype_spatial_spin_Dataset5.csv`**
  （`salience` / `discon_mean` 两种 mode）**与本文件的偏相关列**；两者数值一致，互为交叉验证。
- 本文件的 `r_partial` / `p_partial`（控负荷后的偏相关）与抑制性神经元行
  已于 **2026-10-10 纳入 `verify_manuscript_claims.py` §3.4 核对**（原缺口已闭合）。
- 本文件自 2026-10-10 起**已被 `verify_manuscript_claims.py` 消费**，同时承载正文 §3.4 的
  特异性判据；按用户要求**保留文件与备份**（未删除、未移出本仓库）。
