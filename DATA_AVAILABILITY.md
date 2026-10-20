# Data and code availability

This file contains the availability statements in the form required by NeuroImage, ready to be
pasted into the manuscript. The journal applies **Option C**: shared material must be deposited in
a repository, cited in the text with a link, and — where sharing is not possible — the reason must
be stated concretely. A bare "available upon request" is **not** accepted.

Two placeholders must be filled in before submission:

| Placeholder | What to put there |
|---|---|
| `<REPO_URL>` | Public repository URL (e.g. `https://github.com/<user>/PSF-disconnection-analysis`) |
| `<DOI>` | Zenodo (or equivalent) archive DOI, e.g. `10.5281/zenodo.XXXXXXX` |

---

## 1. English (paste into the manuscript)

### Data availability

> The raw imaging and clinical data contain identifiable patient information and therefore cannot
> be shared publicly, in accordance with the terms of the ethics approval for this study (Ethics
> Review Committee of the Second Affiliated Hospital, Army Medical University; approval no.
> 2024-研第 342-01). De-identified derived data — the regional white-matter disconnection matrix
> (124 subjects × 100 parcels, with study subject labels that map to no hospital identifier),
> PLS loadings, the gene loading ranking, and all cell-type, pathway and receptor results —
> together with the complete analysis code, are openly available at `<REPO_URL>` and archived at `https://doi.org/<DOI>`. Requests for
> access to the raw imaging and clinical data should be addressed to the corresponding author and
> will be reviewed by the ethics committee on a case-by-case basis.

### Code availability

> All analysis scripts, the de-identified derived data, and the third-party reference panels used
> in this study are available at `<REPO_URL>` and archived at `https://doi.org/<DOI>`. The
> repository contains the disconnection estimation (Lesion Quantification Toolkit), behavioral
> PLS, AHBA expression reconstruction and gene ranking, cell-type and pathway enrichment, PET
> receptor mapping, and every specificity analysis (spherical-spin permutation, leave-one-donor-out,
> control vectors) reported in the manuscript and the Supplementary Material.

### Third-party data and resources

> This study used the following openly available resources, whose own availability policies apply:
> the Allen Human Brain Atlas microarray data (Hawrylycz et al., 2012), the normative PET/SPECT
> receptor atlas collection distributed with `neuromaps` (Markello et al., 2022), the Schaefer-Yeo
> 100Parcels 7Networks cortical atlas (Schaefer et al., 2018), the normative cell-type gene sets
> of Seidlitz et al. (2020), and the HCP-842 normative connectome template distributed with the
> Lesion Quantification Toolkit (Griffis et al., 2021).

---

## 2. 中文（贴入中文稿）

### 数据与代码可得性

> 本研究的原始影像与临床数据含有可识别的患者信息，依据本研究的伦理批准条款（陆军军医大学
> 第二附属医院伦理审查委员会，批件号 2024-研第 342-01），不可公开共享。去标识的派生数据
> ——区域白质断连矩阵（124 例 × 100 区，受试者用研究编号，与任何院内标识无映射）、行为 PLS
> 载荷、基因载荷排序，以及全部细胞类型/通路/受体层结果——与完整分析代码一并公开存放于
> `<REPO_URL>`，并归档于 `https://doi.org/<DOI>`。原始影像与临床数据的获取请求请联系通讯
> 作者，由伦理委员会逐案审议。

### 代码可得性

> 本研究的全部分析脚本、去标识派生数据及所用的第三方参考面板均存放于 `<REPO_URL>`，
> 归档于 `https://doi.org/<DOI>`。包内含断连估计（Lesion Quantification Toolkit）、行为 PLS、
> AHBA 表达重建与基因排序、细胞类型与通路富集、PET 受体图谱映射，以及稿件与补充材料中报告的
> 全部特异性检验（球面旋转置换、留一供体、控制向量）脚本。

### 第三方数据与资源

> 本研究使用了以下公开资源，其各自的共享政策适用：艾伦人脑图谱表达数据（Hawrylycz et al.,
> 2012）、随 `neuromaps` 分发的规范化 PET/SPECT 受体图谱数据集（Markello et al., 2022）、
> Schaefer-Yeo 100Parcels 7Networks 皮层图谱（Schaefer et al., 2018）、Seidlitz 等（2020）的
> 规范细胞类型基因集，以及随 Lesion Quantification Toolkit 分发的 HCP-842 规范连接组模板
> （Griffis et al., 2021）。

---

## 3. Why the raw data cannot be shared (concrete reasons)

NeuroImage requires the reason to be specific rather than a generic statement. The grounds here are:

1. **The imaging data are identifiable.** Whole-brain MRI of a rare clinical population recruited
   at a single centre carries a realistic re-identification risk, and the lesion masks are derived
   directly from those scans.
2. **The clinical table contains direct and indirect identifiers** (age, sex, onset-to-scan
   interval, HAMD, and item-level FSS-9 responses) that could be linked to hospital
   records; the clinical table is therefore not distributed.
3. **The ethics approval obtained for this study does not cover public deposition.** It permits
   sharing only after review and approval by the ethics committee.
4. **What is shared instead.** The derived data released with the code let an independent group
   reproduce every statistical result, figure and table in the paper without any patient data;
   only the LQT step (Step 1) requires the lesion masks, which the ethics committee may release on
   request.

## 4. Suggested wording for the manuscript's Methods 2.10

The Chinese and English manuscripts currently state that the reproducible code package "will be
provided with the manuscript". Replace that sentence with a pointer to the repository so that the
Methods and the availability statement agree:

**EN** — *"The analysis scripts corresponding to the figures and tables, the version-controlled
code and the analysis pipeline are openly available at `<REPO_URL>` and archived at
`https://doi.org/<DOI>` (see Data and code availability)."*

**CN** — *"图与表对应的分析脚本、版本化代码及分析流程公开存放于 `<REPO_URL>`，归档于
`https://doi.org/<DOI>`（见数据与代码可得性）。"*
