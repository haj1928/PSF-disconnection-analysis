# Analysis code for: Fatigue-specific white-matter disconnection in acute ischemic stroke covaries with microglial, metabolic and nicotinic receptor maps

**Corresponding author:** [to be completed]
**Code version:** 2026-10-10
**Target journal:** NeuroImage

This repository contains the analysis code and de-identified derived data needed to reproduce
the results, figures and tables of the manuscript above. Raw imaging data and the clinical table
cannot be shared (see [`DATA_AVAILABILITY.md`](DATA_AVAILABILITY.md)).

> **中文读者／验证者**：请直接看 [`README_中文.md`](README_中文.md)。那是逐步操作指南，
> 包含"能验证到什么程度"的分档说明、环境准备、以及一条命令跑完的结论核对。

---

## 0. Quick verification (for reviewers)

```bash
python -m pip install pandas numpy scipy
python verify_manuscript_claims.py
```

`verify_manuscript_claims.py` re-derives **every numeric claim in the Results section**
(96 checks) from the CSVs in `derived_data/` alone — it does not read any author-local path and
does not access the network. Expected result: `PASS 96 | FAIL 0`.

One caveat worth knowing: `receptor_21x21_correlation.csv` stores values already rounded to four
decimals. Re-rounding those to three decimals turns −0.24849594 into −0.249 and produces a
**false** FAIL against the manuscript's −0.248. The checker therefore recomputes correlations
from the full-precision parcel values in `receptor_parcel_values_21x100.csv`.

**Additional tools bundled with the package**

- `repro_smoke_test.py` — ~5-second health check (16 assertions: file set, table shapes, key numeric anchors).
- `check_pkg_hygiene.py` — packaging hygiene (no `.bak`/`__pycache__`/lock files; key files and directories present).
- `CODE_PROVENANCE.md` — provenance & verification matrix: every script classified as verbatim (**V**) / ported (**P**) / self-implemented (**S**) / glue (**G**), with the evidence collected for each class.

---

## 1. What the study does

In 124 patients with acute ischemic stroke (within 14 days of onset) we estimated regional
white-matter disconnection probability from lesion masks with the Lesion Quantification Toolkit,
related it to post-stroke fatigue (FSS-9) using behavioral partial least squares (PLS), and then
asked whether the resulting brain pattern covaries with (i) AHBA cortical gene expression,
(ii) cell-type marker gene sets, (iii) enriched biological pathways, and (iv) normative
neurotransmitter receptor maps. All spatial statistics use spherical-spin permutation nulls.

## 2. Repository layout

```
scripts/
  step1_lqt_disconnection/    lesion masks, HCP-842 template -> 124 x 100 disconnection matrix (MATLAB + R)
  step2_behavioral_pls/       behavioral PLS on the disconnection matrix (Python / pyls)
  step3_ahba_expression/      abagen reconstruction of the AHBA expression matrix (100 x 15,633)
  step4_gene_ranking_celltype/gene loading ranking (R/SIMPLS) and cell-type enrichment
  step5_gsea/                 gseapy preranked GSEA (GO BP, KEGG)
  step6_pet_receptor/         21 normative PET/SPECT atlas correlations + spin nulls
  step7_specificity_controls/ specificity analyses (spin tests, LODO, control vectors)
  step8_figures/              figure generation (Fig. 1-5, Fig. S1); reads CSV only
resources/                    third-party reference panels (cell-type gene sets, HPA filtered list)
PSF_repro_PET/                spin-test code + fsaverage5 surface assets (see section 4)
derived_data/                 de-identified derived data (30 data files, see section 5)
```

## 3. Environment

Python 3.13 with the packages in [`requirements.txt`](requirements.txt); R 4.6.1 with the `pls`
package (SIMPLS) for Step 4; MATLAB with the Lesion Quantification Toolkit (LQT) for Step 1.

| Step | Language | Key dependency |
|---|---|---|
| 1 | MATLAB / R | LQT + HCP-842 normative connectome template; DSI Studio tractography |
| 2 | Python | `pyls` v0.0.1 |
| 3 | Python | `abagen` v0.1.3 |
| 4 | R | `pls` v2.9.0 (SIMPLS) |
| 5 | Python | `gseapy` v1.3.1 |
| 6, 7 | Python | `nibabel` v5.4.2, `neuromaps` atlases |

> **Note on paths.** Scripts use the authors' absolute paths (e.g. `F:/dcmtonii/...`) so that the
> analysis can be traced to the exact files that produced the published numbers. To run them
> elsewhere, change the `BASE`-style constants at the top of each script and place the inputs
> listed in section 4 accordingly. The spin-test module and fsaverage5 surface assets used by
> steps 4, 6 and 7 (and by the two surface figure scripts) ship in the companion folder
> [`PSF_repro_PET/`](PSF_repro_PET/); point `BASE` at the repository root, or `PET` directly at
> that folder, to resolve them.

## 4. Inputs the user must supply

| Input | Where to obtain |
|---|---|
| Lesion masks in MNI space | **Not shared** — subject data. Steps 2 onward can be run directly from `derived_data/disconnection_matrix_N124_cortex100.csv`. |
| HCP-842 normative connectome template | Ships with LQT (`lqt_parcellate`) |
| AHBA microarray data (6 donors) | `abagen` fetches it automatically (`abagen.fetch_microarray()`) |
| Schaefer-Yeo 100Parcels 7Networks atlas | `nilearn.datasets.fetch_atlas_schaefer_2018(n_rois=100, yeo_networks=7)` |
| 21 normative receptor/transporter atlases | `neuromaps` (Markello et al., 2022) |
| Gene sets (GO, KEGG GMT) | Enrichr / MSigDB; full lists in `resources/` |
| Cell-type gene sets | `resources/seidlitz2020_Dataset5_7class_ahba_aligned.{csv,json}` (Seidlitz et al., 2020, Supplementary Dataset 5) |
| Spin-rotation code & fsaverage5 surface assets | **Included** — `PSF_repro_PET/` (ENIGMA-toolbox release; see the README in that folder) |

## 5. Derived data shipped in this repository

All 30 data files are de-identified: they contain no names, dates or medical record numbers.
`disconnection_matrix_N124_cortex100.csv` carries a study subject label (`sub_10002` style) that has
no mapping to any hospital identifier and is not referenced in the manuscript. No clinical
variable is distributed.

| File | Content |
|---|---|
| `disconnection_matrix_N124_cortex100.csv` | 124 × 100 disconnection matrix (main analysis input `X`) |
| `volume_conservation_qc_N124.csv` | Per-patient native/MNI152 lesion-mask volumes and their ratio (`MNI/native`), i.e. the volume-conservation check of Methods 2.3 |
| `lesion_volume_descriptives.csv` | Spearman correlations of lesion volume with clinical variables (`rho`, `p`, `n`); the `FSS9总分` row backs the lesion-volume × FSS-9 check of Results 3.2 (ρ = +0.024, p = 0.791) |
| `pls_LV1_brain_loading_cortex100.csv` | LV1 brain loading + bootstrap ratio (100 parcels) |
| `pls_LV1_bootstrap_CI.csv` | LV1 correlation (the `obs_corr` used above); the file also carries its patient-resampling bootstrap CI, which is **not** a claim of the current manuscript (deleted in the 2026-09 revision) |
| `gene_ranking_w1_cortex100.csv` | Gene loading ranking `w1` (15,633 genes) |
| `celltype_enrichment_Dataset5.csv` | Cell-type enrichment (7 classes, PRAND) |
| `celltype_spatial_spin_Dataset5.csv` | Cell-type expression map correlations + spin *p* |
| `celltype_specificity_contrast.csv` | Cell-type expression-map correlations (LV1 loading; disconnection burden) with the partial correlation controlling the burden map, i.e. the specificity criterion |
| `gsea_GO_w1.csv`, `gsea_KEGG_w1.csv` | Full GSEA output (gseapy weighted KS) |
| `gsea_pathway_spin_survival.csv` | Pathway-level spin test (184 pathways) |
| `receptor_mapping_21atlases.csv` | Spatial correlation of 21 PET/SPECT atlases with LV1 |
| `receptor_parcel_values_21x100.csv` | **Full-precision** regional values for the 21 atlases (21 × 100) — use this to recompute the correlation matrix yourself |
| `receptor_21atlases_full_table.csv` | Full receptor table（Table S3 底表；含 20 种子汇总列） |
| `receptor_partial_correlation.csv` | Partial correlations controlling the other atlases — now carried by **Table 3** in the main text (Fig. 7 was deleted in the 2026-09-22 revision) |
| `receptor_cholinergic_family.csv`, `receptor_cholinergic_pairwise.csv` | Cholinergic family consistency |
| `figS1_receptor_corr.csv` | Atlas–atlas correlation matrix, wide layout — the file the Fig. S1 script reads |
| `receptor_21x21_correlation.csv` | Atlas–atlas correlation matrix, four-decimal rounding (used for cross-checks; see the caveat in section 0) |
| `multiseed_celltype_Dataset5_p.csv` | Cell-type 20-seed PRAND *p* |
| `multiseed_receptor_spin_p_matrix.csv` | 21 atlases × 20 seeds spin *p* |
| `lodo_donor_robustness_report.txt` | Leave-one-donor-out robustness (**Table 4 was deleted** in the 2026-09-22 revision; the per-donor *r* values are now quoted in the running text of §3.3) |
| `lodo_consistency_by_donor.csv` | Per-donor LODO vs all-donor loading correlation — quoted in the running text of §3.3 (formerly Fig. S1, left panel; that supplementary figure was deleted) |
| `lodo_common90_r_by_donor.csv` | Correlation recomputed on the 90 parcels common to all seven matrices (§3.3) |
| `sensitivity_specs.csv`, `core_pathway_results.csv` | Sensitivity analyses（见 §3.2；原图面面板已删） and core pathways |
| `pls_singval_spectrum.csv` | Variance explained by each latent variable (LV1 = 94.0%; reported in section 3.2, no longer a figure panel) |
| `pls_perm_null_singvals.csv` | Permutation null of the singular values, 9 LVs x 5,000 permutations (retained for verification; this panel was removed from Fig. 4) |
| `pls_perm_null_summary.csv` | Observed singular value, null statistics and permutation *p* per LV (retained for verification; this panel was removed from Fig. 4) |

> **Archived, not deleted.** `lodo_top300_jaccard_by_donor.csv` (Jaccard overlap of top-300 gene
> lists) was dropped from the manuscript together with supplementary Fig. S1 on 2026-09-22. It was
> kept in this repository under `derived_data/legacy_removed_20260923/` until 2026-10-10, when it
> was **archived out of the repository** (moved, not deleted) together with the other superseded
> `legacy_removed_20260923/` items. It is no longer a claim of the manuscript and is not checked
> by `verify_manuscript_claims.py`.
>
> **Script file names match the figure numbers** (renamed 2026-09-24):
> `fig1_*` = Fig. 1, `fig2_*` = Fig. 2, `fig3_*` = Fig. 3, `fig4_celltype_*` = Fig. 4,
> `fig5_pathway_enrichment_*` = Fig. 5, `figS1_*` = Fig. S1 (the receptor-correlation matrix figure; earlier internal names used figS2). Scripts belonging to figures that
> were deleted during revision were moved out of the repository on 2026-10-08 (archived, not deleted).

## 6. Analysis order

| Step | Script(s) | Key output |
|---|---|---|
| 1 | `build_lqt_input.py` → LQT (MATLAB/R in `lqt_matlab_*.m`, `psf_sdc_lqt_n124_v2.R`) → `aggregate_lqt_to_matrix.py` | disconnection matrix |
| 2 | `pls_pyls_discon_cov_residY_cortex100.py`, `run_pls_bootCI_cortex100.py` | LV1 loading / BSR |
| 3 | `gen_ahba_100.py` | 100 × 15,633 expression matrix |
| 4 | `gene_ranking_step4b.py` → `run_celltype_enrich_seidlitz_Dataset5.py` | `w1` ranking; cell-type enrichment |
| 5 | `run_gsea_gseapy_w1_cortex100_residY.py` | GO / KEGG GSEA tables |
| 6 | `run_receptor_discon_cortex100_residY.py`, `run_spin_partial_receptor.py` | receptor correlations + partials |
| 7 | `build_discon_control_vector.py`, `build_salience_permutations.py`, `run_sal_specificity_spin_Dataset5.py`, `run_gsea_pathway_spin_official.py`, `run_lodo_leave_one_donor.py` | specificity results |
| 8 | `fig1_lesion_distribution_20260916.py` (Fig. 1), `fig2_behavioral_pls_20260916.py` (Fig. 2), `fig3_robustness_20260916.py` (Fig. 3), `fig4_celltype_20260916.py` (Fig. 4), `fig5_pathway_enrichment_20260916.py` (Fig. 5), `figS_all_20260916.py` (Fig. S1), `figstyle_20260916.py` | Fig. 1–5, Fig. S1 |

`scripts/step8_figures/figstyle_20260916.py` is the shared style module (Elsevier/NeuroImage
geometry: 90 / 140 / 190 mm, Arial, 7–9 pt labels) and **must be importable** when running any
figure script. Figure scripts read CSV only — they contain no hard-coded numbers.

Scripts for figures deleted during revision were moved out of the repository on 2026-10-08
(archived, not deleted); the scripts shipped here all belong to the current manuscript —
Fig. 1–5 and Fig. S1. For the record, the removed set was: `fig1_cohort_*` and the earlier
`fig1_make_editable_pptx_*` (the deleted Fig. 1 — cohort characteristics and participant
flow), `fig8_receptor_atlas_*` (the deleted Fig. 7 — its numbers are now carried by Table 3)
and `figS1_gene_lodo_*` (the deleted supplementary Fig. S1 — donor robustness of the gene
ranking), together with their `data_*` helpers (`data_fig5_gene_lodo_*`,
`data_fig8_receptor_*`). The `fig1_make_editable_pptx_20260916.py` that ships here is the
**current** Fig. 1 export script (a different file that happens to share the name).

## 7. Expected outputs (verification table)

These values are computed **from the CSV files shipped in `derived_data/`**, so a successful
re-run should reproduce them. All *p* values quoted in the manuscript are rounded half-up to
three decimals from these full-precision values.

**Step 1 — disconnection matrix:** 124 subjects × 100 parcels, values 0.0000 – 1.0000
(regional disconnection probability).

**Step 2 — behavioral PLS:** LV1 loading positive in 80/100 parcels;
`obs_corr = 0.3171`.

**Step 4 — cell-type enrichment** (median rank against a PRAND null, BH-FDR over 7 classes):

| Cell type | n genes | median rank | BH-FDR | spatial *r* | spin *p* |
|---|---|---|---|---|---|
| Excitatory neurons | 1,132 | 9,192 | 1.75e-4 | −0.026 | 0.088 |
| **Microglia** | 615 | **6,242** | 1.75e-4 | **+0.290** | **0.044** |
| **Oligodendrocytes** | 783 | **5,232** | 1.75e-4 | **+0.212** | **0.020** |
| Inhibitory neurons | 829 | 9,749 | 1.75e-4 | −0.120 | 0.034 |
| Astrocytes | 777 | 7,157 | 0.027 | +0.122 | 0.309 |
| Endothelial | 839 | 7,237 | 0.038 | +0.206 | 0.070 |
| OPC | 141 | 7,823 | 0.993 | +0.122 | 0.337 |

**Step 5 — GSEA:** 134 / 5,301 GO BP terms and 50 / 312 KEGG terms at FDR *q* < 0.05;
pathway-level spin test performed on 184 pathways.

**Step 6 — PET receptor mapping** (21 atlases; spin *p* from 1,000 rotations) — reported in
**Table 3** of the main text; the corresponding figure (Fig. 7) was deleted on 2026-09-22:

| Atlas | Target | *r* | spin *p* | BH-FDR | partial *r* |
|---|---|---|---|---|---|
| Hillmer 2016 (flubatine) | α4β2-nACh | +0.383 | 0.060 | 0.489 | **+0.301** |
| Kantonen 2020 (carfentanil) | μ-opioid | +0.274 | 0.168 | 0.489 | +0.009 |
| Nørgaard 2021 (flumazenil) | GABA-A | −0.255 | 0.140 | 0.489 | −0.151 |

No atlas survives BH correction across the 21 primary tests (minimum corrected *p* = 0.489);
the three rows above are the most strongly associated atlases — not FDR-significant ones.

## 8. Known limitations

- The AHBA and receptor atlases come from healthy donors, so Step 4–6 results describe spatial
  covariation with normative maps, not measurements in the patients studied here.
- AHBA contains six donors, only two of them with full cortical coverage; leave-one-donor-out
  robustness is reported in Step 7.
- Cell-type enrichment uses the 7-class scheme of Seidlitz et al. (2020); finer cell types are
  not resolvable with this panel.
- The glial spatial associations should be read as one shared signal (microglia and
  oligodendrocytes largely share their spatial pattern; see §3.4). The microglial association is
  the one closest to the significance boundary (spin *p* = 0.044).

## 9. Citation

If you use this code or the derived data, please cite the manuscript above and the third-party
resources listed in `DATA_AVAILABILITY.md`.

## 10. License

The code in this repository is released under the MIT License (see `LICENSE`).
The derived data in `derived_data/` are released for research use only; no
patient-level data are included. Third-party components — the spin-test code and
fsaverage5 surface assets in `PSF_repro_PET/`, and the reference atlases — retain
their own licenses; see `PSF_repro_PET/README.md`.
