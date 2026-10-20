# CODE_PROVENANCE.md — Code provenance and verification matrix

**Code version:** 2026-10-10 · **Package:** `PSF_code_and_data_submission_20260917`

This document classifies **every script** in this package by provenance and lists the
verification evidence collected for it. It answers two questions a reviewer may have:

1. **Which parts are official/third-party code, and which are ours?**
2. **How has each part been checked?**

## Categories

| Class | Meaning | How it is verified |
| --- | --- | --- |
| **V — Verbatim** | Copied verbatim from an official source; **must not be edited locally** (any local modification demotes the copy to *P*). | Line/character comparison against the official source; identical outputs in straight-run comparison. |
| **P — Ported** | An official algorithm re-implemented (e.g. numpy port, vectorized); semantically equivalent. | **Bitwise-equal output** vs the official implementation under identical inputs and seeds. |
| **S — Self-implemented** | No official counterpart exists in our pipeline. | Independent double-implementation, closed-form equivalence, property-based checks (e.g. permutation p-value grids), or in-package assertions. |
| **G — Glue** | Data I/O, reshaping, file management, plotting. | Schema/consistency assertions (shape, NaN alignment, column checks) at runtime. |

## Shared external components

| Component | Source | Used by | Verification evidence |
| --- | --- | --- | --- |
| `centroid_extraction_sphere`, `rotate_parcellation`, `perm_sphere_p` | **enigmatoolbox** (MICA-MNI), *verbatim* in `PSF_repro_PET/spin_test_numpy.py` | all spin/rotation analyses (steps 4–8) | semantic check vs official toolbox; official straight-run vs our pipeline: `max|diff| = 0.000e+00`; single shared source (no per-script re-implementation) |
| `pyls` | Python PLS package | `step2` | official package; outputs cross-checked against archived values (corr = +0.317) by `verify_manuscript_claims.py` |
| `gseapy` (prerank, weighted KS) | gseapy package | `step5` | official implementation defines the primary result; independent recomputation of ES by our vectorized weighted-KS: `max|diff| ≈ 1e-13` (column `es_diff` in `gsea_pathway_spin_survival.csv`) |
| `abagen` | abagen (Markello et al., 2021) | `step3` | default parameters; built-in self-check "first 20 common regions `max|diff| ≈ 0`"; donor archives CRC-verified |
| `Seidlitz et al. (2020)` 7-class gene sets (Dataset 5) | bundled JSON in `resources/` | `step4` | gene-set file bundled with the package; median-rank + PRAND statistics recomputed in-package |
| Lesion Quantification Toolkit (LQT) + HCP-842 tractography template | external (Griffis et al., 2021) | `step1` | LQT outputs aggregated by our scripts; cross-checks inside `crop_verify_cortex100.py` |

> **Seam checks (the interface between official and our code)** are the historically
> error-prone zone: coordinate↔label pairing, parameter order, 0/1-based indexing, LH/RH
> ordering. All `rotate_parcellation` call sites were re-audited on 2026-10-10 and use
> row-index alignment (`enumerate` + `labels == i`), not name–centroid zipping.

---

## Per-script matrix

### PSF_repro_PET/ (shared library)

| Script | Class | Notes | Verification evidence |
| --- | --- | --- | --- |
| `spin_test_numpy.py` | **V** (+G demo) | The three rotation functions are verbatim from enigmatoolbox (MICA-MNI, @saratheriver); header states provenance | official straight-run `max|diff| = 0.000e+00`; single source reused by every spin script |

### resources/supplementary_code/ (sensitivity add-ons)

| Script | Class | Notes | Verification evidence |
| --- | --- | --- | --- |
| `check_fss_floor_addendum.py`（已移出） | — | FSS floor-count addendum，含已于 2026-09-22 从稿件删除的下限效应/二分类口径；**2026-10-10 移出至审计归档**（移动、未删除） | n/a（不在包内） |
| `sens_lesion_volume.py`, `sens_lesion_volume_log.py` | S | lesion-volume sensitivity ladder（`_log` 的杠杆段已于 2026-09-22 删除；2026-10-10 清理残留引用） | NaN guards at runtime; outputs archived under `derived_data/` |
| `sens_gsea_plsw1_input.py` | S | GSEA input cross-check | NES/FDR cross-table vs re-run pipelines |
| `legacy_removed_20260923/`（2 处目录，已移出） | — | **Historical, superseded analyses** removed from the manuscript；**2026-10-10 整体移出至审计归档**（含 `data_fig3_loo_logvol_20260916.py` 同批；移动、未删除） | n/a（不在包内）；not used by any current result |

### scripts/step1_lqt_disconnection/

| Script | Class | Notes | Verification evidence |
| --- | --- | --- | --- |
| `build_lqt_input.py` | G | assembles LQT inputs | runtime shape checks |
| `aggregate_lqt_to_matrix.py` | S | aggregates LQT outputs → 124×100 disconnection matrix | runtime assertions: shape, NaN-pattern equality between matrices |
| `crop_verify_cortex100.py` | G | cortex-100 cropping verification | built-in row/column comparison log |

### scripts/step2_behavioral_pls/

| Script | Class | Notes | Verification evidence |
| --- | --- | --- | --- |
| `pls_pyls_discon_cov_residY_cortex100.py` | G+S | **calls official `pyls`**; our part = age/sex residualization and permutation scheme | `corr = +0.317`, `p_perm = 0.023` reproduced by `verify_manuscript_claims.py` |
| `run_pls_bootCI_cortex100.py` | S | bootstrap CI (1,000×, seed = 0) | `pls_LV1_bootstrap_CI.csv` (obs_corr = 0.3171) retained in package |

### scripts/step3_ahba_expression/

| Script | Class | Notes | Verification evidence |
| --- | --- | --- | --- |
| `gen_ahba_100.py` | G (official `abagen`) | default abagen parameters after Markello et al. (2021) | built-in self-check: first 20 common regions `max|diff| ≈ 0`; donor CRC verification |

### scripts/step4_gene_ranking_celltype/

| Script | Class | Notes | Verification evidence |
| --- | --- | --- | --- |
| `gene_ranking_step4b.py` | S | w1 gene ranking (rank order, ties) | ranked counts and loadings cross-checked by `verify_manuscript_claims.py` |
| `run_celltype_enrich_seidlitz_Dataset5.py` | S | median-rank + PRAND (inclusive two-tailed) + BH within family | p-values match the (1+k)/(N+1) grid; manuscript values reproduced |
| `run_celltype_enrich_Dataset5_hpafilter.py` | S | same as above with HPA cortical filter (15,633 → 12,753 genes) | gene counts and BH-adjusted values reproduced |
| `check_multiseed_Dataset5.py` | S | 20-seed stability of PRAND | per-seed values archived (`multiseed_celltype_Dataset5_p.csv`) |
| `celltype_specificity_contrast.py` | S | closed-form partial correlation (Table S1 column 3 ledger); **2026-10-09 fix: row-index alignment** (former zip pairing was defective); NaN fix (0-fill for the missing parcel) documented inline | closed-form ≡ residual method (err 1.4e-16); call sites re-audited 2026-10-10 |
| `check_layer_independence_gene_overlap.py` | S | gene-overlap layer-independence check (addendum) | Jaccard values archived |

### scripts/step5_gsea/

| Script | Class | Notes | Verification evidence |
| --- | --- | --- | --- |
| `run_gsea_gseapy_w1_cortex100_residY.py` | G (official `gseapy`) | prerank weighted KS — **defines the primary GSEA result** | official implementation by design |
| `run_gsea_control_vectors.py` | S | vectorized weighted KS + label permutation on control input | ES recompute vs gseapy `max|diff| ≈ 1e-13`; p-values on inclusive grid |
| `run_gsea_discon_burden.py` | S | same pipeline on disconnection-burden ranking | same cross-checks |
| `run_gsea_hpafilter.py` | S | same pipeline under HPA filter | same cross-checks |

### scripts/step6_pet_receptor/

| Script | Class | Notes | Verification evidence |
| --- | --- | --- | --- |
| `run_receptor_discon_cortex100_residY.py` | S+P | parcel-value extraction (nearest-neighbour resampling + within-parcel `nanmean`); P = ported spin rotation; **inclusive (1+k)/(N+1) counting, two-direction average** | sandbox rerun 2026-10-10: 21/21 spin p bitwise-identical to manuscript values (0.0604/0.1678/…; BH min 0.4888) |
| `run_spin_partial_receptor.py` | S+P | partial-correlation mirror (other atlases held fixed; both rotation directions averaged) | sandbox rerun 2026-10-10: T3 partial p = 0.0744/0.4685/0.2677 reproduced |
| `run_cholinergic_family.py` | S+P | within-cholinergic-family contrast | sandbox rerun: family values 0.0410/0.2098/0.3017 reproduced |
| `run_multiseed_spin_21atlas.py` | P | vectorized `perm_sphere_p` (equivalence stated inline); inclusive two-direction counting | 20-seed matrix seed0 row bitwise-identical to main analysis (21/21, 2026-10-10) |

### scripts/step7_specificity_controls/

| Script | Class | Notes | Verification evidence |
| --- | --- | --- | --- |
| `build_discon_control_vector.py` | G | builds the disconnection-burden control vector | runtime checks |
| `build_salience_permutations.py` | S+G | builds the 1,000-rotation salience archive | shares `rotate_parcellation` (single source) |
| `crosscheck_gsea_discon_burden.py` | S | cross-check of GSEA under the burden control | archived cross-table |
| `run_gsea_pathway_spin_official.py` | S | vectorized weighted KS + phenotype-level spin; inclusive counting, denominator `K = N+1 = 1001` (verified 2026-10-10) | pathway survival 127/57 reproduced; matches manuscript §3.5(ii) |
| `run_lodo_leave_one_donor.py` | S | leave-one-donor-out analysis | outputs archived (`lodo_*` CSVs) |
| `run_lodo_sensitivity_common94.py` | S | LODO on the common parcels（终稿 90 区；见 `lodo_common90_r_by_donor.csv`） | outputs archived |
| `run_sal_specificity_spin_Dataset5.py` | S+P | salience-specificity spin (cell layer) | p on inclusive grid; values reproduced |
| `run_sal_specificity_spin_proteostasis.py` | S+P | salience-specificity spin (proteostasis) | values reproduced |
| `run_sal_specificity_spin_v2.py` | S+P | salience-specificity spin (main); includes an inline note on the corrected pathway-name matching (former `startswith` fuzzy match removed) | values reproduced |

### scripts/step8_figures/ (figure generation; all Class G unless noted)

| Script | Class | Notes | Verification evidence |
| --- | --- | --- | --- |
| `figstyle_20260916.py` | G | shared style/palette library | palette constants asserted once, imported everywhere |
| `data_fig3_*` (3 scripts) | G+S | figure-3 data（`data_fig3_loo_logvol_20260916.py` 已于 2026-10-10 随杠杆口径移出至审计目录）; `data_fig3_perm_null_20260918.py` recomputes permutation p inclusive `(sum+1)/(N+1)` | inclusive counting visible in source; values match figure |
| `data_fig5_pathway_20260916.py` | G+S | figure-5 data; recomputes display r with `Decimal` ROUND_HALF_UP (single rounding) | inline comment documents the former double-rounding pitfall |
| `data_figS_20260916.py`, `data_unused_receptor_parcels_20260916.py` | G | supplementary-figure data; unused-atlas archival | runtime checks |
| `fig1_build_prob_maps_20260916.py`, `fig1_lesion_distribution_20260916.py`, `fig1_make_editable_pptx_20260916.py` | G | figure 1 (probability maps, lesion distribution, editable pptx) | shape/range checks |
| `fig2_behavioral_pls_20260916.py` | G | figure 2 (PLS results; BSR filtering by `|BSR|>2`) | filter counts reproduced |
| `fig3_robustness_20260916.py`, `fig4_celltype_20260916.py`, `fig5_pathway_enrichment_20260916.py`, `figS_all_20260916.py` | G | figures 3–5, S1 | pixel-stability checks via PNG byte comparison |

---

## Reproducibility conventions

- **Randomness:** every stochastic analysis records and fixes its seed (all `numpy`
  `default_rng(seed)`); seeds are printed in the logs and stored in output CSVs where relevant.
- **Rounding:** intermediate values are stored at full precision; rounding to report precision
  (three decimals, half-up) is applied **once**, at the last step.
- **Permutation p-values:** inclusive counting `(1+k)/(N+1)` (the observed value is part of the
  null distribution) — identical across layers.
- **Spin tests:** spherical-spin null (Alexander-Bloch et al., 2018), 1,000 rotations,
  seed = 0, 100-parcel Schaefer-Yeo space, two-direction average.
- **Standing checks:** `repro_smoke_test.py` (fast health check) and
  `verify_manuscript_claims.py` (all manuscript claims, from `derived_data/` only).
