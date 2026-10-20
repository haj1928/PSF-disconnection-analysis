# Step 4 — Gene loading ranking and cell-type enrichment

**Purpose.** Rank genes by their spatial loading on the fatigue-related disconnection pattern
(`w1`), then test whether canonical cell-type marker gene sets sit at the front or back of that
ranking.

**Inputs.** Step 3 expression matrix; Step 2 LV1 loading.

**Scripts.**
| File | Role |
|---|---|
| `gene_ranking_step4b.py` | per-gene PLS loading ranking `w1` (R, `pls` package, SIMPLS) |
| `run_celltype_enrich_seidlitz_Dataset5.py` | cell-type enrichment, primary (PRAND null) |
| `run_celltype_enrich_Dataset5_hpafilter.py` | same, after HPA brain-expression filtering |
| `check_multiseed_Dataset5.py` | 20-seed robustness of the PRAND *p* values |
| `check_layer_independence_gene_overlap.py` | gene-set overlap between cell and pathway layers + hypergeometric null |
| `celltype_specificity_contrast.py` | specificity contrast between cell types |

**Outputs.** `derived_data/gene_ranking_w1_cortex100.csv`,
`derived_data/celltype_enrichment_Dataset5.csv`,
`derived_data/celltype_spatial_spin_Dataset5.csv`,
`derived_data/celltype_specificity_contrast.csv`,
`derived_data/multiseed_celltype_Dataset5_p.csv`.

**Notes.** Each gene is z-scored across parcels before ranking. Cell-type gene sets are the
7-class scheme of Seidlitz et al. (2020), shipped in `resources/`. Expected values: main README
section 7.
