# Step 5 — Pathway enrichment (GSEA)

**Purpose.** Test which GO Biological Process and KEGG pathways are enriched along the `w1` gene
ranking, using the official weighted Kolmogorov–Smirnov statistic.

**Inputs.** `derived_data/gene_ranking_w1_cortex100.csv`; GO / KEGG GMT gene sets.

**Scripts.**
| File | Role |
|---|---|
| `run_gsea_gseapy_w1_cortex100_residY.py` | primary GSEA (`gseapy` preranked, weighted KS) |
| `run_gsea_hpafilter.py` | same, restricted to HPA cortex-detected genes |
| `run_gsea_discon_burden.py` | burden control: ranking derived from total disconnection |
| `run_gsea_control_vectors.py` | further control vectors |

**Outputs.** `derived_data/gsea_GO_w1.csv`, `derived_data/gsea_KEGG_w1.csv`.

**Notes.** The official `gseapy` weighted-KS result is the single primary GSEA readout; a
self-implemented mean-ES (Tian 2005) variant is retained only as a sensitivity analysis.
Expected values: main README section 7.
