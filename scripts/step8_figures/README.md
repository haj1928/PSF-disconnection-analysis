# Step 8 — Figures

**Purpose.** Reproduce Fig. 1–5 and Fig. S1 of the manuscript.

**Scripts.** `fig1_lesion_distribution_*`, `fig2_behavioral_pls_*`, `fig3_robustness_*`,
`fig4_celltype_*` and `fig5_pathway_enrichment_*` draw Fig. 1–5;
`figS_all_20260916.py` draws Fig. S1 (receptor–receptor correlation matrix); the
`EMIT_UNUSED` switch is `False`, so material not used in the manuscript is not regenerated;
`data_fig*_*.py` extract the numbers from `derived_data/` into small CSVs that the plotting
scripts read.

**`figstyle_20260916.py` is required.** It is the shared style module — figure scripts import it
directly, so it must stay in the same directory.

**Notes.** Plotting scripts **read CSV only and contain no hard-coded numbers**, so a change in
the analysis propagates automatically. Geometry follows the Elsevier/NeuroImage specification
(90 / 140 / 190 mm column widths, Arial, 7–9 pt labels, 0.6 pt lines). Numeric labels are rounded
half-up through `figstyle.f3()` for consistency with the manuscript.
