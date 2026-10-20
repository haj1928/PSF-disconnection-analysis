# Step 2 — Behavioral PLS: disconnection–fatigue covariation

**Purpose.** Extract the latent variable linking the disconnection matrix to FSS-9, with the
behavioral side residualised on age and sex only (X left untouched).

**Inputs.** `derived_data/disconnection_matrix_N124_cortex100.csv`; FSS-9 item scores and
age/sex (not shared — patient data).

**Scripts.**
| File | Role |
|---|---|
| `pls_pyls_discon_cov_residY_cortex100.py` | main behavioral PLS (`pyls`), LV1 loading + BSR (Y residualised on age and sex; X left untouched) |
| `run_pls_bootCI_cortex100.py` | bootstrap CI of the LV1 correlation (1,000 resamples, seed 0) |

**Outputs.** `derived_data/pls_LV1_brain_loading_cortex100.csv`,
`derived_data/pls_LV1_bootstrap_CI.csv`.

**Notes.** The covariate strategy (residualise Y only) follows Winkler et al. (2020) and
Petrican et al. (2022); see Methods 2.5. Expected values: main README section 7.
