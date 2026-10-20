# Step 6 — PET receptor atlas mapping

**Purpose.** Correlate the LV1 disconnection pattern with 21 normative receptor/transporter
atlases, calibrate significance with spherical-spin nulls, and separate the independent
contribution of each atlas.

**Inputs.** `derived_data/pls_LV1_brain_loading_cortex100.csv`; 21 atlases from `neuromaps`
(Markello et al., 2022).

**Scripts.**
| File | Role |
|---|---|
| `run_receptor_discon_cortex100_residY.py` | marginal spatial correlation + spin null + within-family BH |
| `run_spin_partial_receptor.py` | partial correlations controlling the other atlases, and controlling the atlas PC1 |
| `run_cholinergic_family.py` | consistency within the three cholinergic atlases |
| `run_multiseed_spin_21atlas.py` | 20 independent rotation seeds (1,000 spins each) |

**Outputs.** `derived_data/receptor_mapping_21atlases.csv`,
`derived_data/receptor_partial_correlation.csv`,
`derived_data/receptor_partial_control_PC1.csv`,
`derived_data/receptor_cholinergic_family.csv`,
`derived_data/receptor_21x21_correlation.csv`,
`derived_data/multiseed_receptor_spin_p_matrix.csv`.

**Notes.** Three of 21 atlases survive within-family BH correction; only α4β2 remains after
controlling the others, and it is unchanged after controlling the atlas PC1. The 21 atlases are
collinear, so the "single independent receptor signal" claim rests on the partial correlations.
Expected values: main README section 7.
