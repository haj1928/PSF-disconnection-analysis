# Step 7 — Specificity controls

**Purpose.** Establish that the disconnection–fatigue pattern is specific to fatigue rather than to
disconnection burden, whole-brain gradients, or single AHBA donors.

**Scripts.**
| File | Role |
|---|---|
| `build_discon_control_vector.py` | builds the disconnection-burden control vector |
| `build_salience_permutations.py` | phenotype permutation sets for spatial nulls |
| `run_sal_specificity_spin_Dataset5.py` | cell-type expression map vs pattern, spin null |
| `run_sal_specificity_spin_v2.py` | same, with explicit gene-set retrieval by GO ID |
| `run_sal_specificity_spin_proteostasis.py` | proteostasis pathway spin test |
| `run_gsea_pathway_spin_official.py` | pathway-level spin survival test |
| `crosscheck_gsea_discon_burden.py` | cross-check of the burden GSEA control |
| `run_lodo_leave_one_donor.py` | leave-one-donor-out robustness |
| `run_lodo_sensitivity_common94.py` | LODO restricted to the parcels common to all donors (final caliber: the 90 common parcels; see `lodo_common90_r_by_donor.csv`) |

**Outputs.** `derived_data/gsea_pathway_spin_survival.csv`,
`derived_data/lodo_donor_robustness_report.txt`, `derived_data/sensitivity_specs.csv`.

**Notes.** Gene sets are always retrieved **by explicit GO/KEGG ID**, never by name prefix
matching — several gene sets share name prefixes in the GMT files. Expected values: main README
section 7.
