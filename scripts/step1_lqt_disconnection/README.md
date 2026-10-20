# Step 1 — Regional white-matter disconnection estimation

**Purpose.** Convert each patient's lesion mask into a parcel-level disconnection probability
profile using the Lesion Quantification Toolkit (LQT) with the HCP-842 normative connectome.

**Inputs.** Lesion masks in MNI space (not shared — patient data); HCP-842 template (ships with LQT).

**Scripts.**
| File | Role |
|---|---|
| `build_lqt_input.py` | assembles the LQT input list (124 subjects, new-mask cohort) |
| `lqt_matlab_full_N124.m` | main LQT run |
| `lqt_matlab_crossval.m` | cross-validation of the parcel damage estimates |
| `psf_sdc_lqt_n124_v2.R` | R driver for the SDC/LQT pipeline |
| `util_get_parcel_damage.m` | parcel-level damage fraction utility |
| `crop_verify_cortex100.py` | crops the 1,000-parcel output to the cortical 100-parcel space |
| `aggregate_lqt_to_matrix.py` | aggregates LQT output into the 124 × 100 matrix |

**Output.** `derived_data/disconnection_matrix_N124_cortex100.csv` (shipped).

**Notes.** LQT requires MATLAB plus the DSI Studio tractography assets distributed with the toolkit.
Steps 2–8 can be run directly from the shipped matrix, so this step is only needed for a full
end-to-end reproduction. Expected values: main README section 7.
