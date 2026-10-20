# PSF_repro_PET — spin-test code and fsaverage5 surface assets

This companion folder ships the **code and surface assets** required by the spherical-spin
permutation tests in:

- `scripts/step4_gene_ranking_celltype/` (cell-type specificity contrasts)
- `scripts/step6_pet_receptor/` (receptor spin nulls)
- `scripts/step7_specificity_controls/` (salience permutations, pathway spin tests)
- the two surface-rendering figure scripts (`fig2_behavioral_pls_20260916.py` and
  `figS_all_20260916.py`, through the shared `ANNOT`/aparc paths)

The scripts locate this folder through a `PET`-style constant, e.g. `PET = BASE + '/PSF_repro_PET'`.
Point `BASE` at the repository root (or `PET` directly at this folder) to resolve it.

## Contents

| File | Role |
|---|---|
| `spin_test_numpy.py` | vtk-free spin-rotation implementation (`rotate_parcellation`, `perm_sphere_p`), verbatim from the ENIGMA toolbox (MICA-MNI), with a fixed `seed=0` so the 1,000-rotation null is deterministic across runs |
| `assets/fsa5_sphere_lh.gii`, `assets/fsa5_sphere_rh.gii` | fsaverage5 sphere coordinates (10,242 vertices per hemisphere) |
| `assets/fsa5_lh_aparc.annot`, `assets/fsa5_rh_aparc.annot` | fsaverage5 Desikan-Killiany parcellation |
| `assets/schaefer100_fsa5/lh.Schaefer2018_100Parcels_7Networks_order.annot`, `assets/schaefer100_fsa5/rh.…` | Schaefer et al. (2018) 100-parcel / 7-network annotation on fsaverage5 (Thomas Yeo lab distribution) |

## Provenance

- The three rotation functions in `spin_test_numpy.py` are verbatim copies of
  `enigmatoolbox.permutation_testing` (MICA-MNI/ENIGMA), reused under its release terms; see
  Alexander-Bloch et al. 2018, *NeuroImage* 178:540-51 and Vasa et al. 2018, *Cereb Cortex*
  28(1):281-294.
- The `.gii` / `.annot` assets are the official ENIGMA-toolbox distribution (MICA-MNI/ENIGMA),
  identical to the files used to produce the published spin *p*-values.
- The full PRAND/PET reproduction bundle (cached `neuromaps` atlas collection, ~148 MB) is not
  shipped here; only the code and surface assets needed to re-run the spin tests.

## Licensing of third-party components

The MIT licence at the repository root covers this repository's own code and documentation. It
does **not** relicense the third-party components below, which are redistributed unmodified under
their own terms:

| Component | Source | Terms |
|---|---|---|
| `spin_test_numpy.py` (rotation functions) | MICA-MNI/ENIGMA toolbox (`enigmatoolbox.permutation_testing`), v2.x | BSD-3-Clause; please cite Larivière et al. 2021, *Nat Methods* 18:698–700 |
| `assets/fsa5_sphere_*.gii`, `assets/fsa5_*_aparc.annot` | ENIGMA toolbox distribution of the FreeSurfer fsaverage5 surface | BSD-3-Clause as distributed by ENIGMA; fsaverage surfaces were originally generated with FreeSurfer (Fischl 2012), whose own terms apply to FreeSurfer-derived data |
| `assets/schaefer100_fsa5/*.annot` | Schaefer et al. 2018, 100-parcel / 7-network annotation on fsaverage5 (Thomas Yeo lab distribution) | Cite Schaefer et al. 2018, *Cereb Cortex* 28:3095–3114 |

If you redistribute this folder, keep this table and the upstream citations with it. If you prefer
not to rely on the bundled surfaces, `nilearn.datasets.fetch_surf_fsaverage('fsaverage5')` provides
equivalent fsaverage5 geometry and can be substituted by pointing the `ANNOT`/`sphere` paths at the
nilearn cache.
