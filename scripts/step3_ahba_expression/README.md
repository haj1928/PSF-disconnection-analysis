# Step 3 — AHBA expression matrix reconstruction

**Purpose.** Reconstruct the 100-parcel × 15,633-gene cortical expression matrix from the Allen
Human Brain Atlas with `abagen`.

**Inputs.** AHBA microarray data, fetched automatically by abagen (6 donors).

**Scripts.** `gen_ahba_100.py` — abagen pipeline (normalisation, donor selection, parcel
aggregation) aligned to the Schaefer 100Parcels atlas.

**Output.** `PSF_ahba_100_expression.csv` (100 × 15,633). Not shipped: it is large and fully
reproducible from the public AHBA data with this script.

**Notes.** One parcel (`RH_SomMot_2`) has no AHBA coverage and is NaN throughout, leaving 99
parcels for all transcriptomic analyses. Daily reproduction does not require re-running this step.
