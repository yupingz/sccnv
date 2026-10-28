# results/reference

Result tables and figures behind `analysis/testing_summary.Rmd`. They are small summaries; raw counts and patient-level objects are not included.

| Files | Produced by | Notes |
|---|---|---|
| `depth_experiment_*.csv`, `depth_experiment_summary.csv`, `fig_depth_sensitivity.png`, `fig_wong_primary_profiles.png`, `wong_primary_consensus_profiles.csv`, `wong_lineage_offset_negative_control.csv` | `scripts/10` to `12`, `33`, `40` | cell-line and WA71 depth tables regenerated from this repository (full runs, 2026-10-24) and identical to the first analysis run at every row; the primary-tumour depth results use a stable per-patient seed and differ slightly from the first run (see note below) |
| `wa71_healthy_offset_results.csv`, `healthy_lodo_negative_control.csv` | `scripts/22`, `21` | re-run from this repository; identical to the first run (checked to 4 decimals) |
| `reproducibility_check.csv`, `thinning_replicate_check.csv` | `scripts/34` | differences between the regenerated tables and the first analysis run; effect of a second binomial-thinning realisation |
| `wa71_external_offset_results.csv` | `scripts/legacy/offset_wa71_wong_benign.py` | superseded: offset built from benign-LABELLED libraries |
| `wong_group_level_calls.csv`, `healthy_group_level_calls.csv`, `group_amplitude_healthy_vs_wong.csv`, `group_arm_means.csv`, `wong_epi_metacell_zrms.csv`, `arm_library_check.csv`, `fig_wong_metacell_heatmap_*.png`, `fig_healthy_offset_negative_control.png` | `scripts/30`, `31`, `32`, `40` | regenerated from this repository (see note below) |
| `fig_wa71_dna_vs_expression.png` | `scripts/40` | WA71 vs tumour DNA |
| `pytest_report.txt` | `pytest -v` | 18 simulated-data tests |

Note on seeds: the first analysis seeded per-patient subsampling with Python's `hash(patient)`, which changes between processes. The repository uses a stable
CRC32 of the patient label instead, so the Wong tables regenerated here differ slightly from the first run (for example metacell counts and
the third decimal of some medians) while giving the same picture. Everything else is seeded and reproduces exactly.
