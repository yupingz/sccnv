# sccnv

Expression-only copy-number inference from single-cell RNA-seq (10x UMI counts), built around three ideas:

1. **Lineage-aware reference.** Copy number is read against stromal and/or immune cells, not against an average of all cells, and an optional per-bin epithelium-minus-reference offset (estimated inside healthy donors) can be subtracted from query cells.
2. **Count-split metacells.** Cells are grouped within a lineage/sample on one half of each cell's UMIs, and copy number is measured on the other half, so the grouping noise is not tested as signal.
3. **A reference-calibrated null.** Reference metacells (leave-one-out) set the segmentation penalty and the call threshold, instead of fixed cut-offs.

The method is described in [`docs/METHODS.md`](docs/METHODS.md); the tests and their limits are in [`analysis/testing_summary.md`](analysis/testing_summary.md) (source: `analysis/testing_summary.Rmd`).

## Status

Research code from one analysis project. 18 unit tests pass on simulated counts. On real data it has been tested on six prostate cancer cell lines (CCLE DNA copy number; the same object also holds PrEC, the normal reference, and PNT2, which has no DNA truth and is not scored), one lymph-node metastasis with tumour DNA (WA71), three healthy prostate donors, and seven primary tumours without DNA.
No licence has been chosen yet; add one before making the repository public.

## What the tests show

* **Metacells keep sensitivity at low depth.** At the lowest downsampled depth, gain sensitivity at a 10% per-cell false-positive rate rises from single cells to 20-cell metacells in all three sets: cell lines 0.25 to 0.48 (443 UMI/cell), WA71 0.23 to 0.61 (189 UMI/cell), primary tumours 0.20 to 0.45 (504 UMI/cell).
  In WA71 the 20-cell metacells at 189 UMI/cell match single cells at full depth (0.61 against 0.56 for this pipeline's single cells and 0.60 for the inferCNV-style baseline).
* **That sensitivity metric uses the truth to set thresholds.** It equalises false-positive rates across methods but is not a deployable threshold. Each depth uses one thinning realisation; a second realisation moved mean sensitivities by at most 0.04.
* **WA71 agreement with DNA is moderate** (pooled-profile Pearson r = 0.58 across bins). Loss detection is weak (best 0.36), probably because most bins are CN 3 or 4.
* **Healthy-donor lineage offsets transfer only partly.** They reduce spurious signal in diploid epithelium of held-out donors (RMS log2 ratio 0.345 to 0.28; metacells with any call 5.0% to 2.7%), but did not improve detection in WA71. Only 3 donors were available.
* **Primary tumours have no DNA, and identity resolution is not solved.** Pooled amplitudes of healthy and tumour epithelial groups overlap; only strong, patient-specific events separate. The primary-tumour depth "truth" is a pooled expression profile, so those sensitivities do not show detection of true events.

Full numbers, caveats and the reproducibility check are in the report.

## Layout

| Path | Content |
|---|---|
| `sccnv/` | package: `core.py` (pipeline), `offsets.py`, `benchmark.py`, `experiments.py`, `datasets.py`, `loaders.py`, `simulate.py`, `figures.py`, `heatmaps.py` |
| `scripts/` | numbered scripts that regenerate every table and figure (below); `legacy/` holds superseded exploratory scripts |
| `tests/` | simulated-data unit tests |
| `examples/simulated_quickstart.py` | end-to-end example that needs no data files |
| `analysis/` | R Markdown testing report and its rendered `.md` / `.html` |
| `results/reference/` | result tables and figures behind the report, with a provenance note |
| `docs/METHODS.md` | method description |

## Install

```bash
conda env create -f environment.yml && conda activate sccnv
pip install -e .
python examples/simulated_quickstart.py
pytest
```

Minimal use:

```python
from sccnv import Config, run, select_genes_and_bins
cfg = Config(genes_per_bin=100, metacell_cells=20)
sel = select_genes_and_bins(X, gene_symbols, gene_table, centromeres, cfg, groups, reference_mask)
res = run(X, gene_symbols, gene_table, centromeres, groups, {"reference": reference_mask}, cfg, sel=sel)
res.cell_matrix("seg")      # cells x bins
res.segments                # metacell segment table
```

`X` is a cells x genes sparse UMI matrix (integers; round ambient-corrected counts first with `sccnv.core.stochastic_round(X, rng)`), `groups` a label per cell (lineage and cluster), and `reference_mask` marks reference cells.

## Reproducing the results

The input objects are not included. Put these files in `$SCCNV_DATA` (default `data/`): `pca.with.annotation.rds` (cell-line mixture), `wa71.lnmet.rds`, `WA_71_segment_data.csv`, `Wong.merge.rds`, `P16_2018_CellRep_Henry_obj.split.list.rds`. Check the terms of each dataset before sharing them.

```bash
export SCCNV_DATA=/path/to/data SCCNV_WORK=work SCCNV_RESULTS=results/run
python scripts/00_prepare_reference.py        # Ensembl 98 gene table and arm boundaries (downloads the GTF)
python scripts/01_fetch_ccle_truth.py         # CCLE copy number for six prostate lines (cBioPortal API)
bash   scripts/03_export_all.sh               # Seurat -> binary sparse pieces (R + Seurat)
python scripts/10_depth_celllines.py          # add --quick for a reduced run
python scripts/11_depth_wa71.py stroma        # and: both
python scripts/12_depth_wong_benign_library.py
python scripts/20_healthy_offset_qc.py; python scripts/21_healthy_lodo_negative_control.py; python scripts/22_wa71_healthy_offset.py
python scripts/30_wong_library_agnostic.py; python scripts/31_healthy_on_wong_bins.py; python scripts/32_group_level_summaries.py
python scripts/33_summarise_depth.py; python scripts/34_reproducibility_check.py
python scripts/40_make_figures.py
```

Run times on a 48-core machine with four of these jobs running concurrently: `scripts/10` about 15 minutes, `scripts/11` about 19 (stroma) and 22 (both) minutes, `scripts/12` about 16 minutes. `scripts/34` compares fresh outputs with `results/reference/`.
To re-knit the report: `Rscript -e 'rmarkdown::render("analysis/testing_summary.Rmd", output_format = c("github_document","html_document"))'`.

## Known limitations

* Expression-only: ploidy, copy-neutral LOH and balanced gains are not identifiable; truth is mode-centred before comparison.
* One metastasis has DNA truth; six cell lines are a mixture, not a tumour microenvironment.
* The cell-line object has eight labels: six cancer lines (scored), PrEC (reference) and PNT2 (1,269 cells). PNT2 is an immortalised non-tumour line with no CCLE copy-number entry in the study queried (`ccle_broad_2025`). Its cells pass through the pipeline as an unscored query group and take part in gene selection and downsampling, but no metric uses them. Its copy-number state is not established, so it is not assumed diploid.
* `scripts/12_depth_wong_benign_library.py` treats benign-labelled libraries as normal. Benign libraries can contain tumour cells, so the label-free analysis (`scripts/30`) is the one to use for identity questions.
* Metacell profiles are shared by their cells, so events smaller than a metacell are blurred (not tested).
* Count splitting is expected to protect against false positives on theoretical grounds; this was not measured directly.
* No confidence intervals were computed.
