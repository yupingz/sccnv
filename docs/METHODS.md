# Method

`sccnv.run` infers relative copy number from a cells x genes UMI count matrix, using only expression. All log-ratios are log2.
Code: `sccnv/core.py` (pipeline), `sccnv/offsets.py` (lineage offsets), `sccnv/benchmark.py` (evaluation).

## Steps

1. **Groups.** Labels per cell (curated clusters or annotation; `cluster_cells()` gives Leiden clusters on highly variable genes when none exist).
   Metacells are never built across groups, so a metacell stays within one lineage and one sample. Centering is never done per group: a clonal group would absorb its own
   shared copy-number changes into its baseline.
2. **Genes and bins.** Genes must be detected in at least 1% of cells and must not match the blacklist (mitochondrial, ribosomal, histone, immunoglobulin / T-cell receptor and HLA genes).
   A gene's weight is capped so that no gene exceeds 15% of a bin in any group, and weights are tempered (`weight_power = 0.5`: weight proportional to the gene's expression relative to its bin, raised to -0.5). Genes are ordered by midpoint along the genome and split
   into bins of about 100 genes inside each chromosome arm (`genes_per_bin`); arms with fewer than 20 expressed genes are dropped.
3. **Count splitting.** Each cell's UMIs are thinned binomially into fold A (used only to build metacells) and fold B (used only for the bin counts that are tested).
   This avoids testing the same noise that was used to group the cells (`thin = True`).
4. **Metacells.** Within each group, the fold-A bin profiles (log-ratio to the group mean) are reduced to their top 10 principal components and clustered by k-means
   with a capacity cap, so metacells have similar sizes. Size is a cell number (`metacell_cells`, default 20) or a UMI target for the signal fold (`target_umi`; depth-adaptive).
   Groups smaller than `min_metacell_cells` (8) are dropped. If the smallest reference set ends up with fewer than `min_ref_metacells` (12) metacells, the metacell size is halved (up to four times).
5. **Reference.** Any number of reference lineages (stroma, immune, normal epithelium). The baseline is pooled or the per-bin median across reference sets (`ref_combine`).
   Reference metacells are scored leave-one-out, so they provide an empirical null for their own lineage.
6. **Null.** Per-bin variance = Poisson term from the bin counts + extra per-bin variance estimated from the reference residuals + an arm-level correlated component
   (`lambda`, `tau`, fitted on reference window means) + a group-specific floor for non-copy-number expression mismatch between a test group and the reference + the variance of an optional lineage offset.
7. **Segmentation and calls.** Penalised dynamic programming within each arm, after centering each profile on its length-weighted modal segment level (so aneuploid genomes are centered on their own baseline).
   Two settings are calibrated on the leave-one-out reference metacells: the segmentation penalty is the smallest value in `penalty_grid` at which the reference profiles average at most 1.05 segments per arm,
   and the |z| threshold is the smallest value at or above `min_z` (4.0) at which no more than `fp_target` (2%) of reference metacells carry any call.
   A segment is called when |log2 ratio| >= `min_effect` (0.10), it spans at least `min_seg_bins` (2) bins and |z| passes the threshold; the call is doubled (+-2) when |log2 ratio| >= 0.5.
8. **Pooled profiles.** `CNVResult.pool(groups)` combines metacells by precision weighting and re-segments (used for group-level profiles).
9. **Lineage offset (optional).** `lineage_offset` and `offset_var` (an array or `{group label: array}`) subtract a per-bin epithelium-minus-reference offset from query metacells only.
   `sccnv.offsets` estimates it inside each healthy donor (so donor-wide effects cancel), combines donors by the median and carries a robust between-donor variance.

## Output (`CNVResult`)

| Field | Content |
|---|---|
| `bins`, `genes` | bin table (chromosome, arm, start, end, genes) and the kept genes with bin and weight |
| `meta_table` | one row per metacell: group, cells, UMIs, reference flag, group floor |
| `cell_meta` | metacell index of every cell (-1 = dropped) |
| `r`, `sd`, `z` | centered log2 ratio, null SD and z-score, metacells x bins |
| `seg`, `call` | fitted segment value and call (-1, 0, +1; 2 for high-level) per metacell and bin |
| `segments`, `group_segments` | segment tables for metacells and for groups |
| `cell_matrix(which)` | cells x bins matrix of `r`, `seg` or `call` (each cell inherits its metacell's profile) |
| `pool(groups)` | pooled `(r, seg, call)` over the metacells of the given groups |

## Evaluation conventions (`sccnv/benchmark.py`)

* Truth is mode-centered, because expression cannot give absolute ploidy. Gain / loss / neutral bins are log2 >= +0.25 / <= -0.25 / |.| < 0.10 relative to the modal level (arguments `gain`, `loss`, `neutral` of `eval_profiles`).
* `sens_gain_fpr10`, `sens_loss_fpr10`: for each cell, the call threshold is the 90th (gain) or 10th (loss) percentile of that cell's own truth-neutral bins.
  This uses the truth to equalise false-positive rates across methods; it is not a deployable threshold.
* `auc_gain`, `auc_loss`: directional AUC of gain (loss) bins against neutral bins. `pearson`: correlation of the cell's profile with the truth over bins.
* Downsampling is binomial thinning of every cell's UMIs; one thinning realisation per depth.
* Cells of one metacell share the metacell profile.
