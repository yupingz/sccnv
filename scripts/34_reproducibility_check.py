"""Regression check: compare a fresh run (SCCNV_RESULTS, default results/run) with the committed reference tables (results/reference).

    python scripts/34_reproducibility_check.py [reference_dir]

Writes reproducibility_check.csv to SCCNV_RESULTS: per file, rows matched and the largest / mean absolute difference per metric.
If *_quick.csv depth outputs exist (scripts/10 and 11 with --quick), also writes thinning_replicate_check.csv: the same experiment with a different
binomial-thinning realisation at depths shared with the full run, which shows how much a single thinning draw moves the mean sensitivities."""
import sys, os
import numpy as np, pandas as pd
from sccnv import config

REF = sys.argv[1] if len(sys.argv) > 1 else 'results/reference'
RUN = str(config.RESULTS)
DEPTH_KEYS = ['frac', 'method', 'profile', 'key']; DEPTH_MET = ['sens_gain_fpr10', 'sens_loss_fpr10', 'auc_gain', 'auc_loss', 'pearson']
CHECKS = [('depth_experiment_celllines.csv', DEPTH_KEYS, DEPTH_MET),
          ('depth_experiment_wa71_stroma.csv', DEPTH_KEYS, DEPTH_MET),
          ('depth_experiment_wa71_both.csv', DEPTH_KEYS, DEPTH_MET),
          ('depth_experiment_summary.csv', ['dataset', 'frac', 'method'], DEPTH_MET),
          ('wa71_healthy_offset_results.csv', ['trial', 'profile'], ['sens_gain_fpr10', 'sens_loss_fpr10', 'auc_gain', 'auc_loss', 'pearson', 'pooled_pearson']),
          ('healthy_lodo_negative_control.csv', ['donor', 'variant'], ['frac_meta_any_call', 'rms_r_epi', 'floor_epi']),
          ('wong_lineage_offset_negative_control.csv', ['patient', 'lineage_offset'], ['benign_frac_meta_called', 'benign_rms_r', 'n_meta']),
          ('group_amplitude_healthy_vs_wong.csv', ['who', 'label', 'group'], ['seg_rms', 'max_abs_seg', 'n_meta'])]
rows = []
for f, keys, mets in CHECKS:
    a, b = os.path.join(REF, f), os.path.join(RUN, f)
    if not (os.path.exists(a) and os.path.exists(b)): print('skip', f); continue
    A, B = pd.read_csv(a), pd.read_csv(b); m = A.merge(B, on=keys, suffixes=('_ref', '_run'))
    for c in mets:
        if c + '_ref' not in m: continue
        d = (m[c + '_ref'] - m[c + '_run']).abs()
        rows.append(dict(file=f, metric=c, rows_ref=len(A), rows_run=len(B), rows_matched=len(m), max_abs_diff=d.max(), mean_abs_diff=d.mean()))
out = pd.DataFrame(rows); out.to_csv(config.results_path('reproducibility_check.csv'), index=False)
print(out.round(4).to_string(index=False))

reps = []
for tag, full in [('celllines', 'depth_experiment_celllines.csv'), ('wa71_stroma', 'depth_experiment_wa71_stroma.csv')]:
    q, f = os.path.join(RUN, f'depth_experiment_{tag}_quick.csv'), os.path.join(RUN, full)
    if not (os.path.exists(q) and os.path.exists(f)): continue
    Q, F = pd.read_csv(q), pd.read_csv(f); keep = ['single-cell (pipeline)', 'metacell-20 (pipeline)']
    for D in (Q, F): D.drop(D[~(D.method.isin(keep) & (D.profile == 'seg'))].index, inplace=True)
    gq = Q.groupby(['frac', 'method'])[['sens_gain_fpr10', 'sens_loss_fpr10']].mean(); gf = F.groupby(['frac', 'method'])[['sens_gain_fpr10', 'sens_loss_fpr10']].mean()
    j = gq.join(gf, lsuffix='_replicate', rsuffix='_full', how='inner').reset_index(); j.insert(0, 'dataset', tag)
    j['abs_diff_gain'] = (j.sens_gain_fpr10_replicate - j.sens_gain_fpr10_full).abs(); j['abs_diff_loss'] = (j.sens_loss_fpr10_replicate - j.sens_loss_fpr10_full).abs(); reps.append(j)
if reps:
    r = pd.concat(reps); r.to_csv(config.results_path('thinning_replicate_check.csv'), index=False); print(r.round(3).to_string(index=False))
