"""Collect the depth-experiment outputs into the summary tables used by the report.
Reads (SCCNV_RESULTS): depth_experiment_celllines.csv, depth_experiment_wa71_{stroma,both}.csv, wong_depth_results.pkl
Writes: depth_experiment_summary.csv, wong_lineage_offset_negative_control.csv, wong_primary_consensus_profiles.csv"""
import pickle
import pandas as pd
from sccnv import config

R = lambda f: config.results_path(f)
MET = ['sens_gain_fpr10', 'sens_loss_fpr10', 'auc_gain', 'auc_loss', 'pearson']
W = pickle.load(open(R('wong_depth_results.pkl'), 'rb'))
wo = pd.concat([o['df'] for o in W['out']])
sets = [('cell lines (CCLE DNA)', pd.read_csv(R('depth_experiment_celllines.csv'))),
        ('WA71 metastasis, stroma ref (tumour DNA)', pd.read_csv(R('depth_experiment_wa71_stroma.csv'))),
        ('WA71 metastasis, stroma+immune ref (tumour DNA)', pd.read_csv(R('depth_experiment_wa71_both.csv'))),
        ('primary tumours, split-half consensus', wo)]
rows = []
for name, d in sets:
    d = d[((d.profile == 'seg') & d.method.str.contains('pipeline|no count')) | (d.profile == 'baseline')]
    g = d.groupby(['frac', 'method'])[MET + ['median_umi']].mean().reset_index(); g.insert(0, 'dataset', name); rows.append(g)
pd.concat(rows).rename(columns={'median_umi': 'median_umi_per_cell'}).to_csv(R('depth_experiment_summary.csv'), index=False)
pd.concat([o['lin'] for o in W['out']]).to_csv(R('wong_lineage_offset_negative_control.csv'), index=False)
b = W['bins']; prof = []
for o in W['out']:
    d = b[['bin', 'chr', 'arm', 'start', 'end', 'n_genes']].copy(); d['patient'] = o['patient']
    d['log2_ratio'] = o['truth_r']; d['log2_ratio_segmented'] = o['truth_seg']; d['call'] = o['truth_call']; prof.append(d)
pd.concat(prof).to_csv(R('wong_primary_consensus_profiles.csv'), index=False)
print('written depth_experiment_summary.csv and Wong tables')
