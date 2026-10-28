"""Rebuild the figures from the result files in SCCNV_RESULTS (missing inputs are skipped with a message).
   python scripts/40_make_figures.py"""
import os, pickle
import numpy as np, matplotlib.pyplot as plt
from sccnv import config, datasets, figures, heatmaps, benchmark as bench
from sccnv.core import Config, run, select_genes_and_bins, bin_truth

R = str(config.RESULTS); P = lambda f: os.path.join(R, f)

def have(*fs):
    miss = [f for f in fs if not os.path.exists(P(f))]
    if miss: print('skip (missing):', ', '.join(miss))
    return not miss

if have('depth_experiment_celllines.csv', 'depth_experiment_wa71_stroma.csv', 'wong_depth_results.pkl'):
    figures.fig_depth(figures.load_depth(R), P('fig_depth_sensitivity.png')); print('fig_depth_sensitivity.png')
if have('wong_depth_results.pkl'):
    figures.fig_wong_profiles(P('wong_depth_results.pkl'), P('fig_wong_primary_profiles.png')); print('fig_wong_primary_profiles.png')
if have('healthy_lodo_negative_control.csv', 'group_arm_means.csv'):
    figures.fig_healthy_offset(P('healthy_lodo_negative_control.csv'), P('group_arm_means.csv'), P('fig_healthy_offset_negative_control.png')); print('fig_healthy_offset_negative_control.png')
if have('wong_resolved.pkl'):
    D = pickle.load(open(P('wong_resolved.pkl'), 'rb')); outs = {o['patient']: o for o in D['outs']}
    for which, tag, vmax, title in [('seg', 'seg', 0.6, 'fitted segments'), ('r', 'raw', 0.8, 'unsegmented bin ratios')]:
        fig, _ = heatmaps.wong_heatmap(outs, D['bins'], 'healthy_matched', which, P(f'fig_wong_metacell_heatmap_{tag}.png'), vmax=vmax,
                              title=f'Metacell copy-number profiles, primary prostate tumours ({title})'); plt.close(fig); print(f'fig_wong_metacell_heatmap_{tag}.png')
if os.path.exists(config.EXPORT / 'wa71' / 'meta.csv'):
    S = datasets.wa71(); cfg = Config(target_umi=80000); ref = {'stroma': S['is_str'] & ~S['ho']}
    sel = select_genes_and_bins(S['X'], S['genes'], S['table'], S['cen'], cfg, S['grp'], ref['stroma'])
    res = run(S['X'], S['genes'], S['table'], S['cen'], S['grp'], ref, cfg, verbose=False, sel=sel)
    raw = bin_truth(res.bins, res.genes, bench.segments_to_gene_log2(S['seg'], S['table'], list(res.genes.symbol)))
    figures.fig_wa71(res, datasets.wa71_truth(S, sel), raw, P('fig_wa71_dna_vs_expression.png')); print('fig_wa71_dna_vs_expression.png')
