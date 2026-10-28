"""Figure functions for the depth, WA71, primary-tumour and healthy-offset results (matplotlib)."""
import pickle
import numpy as np, pandas as pd, matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import TwoSlopeNorm
from matplotlib.ticker import NullFormatter
from .plotting import _chr_ticks

C_EXPR, C_DNA, C_SEG = '#1b6ca8', '#555555', '#d1495b'


def style(small=6, medium=7, large=8):
    mpl.rcParams.update({'font.size': large, 'axes.labelsize': large, 'xtick.labelsize': small, 'ytick.labelsize': small, 'legend.fontsize': medium,
                         'axes.spines.top': False, 'axes.spines.right': False})


# ---------------------------------------------------------------------------- depth curves
STYLE = {'single-cell (pipeline)': ('Single cells, same pipeline', '#8a8a8a', '--', 1.2),
         'inferCNV-style (single cell)': ('Single cells, inferCNV-style', '#c4c4c4', '-', 1.2),
         'metacell-20 (pipeline)': ('Metacells, 20 cells', C_EXPR, '-', 2.0),
         'metacell-50 (pipeline)': ('Metacells, 50 cells', '#6fb1de', '-', 1.4),
         'metacell-80k UMI (pipeline)': ('Metacells, 80k UMI (depth-adaptive)', '#0b3c63', '-', 1.4)}
LIMS = {'Prostate cell lines (6 lines, CCLE DNA)': (300, 25000), 'Metastasis WA71 (tumour DNA)': (130, 5500), 'Primary tumours (7 patients)': (330, 8000)}


def load_depth(results_dir, wong_pkl='wong_depth_results.pkl'):
    r = lambda f: pd.read_csv(f'{results_dir}/{f}')
    wo = pd.concat([o['df'] for o in pickle.load(open(f'{results_dir}/{wong_pkl}', 'rb'))['out']])
    return {'Prostate cell lines (6 lines, CCLE DNA)': r('depth_experiment_celllines.csv'),
            'Metastasis WA71 (tumour DNA)': r('depth_experiment_wa71_stroma.csv'), 'Primary tumours (7 patients)': wo}


def _curve(d, metric, method):
    prof = 'baseline' if method.startswith('inferCNV') else 'seg'
    sub = d[(d.method == method) & (d.profile == prof)]
    y = sub.groupby('frac')[metric].mean(); x = sub.groupby('frac').median_umi.mean()   # average over lines/patients per downsampling fraction
    y.index = x.loc[y.index].values
    return y.sort_index()


def fig_depth(data, fname, metrics=('sens_gain_fpr10', 'sens_loss_fpr10'),
              ylabels=('Gain detection\n(sensitivity at 10% FPR)', 'Loss detection\n(sensitivity at 10% FPR)')):
    style()
    fig, axs = plt.subplots(len(metrics), len(data), figsize=(7.2, 2.1 * len(metrics) + 0.5), sharey='row', squeeze=False)
    for j, (title, d) in enumerate(data.items()):
        for i, m in enumerate(metrics):
            ax = axs[i, j]
            for meth, (lab, col, ls, lw) in STYLE.items():
                if meth not in set(d.method): continue
                y = _curve(d, m, meth)
                if y.isna().all(): continue
                ax.plot(y.index, y.values, ls, color=col, lw=lw, marker='o', ms=2.5, label=lab, zorder=3 if 'metacell-20' in meth else 2)
            ax.set_xscale('log'); ax.xaxis.set_minor_formatter(NullFormatter()); ax.margins(x=0.05); ax.set_ylim(0, 1.0)
            ax.set_xticks([200, 500, 1000, 2000, 5000, 10000, 20000]); ax.set_xticklabels(['200', '500', '1k', '2k', '5k', '10k', '20k']); ax.set_xlim(*LIMS[title])
            if j == 0: ax.set_ylabel(ylabels[i])
            if i == 0: ax.set_title(title, loc='left', fontsize=7)
            if i == len(metrics) - 1 and j == 1: ax.set_xlabel('Median UMIs per cell after downsampling')
    h, l = axs[0, 1].get_legend_handles_labels()
    fig.legend(h, l, loc='lower center', ncol=3, frameon=False, fontsize=7, bbox_to_anchor=(0.5, 0.0))
    fig.subplots_adjust(left=0.09, right=0.99, top=0.93, bottom=0.17, wspace=0.12, hspace=0.35)
    fig.savefig(fname, dpi=200); plt.close(fig)


# ---------------------------------------------------------------------------- WA71: metacell heatmap, pooled tumour profile vs DNA, bin-level scatter
def fig_wa71(res, truth_bins, raw_log2_cn, fname, n_tumour=40, n_stroma=25, n_immune=12, seed=0):
    """res: CNVResult from sccnv.run on WA71 with a stroma reference; truth_bins: mode-centred log2 DNA CN per bin; raw_log2_cn: uncentred log2 CN per bin."""
    style(); bins = res.bins; x = np.arange(len(bins)); bnd, mid, lab = _chr_ticks(bins); mt = res.meta_table.copy()
    tg = [g for g in mt.group.unique() if g.startswith('Tumor')]; r_t, seg_t, _ = res.pool(tg); ok = ~np.isnan(truth_bins)
    def panel(g):
        return ('Tumor' if g.startswith('Tumor') else 'Endothelial / SMC (reference)' if g.split('|')[0] in ('Endo_1', 'Endo_2', 'Endo_3', 'SMC')
                else 'Immune (not in reference)')
    mt['panel'] = mt.group.map(panel); rng = np.random.default_rng(seed); rows, labs = [], []
    for pn, mr in [('Tumor', n_tumour), ('Endothelial / SMC (reference)', n_stroma), ('Immune (not in reference)', n_immune)]:
        ii = np.where((mt.panel == pn).values)[0]; ii = np.sort(rng.choice(ii, min(mr, len(ii)), replace=False))
        u = np.linalg.svd(res.seg[ii] - res.seg[ii].mean(0), full_matrices=False)[0][:, 0]; rows.append(ii[np.argsort(u)]); labs.append((pn, len(ii)))
    order = np.concatenate(rows)
    fig = plt.figure(figsize=(7.2, 7.4)); gs = fig.add_gridspec(3, 1, height_ratios=[2.2, 1.3, 1.6], hspace=0.55)
    ax = fig.add_subplot(gs[0]); im = ax.imshow(res.seg[order], aspect='auto', cmap='RdBu_r', norm=TwoSlopeNorm(0, -0.6, 0.6), interpolation='nearest'); pos = 0
    for pn, n in labs:
        if pos > 0: ax.axhline(pos - .5, color='k', lw=.8)
        ax.text(-2, pos + n / 2, pn.replace(' (', '\n('), ha='right', va='center', fontsize=6); pos += n
    for b in bnd: ax.axvline(b, color='k', lw=.3, alpha=.5)
    ax.set_xticks(mid); ax.set_xticklabels(lab, fontsize=6); ax.tick_params(axis='x', length=0); ax.set_yticks([])
    ax.set_title('Tumour metacells carry large copy-number changes; reference-lineage metacells stay flat', loc='left', fontsize=8)
    cb = fig.colorbar(im, ax=ax, fraction=0.018, pad=0.01); cb.set_label('log2 ratio vs reference', fontsize=6); cb.ax.tick_params(labelsize=6)
    ax2 = fig.add_subplot(gs[1])
    ax2.plot(x, np.clip(truth_bins, -1.3, 1.3), color=C_DNA, lw=1.4, drawstyle='steps-mid', label='DNA copy number (log2, modal level = 0)')
    ax2.plot(x, r_t, '.', ms=3, color=C_EXPR, alpha=.7, label='Expression: pooled tumour bins'); ax2.plot(x, seg_t, color=C_SEG, lw=1.6, label='Expression: segmented')
    ax2.axhline(0, color='k', lw=.4)
    for b in bnd: ax2.axvline(b, color='k', lw=.3, alpha=.4)
    ax2.set_ylim(-1.3, 1.3); ax2.set_xlim(-0.5, len(bins) - 0.5); ax2.set_ylabel('log2 ratio'); ax2.set_xticks(mid); ax2.set_xticklabels(lab, fontsize=6); ax2.tick_params(axis='x', length=0)
    ax2.set_title('Pooled expression profile follows the WA71 DNA segments', loc='left', fontsize=8, pad=14)
    ax2.legend(loc='lower left', bbox_to_anchor=(0, 0.98), ncol=3, fontsize=6, frameon=False, handlelength=1.2, columnspacing=1.0, borderaxespad=0)
    ax3 = fig.add_subplot(gs[2]); dfp = pd.DataFrame({'cn': np.round(2 ** raw_log2_cn[ok]), 'r': r_t[ok]})
    for c_, g in dfp.groupby('cn'):
        ax3.scatter(c_ + np.random.default_rng(1).uniform(-.18, .18, len(g)), g.r, s=7, color=C_EXPR, alpha=.6, lw=0); ax3.plot([c_ - .3, c_ + .3], [g.r.median()] * 2, color=C_SEG, lw=1.6)
    ax3.axhline(0, color='k', lw=.4); ax3.set_xlabel('DNA copy number of bin (WA71 segments)'); ax3.set_ylabel('Expression log2 ratio')
    ax3.set_xticks(sorted(dfp.cn.unique())); ax3.set_xticklabels([f'{int(v)}' for v in sorted(dfp.cn.unique())])
    ax3.set_title(f'Expression ratio rises with DNA copy number (Pearson r = {np.corrcoef(r_t[ok], truth_bins[ok])[0, 1]:.2f})', loc='left', fontsize=8)
    for a_, l_ in zip((ax, ax2, ax3), 'abc'): a_.text(-0.17, 1.04, l_, transform=a_.transAxes, fontsize=10, fontweight='bold', va='bottom')
    fig.subplots_adjust(left=0.2, right=0.93, top=0.95, bottom=0.07); fig.savefig(fname, dpi=200); plt.close(fig)


# ---------------------------------------------------------------------------- primary tumours: split-half consensus profiles per patient
def fig_wong_profiles(wong_depth_pkl, fname):
    style(); W = pickle.load(open(wong_depth_pkl, 'rb')); wb = W['bins']; outs = {o['patient']: o for o in W['out']}; pats = sorted(outs)
    bnd, mid, lab = _chr_ticks(wb); xx = np.arange(len(wb))
    fig, axs = plt.subplots(len(pats), 1, figsize=(7.2, 6.6), sharex=True, sharey=True)
    for ax, p in zip(axs, pats):
        o = outs[p]; ax.plot(xx, o['truth_r'], '.', ms=2.5, color=C_EXPR, alpha=.55); ax.plot(xx, o['truth_seg'], color=C_SEG, lw=1.4); ax.axhline(0, color='k', lw=.4)
        for b in bnd: ax.axvline(b, color='k', lw=.3, alpha=.4)
        ax.set_ylabel(p, rotation=0, ha='right', va='center', fontsize=7); ax.set_ylim(-1.1, 0.9); ax.tick_params(axis='y', labelsize=6); ax.yaxis.set_label_coords(-0.06, 0.5)
    axs[-1].set_xticks(mid); axs[-1].set_xticklabels(lab, fontsize=6); axs[-1].tick_params(axis='x', length=0); axs[-1].set_xlim(-0.5, len(wb) - 0.5)
    axs[0].set_title('Primary tumours carry modest, patient-specific copy-number changes (ICC3 is nearly flat)', loc='left', fontsize=8)
    fig.text(0.012, 0.5, 'log2 ratio vs stromal / immune reference', rotation=90, va='center', fontsize=7)
    fig.subplots_adjust(left=0.12, right=0.99, top=0.94, bottom=0.05, hspace=0.12); fig.savefig(fname, dpi=200); plt.close(fig)


# ---------------------------------------------------------------------------- healthy offsets: leave-one-donor-out control and overlap with tumour groups
def fig_healthy_offset(lodo_csv, arm_means_csv, fname):
    style(); lodo = pd.read_csv(lodo_csv); AM = pd.read_csv(arm_means_csv, index_col=0)
    Hm = AM[AM.index.str.startswith('Healthy')]; Wm = AM[AM.index.str.startswith('Wong')]
    fig, axs = plt.subplots(1, 2, figsize=(7.2, 3.0), gridspec_kw=dict(width_ratios=[1, 1.15]))
    ax = axs[0]
    vs = ['no offset', 'other donors, all-epithelium offset', 'other donors, matched-subtype offset', 'control: all-epithelium offset, bins permuted']
    labs = ['None', 'Healthy\nall epi.', 'Healthy\nmatched', 'Permuted\ncontrol']; cols = ['#9a9a9a', C_EXPR, '#0b3c63', C_SEG]
    for j, (v, c) in enumerate(zip(vs, cols)):
        y = lodo[lodo.variant == v].set_index('donor').rms_r_epi
        ax.bar(j, y.mean(), color=c, alpha=.85, width=.65); ax.scatter(np.full(len(y), j) + np.linspace(-.12, .12, len(y)), y.values, color='k', s=8, zorder=3)
    ax.set_xticks(range(4)); ax.set_xticklabels(labs, fontsize=6); ax.set_ylabel('RMS log2 ratio of healthy epithelial\nmetacells vs own stroma/immune')
    ax.set_title('Offsets from other donors cut spurious signal\nin diploid epithelium (3 donors, dots)', loc='left', fontsize=7.5)
    ax = axs[1]; hv_, wv_ = Hm.abs().max(1), Wm.abs().max(1); rng = np.random.default_rng(0)
    ax.scatter(rng.uniform(-.15, .15, len(hv_)), hv_.values, s=12, color=C_EXPR, alpha=.7, lw=0)
    ax.scatter(1 + rng.uniform(-.15, .15, len(wv_)), wv_.values, s=12, color=C_SEG, alpha=.8, lw=0)
    for j, v in enumerate([hv_, wv_]): ax.plot([j - .25, j + .25], [v.median()] * 2, color='k', lw=1.2)
    ax.set_xticks([0, 1]); ax.set_xticklabels([f'Healthy\n(3 donors, n={len(hv_)})', f'Primary tumours\n(7 patients, n={len(wv_)})'], fontsize=6.5); ax.set_xlim(-.6, 1.6)
    ax.set_ylabel('Largest arm-level |log2 ratio|\nof pooled group profile'); ax.set_title('Pooled profiles overlap: healthy epithelium\nhas a sizeable non-CNV floor', loc='left', fontsize=7.5)
    fig.subplots_adjust(left=0.1, right=0.99, top=0.82, bottom=0.2, wspace=0.35); fig.savefig(fname, dpi=200); plt.close(fig)
