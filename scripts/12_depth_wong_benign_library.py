"""Primary prostate tumours (Wong et al.): depth-downsampling experiment, one process per patient.

NOTE (see README, "Known limitations"): this script uses the benign/tumour LIBRARY labels, i.e. it assumes benign-library epithelium is
normal. It builds the lineage offset from other patients' benign-labelled epithelium and scores tumour-library epithelium against a
split-half pooled-expression profile (no DNA exists for these samples). The CNV-label-free analysis is scripts/30_wong_library_agnostic.py.
   python scripts/12_depth_wong_benign_library.py [patient ...|test]"""
import os
os.environ.setdefault('OMP_NUM_THREADS', '4'); os.environ.setdefault('OPENBLAS_NUM_THREADS', '4'); os.environ.setdefault('NUMBA_NUM_THREADS', '4')
import numpy as np, pandas as pd, sys, pickle, time, zlib, multiprocessing as mp
from sccnv import config, datasets, offsets as O
from sccnv import benchmark as bench, experiments, loaders as data_io
from sccnv.core import *

EPI = ['luminal', 'basal', 'club', 'hillock']; STRO = ['endothelial', 'fibroblast', 'smooth_muscle', 'schwann']
IMM = ['b_cells', 't_cells', 'nk_cells', 'macrophage', 'monocyte', 'mast', 'plasma', 'dendritic_plasmacytoid', 'granulocytes']
X, genes, cells, meta = data_io.load_export(str(config.EXPORT / 'wong'))
table, cen = datasets.reference_tables()
anno = meta.anno.values.astype(str); lib = meta.tumor_normal.values.astype(str); pat = meta.patient.values.astype(str)
clus = meta.seurat_clusters.astype(int).astype(str).values
cls = np.array(['epi' if a in EPI else 'stroma' if a in STRO else 'immune' if a in IMM else 'other' for a in anno])
BASE = dict(genes_per_bin=100, weight_power=0.5, ref_combine='median')
cfg_sel = Config(**BASE)
ref_all = np.isin(cls, ['stroma', 'immune'])
sel = select_genes_and_bins(X, genes, table, cen, cfg_sel, np.array([f'{a}|{l}' for a, l in zip(anno, lib)]), ref_all)
print('[global] genes', len(sel['genes']), 'bins', len(sel['bins']), flush=True)

# ---- bin counts per cell (weighted) for lineage-offset estimation
Wd = sp.diags(sel['w'].astype(np.float32))
Yb = np.asarray(((X[:, sel['genes'].col.values] @ Wd) @ sel['M']).todense())

def lp(Ysum, pc=0.5):
    a = Ysum + pc; return np.log2(a / a.sum())

def patient_offset(q):
    m = pat == q
    l_s = lp(Yb[m & (cls == 'stroma')].sum(0)); l_i = lp(Yb[m & (cls == 'immune')].sum(0))
    base = np.median(np.vstack([l_s, l_i]), 0); base = np.log2(2 ** base / (2 ** base).sum())
    epi = lp(Yb[m & (cls == 'epi') & (lib == 'B')].sum(0))
    return epi - base

patients = sorted(set(pat))
n_benign = {q: int(((pat == q) & (cls == 'epi') & (lib == 'B')).sum()) for q in patients}
OFF = {q: patient_offset(q) for q in patients}
def loo_offset(p):
    donors = [q for q in patients if q != p and n_benign[q] >= 150]
    D = np.vstack([OFF[q] for q in donors]); off = np.median(D, 0)
    var = (1.4826 * np.median(np.abs(D - off), 0)) ** 2
    return off, var, donors

CAP = 1500
def process(p, fracs=(1.0, 0.5, 0.25, 0.1), methods=None):
    t00 = time.time(); rng = np.random.default_rng(zlib.crc32(p.encode()))      # stable across runs (str hash() is salted per process)
    idx = np.where(pat == p)[0]; take = []
    for c in ['stroma', 'immune', 'epi']:
        for l in ['B', 'T']:
            ii = idx[(cls[idx] == c) & (lib[idx] == l)]
            take.append(rng.choice(ii, CAP, replace=False) if len(ii) > CAP else ii)
    idx = np.sort(np.concatenate(take)); Xp = X[idx]
    a, l, c_, k = anno[idx], lib[idx], cls[idx], clus[idx]
    ho = (rng.random(len(idx)) < 0.2) & np.isin(c_, ['stroma', 'immune'])
    grp = np.array([f'{a_}|{l_}|c{k_}' + ('_ho' if h else '') for a_, l_, k_, h in zip(a, l, k, ho)], dtype=object)
    ref_masks = {'stroma': (c_ == 'stroma') & ~ho, 'immune': (c_ == 'immune') & ~ho}
    off, ovar, donors = loo_offset(p)
    is_tepi = (c_ == 'epi') & (l == 'T')
    eval_half = is_tepi & (rng.random(len(idx)) < 0.5)          # these cells are scored; they never contribute to the truth
    cell_key = np.where(eval_half, 'T_epi', '')
    norm_mask = (c_ == 'epi') & (l == 'B')
    base = dict(BASE)
    # ---- consensus profile of the *other* half of tumour-library epithelium at full depth (pipeline, metacell-20, offset applied)
    keep = ~eval_half
    res = run(Xp[keep], genes, table, cen, grp[keep], {k: v[keep] for k, v in ref_masks.items()}, Config(**base, metacell_cells=20),
              verbose=False, sel=sel, lineage_offset=off, offset_var=ovar)
    tgroups = [g for g in res.meta_table.group.unique() if g.split('|')[0] in EPI and g.split('|')[1] == 'T']
    r_pool, seg_pool, call_pool = res.pool(tgroups)
    # ---- lineage-offset effect on benign epithelium (negative control) at full depth
    res0 = run(Xp[keep], genes, table, cen, grp[keep], {k: v[keep] for k, v in ref_masks.items()}, Config(**base, metacell_cells=20), verbose=False, sel=sel)
    def fp(r):
        mt = r.meta_table; b = [g for g in mt.group.unique() if g.split('|')[0] in EPI and g.split('|')[1] == 'B']
        rows = np.where(mt.group.isin(b).values)[0]
        return dict(benign_frac_bins_called=float((r.call[rows] != 0).mean()), benign_frac_meta_called=float((r.call[rows] != 0).any(1).mean()),
                    benign_rms_r=float(np.sqrt((r.r[rows] ** 2).mean())), n_meta=len(rows))
    lin = pd.DataFrame([dict(patient=p, lineage_offset='LOO benign-epithelium offset', **fp(res)), dict(patient=p, lineage_offset='none', **fp(res0))])
    # ---- depth experiment (offset applied for pipeline methods)
    truth = {'T_epi': seg_pool}
    df = experiments.run_depth_experiment(Xp, genes, table, cen, grp, ref_masks, sel, cell_key, truth, norm_mask, list(fracs), methods, base, seed=1,
                                          run_kw=dict(lineage_offset=off, offset_var=ovar), eval_kw=dict(alpha=0.025),
                                          log=lambda s: None)
    df['patient'] = p
    print(f'[{p}] done {time.time()-t00:.0f}s  n_cells={len(idx)}  truth bins gain/loss: {(seg_pool>=0.25).sum()}/{(seg_pool<=-0.25).sum()}', flush=True)
    return dict(patient=p, df=df, lin=lin, truth_seg=seg_pool, truth_r=r_pool, truth_call=call_pool, n=len(idx), donors=donors, off=off)

if __name__ == '__main__':
    which = sys.argv[1:] or patients
    if which == ['test']:
        out = [process('ICC2', fracs=(1.0, 0.25), methods={k: v for k, v in experiments.DEFAULT_METHODS.items() if k in ('metacell-20 (pipeline)', 'inferCNV-style (single cell)')})]
    else:
        with mp.get_context('fork').Pool(len(which)) as pool: out = pool.map(process, which)
    pickle.dump(dict(out=out, bins=sel['bins']), open(config.results_path('wong_depth_results.pkl' if which != ['test'] else 'wong_test.pkl'), 'wb'))
    df = pd.concat([o['df'] for o in out]); lin = pd.concat([o['lin'] for o in out])
    print(lin.round(4).to_string(index=False))
    print(df.groupby(['frac', 'method', 'profile'])[['sens_gain_fpr10', 'sens_loss_fpr10', 'auc_gain', 'auc_loss', 'auc_burden', 'pearson']].mean().round(3).to_string())
