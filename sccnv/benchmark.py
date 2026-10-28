"""Benchmark helpers: compare inferred profiles against orthogonal DNA copy-number ground truth."""
import numpy as np, pandas as pd
from scipy.stats import spearmanr, pearsonr
from .core import bin_truth, _mode_level

def mode_center(v, tol=0.1):
    """centre a bin-level log2 profile on its modal level (aneuploid-safe baseline)."""
    ok = ~np.isnan(v); x = v[ok]
    best, lvl = -1, 0.0
    for c in x[::max(1, len(x) // 400)]:
        s = np.abs(x - c) <= tol
        if s.sum() > best: best = s.sum(); lvl = x[s].mean()
    return v - lvl

def truth_classes(t, gain=0.25, loss=-0.25, neutral=0.10):
    c = np.full(len(t), np.nan)
    c[np.abs(t) < neutral] = 0; c[t >= gain] = 1; c[t <= loss] = -1
    return c

def compare(pred_seg, pred_call, truth_bin, label=''):
    """pred_seg: fitted log2 profile; pred_call: -2..2 per bin; truth_bin: mode-centred log2 truth per bin."""
    ok = ~np.isnan(truth_bin)
    t, p = truth_bin[ok], pred_seg[ok]
    out = dict(label=label, n_bins=int(ok.sum()), spearman=spearmanr(p, t)[0], pearson=pearsonr(p, t)[0])
    out['slope'] = np.polyfit(t, p, 1)[0]
    tc = truth_classes(t); pc = np.sign(pred_call[ok])
    m = ~np.isnan(tc)
    for name, v in (('gain', 1), ('loss', -1), ('neutral', 0)):
        tp = ((pc == v) & (tc == v) & m).sum(); fn = ((pc != v) & (tc == v) & m).sum(); fp = ((pc == v) & (tc != v) & m).sum()
        out[f'sens_{name}'] = tp / max(tp + fn, 1); out[f'prec_{name}'] = tp / max(tp + fp, 1); out[f'n_true_{name}'] = int(((tc == v) & m).sum())
    out['balanced_acc'] = np.mean([out['sens_gain'], out['sens_loss'], out['sens_neutral']])
    return out


# ---------------------------------------------------------------- profile-level evaluation (method-agnostic)
from scipy.stats import rankdata

def _dir_auc(P, pos_cols, neg_cols, sign):
    cols = np.r_[pos_cols, neg_cols]; Rk = rankdata(sign * P[:, cols], axis=1)
    npos, nneg = len(pos_cols), len(neg_cols)
    return (Rk[:, :npos].sum(1) - npos * (npos + 1) / 2) / (npos * nneg)

def eval_profiles(P, cell_key, truth_by_key, norm_mask, alpha=0.025, max_cells=300, seed=0, gain=0.25, loss=-0.25, neutral=0.10):
    """
    P: cells x bins (rows with NaN are ignored). cell_key: truth key per cell (None/'' = not evaluated).
    norm_mask: held-out normal cells; per-method thresholds are the alpha / 1-alpha quantiles of their values, so all
    methods are compared at the same per-bin false-positive rate on diploid cells.
    Returns one dict per truth key: sensitivity (gain/loss) at that FPR, FPR on truth-neutral bins of tumour cells,
    directional AUC (gain vs neutral, loss vs neutral), Pearson with truth, and CNV-burden AUC (tumour vs normal).
    """
    rng = np.random.default_rng(seed); valid = ~np.isnan(P).any(1)
    Pn = P[norm_mask & valid]
    qhi, qlo = np.quantile(Pn, 1 - alpha), np.quantile(Pn, alpha)
    burden_n = np.abs(Pn).mean(1)
    out = []
    for key, t in truth_by_key.items():
        idx = np.where((cell_key == key) & valid)[0]
        if len(idx) == 0: continue
        if len(idx) > max_cells: idx = rng.choice(idx, max_cells, replace=False)
        tc = truth_classes(t, gain, loss, neutral); g = np.where(tc == 1)[0]; l = np.where(tc == -1)[0]; n = np.where(tc == 0)[0]
        Pt = P[idx]
        d = dict(key=key, n_cells=len(idx), n_gain=len(g), n_loss=len(l), n_neutral=len(n))
        d['sens_gain'] = (Pt[:, g] > qhi).mean() if len(g) else np.nan
        d['sens_loss'] = (Pt[:, l] < qlo).mean() if len(l) else np.nan
        d['fpr_tumor_neutral'] = ((Pt[:, n] > qhi) | (Pt[:, n] < qlo)).mean() if len(n) else np.nan
        d['auc_gain'] = _dir_auc(Pt, g, n, 1).mean() if len(g) >= 5 and len(n) >= 5 else np.nan
        d['auc_loss'] = _dir_auc(Pt, l, n, -1).mean() if len(l) >= 5 and len(n) >= 5 else np.nan
        if len(n) >= 20:      # sensitivity at a 10% per-cell false-positive rate on truth-neutral bins (threshold-free comparison)
            thr_hi = np.quantile(Pt[:, n], 0.90, axis=1); thr_lo = np.quantile(Pt[:, n], 0.10, axis=1)
            d['sens_gain_fpr10'] = (Pt[:, g] > thr_hi[:, None]).mean() if len(g) else np.nan
            d['sens_loss_fpr10'] = (Pt[:, l] < thr_lo[:, None]).mean() if len(l) else np.nan
        ok = ~np.isnan(t); A = Pt[:, ok] - Pt[:, ok].mean(1, keepdims=True); B = t[ok] - t[ok].mean()
        d['pearson'] = ((A @ B) / np.sqrt((A ** 2).sum(1) * (B ** 2).sum() + 1e-12)).mean()
        b_t = np.abs(Pt).mean(1); allb = np.r_[b_t, burden_n]; Rk = rankdata(allb)
        d['auc_burden'] = (Rk[:len(b_t)].sum() - len(b_t) * (len(b_t) + 1) / 2) / (len(b_t) * len(burden_n))
        out.append(d)
    return out


# ---------------------------------------------------------------- DNA segment truth
def load_wa71_segments(path):
    d = pd.read_csv(path)
    rg = d['Region'].str.extract(r'chr([^:]+):(\d+)-(\d+)')
    return pd.DataFrame(dict(chr=rg[0], start=rg[1].astype(int), end=rg[2].astype(int), cn=d['Copy Number'].astype(float), status=d['Status']))

def segments_to_gene_log2(seg, table, genes):
    """gene-level log2(absolute copy number); genes outside every segment -> NaN. Later segments win at overlaps."""
    out = pd.Series(np.nan, index=genes, dtype=float)
    seg = seg.sort_values(['chr', 'start'])
    for c, g in seg.groupby('chr'):
        tg = table[table.chr == c]
        tg = tg[tg.index.isin(genes)]
        pos = tg.pos.values
        for r in g.itertuples():
            m = (pos >= r.start) & (pos <= r.end)
            out.loc[tg.index[m]] = np.log2(max(r.cn, 0.2))
    return out
