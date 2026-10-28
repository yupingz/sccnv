"""Lineage offsets from healthy-donor data: per-bin log2 difference between epithelium and the reference lineage,
estimated within each donor (dataset-level effects cancel). Bins/weights come from the query dataset's selection."""
import numpy as np, pandas as pd, scipy.sparse as sp

EPI = ['luminal', 'basal', 'club', 'hillock']
H_STROMA = ['endothelial', 'fibroblast', 'smooth_muscle', 'schwann']
H_IMMUNE = ['t_cells', 'macrophage', 'granulocytes']

def bin_counts(X, genes, sel):
    """cells x bins weighted counts on the query dataset's bins (genes matched by symbol)."""
    idx = {g: i for i, g in enumerate(genes)}
    miss = [s for s in sel['genes'].symbol if s not in idx]
    assert not miss, f"{len(miss)} bin genes missing, e.g. {miss[:5]}"
    cols = [idx[s] for s in sel['genes'].symbol]
    Wd = sp.diags(sel['w'].astype(np.float32))
    return np.asarray(((X[:, cols] @ Wd) @ sel['M']).todense())

def lp(v, pc=0.5):
    a = v + pc
    return np.log2(a / a.sum())

def ref_lp(Y_sets, combine='median'):
    """log2 reference proportions from one or several reference sets (matches sc_cnv.baseline_logp)."""
    lps = [lp(Y.sum(0)) for Y in Y_sets]
    if len(lps) == 1: l = lps[0]
    elif combine == 'pooled': l = lp(sum(Y.sum(0) for Y in Y_sets))
    else: l = np.median(np.vstack(lps), 0)
    p = 2 ** l
    return np.log2(p / p.sum())

def offset_table(Yb, ctype, unit, epi_types, ref_sets, combine='median'):
    """
    Yb: cells x bins; ctype: cell-type label per cell; unit: donor/library label per cell.
    epi_types: list of cell types forming the epithelial contrast; ref_sets: list of lists of cell types (each its own reference set).
    Returns (units, D) with D[u] = log2 prop(epithelium) - log2 baseline of the reference, per bin.
    """
    ctype = np.asarray(ctype); unit = np.asarray(unit); units = sorted(set(unit)); D = []; used = []
    for u in units:
        m = unit == u
        e = m & np.isin(ctype, epi_types)
        sets = [Yb[m & np.isin(ctype, rs)] for rs in ref_sets]
        if e.sum() < 30 or any(len(s) < 30 for s in sets): continue
        D.append(lp(Yb[e].sum(0)) - ref_lp(sets, combine)); used.append(u)
    return used, (np.vstack(D) if D else np.zeros((0, Yb.shape[1])))

def combine_units(D, min_units=3):
    """median over units + robust between-unit variance (floor: plain variance when <3 units)."""
    off = np.median(D, 0)
    if len(D) >= min_units: var = (1.4826 * np.median(np.abs(D - off), 0)) ** 2
    elif len(D) == 2: var = (D[0] - D[1]) ** 2 / 4
    else: var = np.zeros(D.shape[1])
    return off, var
