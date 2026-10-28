"""
sc_cnv: expression-only copy-number inference from scRNA-seq counts.

Pipeline (see README):
  1. cluster cells (or use supplied labels)                          -> groups
  2. bin genes along the genome (fixed gene count, per chromosome arm)
  3. metacells within each group, built on a *thinned* count fold (A) and
     quantified on the complementary fold (B)  -> no double dipping
  4. reference lineages (stromal / immune / normal epithelium ...) give a
     bin-level baseline + empirical null; optional lineage offset
  5. depth-aware null SD, per-metacell centring, penalised segmentation,
     gain / loss calls calibrated on held-out reference metacells.

Conventions: counts are cells x genes (scipy CSR); all log-ratios are log2.
"""
from __future__ import annotations
import re, json
from dataclasses import dataclass, field, asdict, replace
import numpy as np, pandas as pd, scipy.sparse as sp
from numba import njit

LN2 = np.log(2.0)

# --------------------------------------------------------------------------- config
@dataclass
class Config:
    genes_per_bin: int = 100
    min_arm_genes: int = 20            # arms with fewer expressed genes are dropped
    min_cells_frac: float = 0.01       # gene must be detected in >= this fraction of cells
    include_x: bool = False
    blacklist_regex: str = r'^(MT-|RPL\d|RPS\d|HIST|H1-|H2A|H2B|H3-|H4-|IG[HKL][VJC]|TR[ABDG][VJC]|HLA-)'
    max_gene_share: float = 0.15       # cap on one gene's share of a bin (any group)
    # metacells
    metacell_cells: int = 20           # target cells per metacell (1 = single-cell mode)
    target_umi: float | None = None    # alternatively: target signal-fold UMIs per metacell
    min_metacell_cells: int = 8
    min_ref_metacells: int = 12        # shrink metacells until the reference has at least this many
    n_pcs: int = 10
    thin: bool = True                  # count splitting: structure fold A / signal fold B
    split_frac: float = 0.5
    # log-ratio / null
    pseudocount: float = 0.5
    ref_combine: str = 'pooled'        # 'pooled' | 'median' (across reference sets)
    var_shrink: float = 0.5            # shrink per-bin extra variance toward its median
    n_ref_pcs: int = 0                 # remove top PCs of reference-metacell variation (cross-fitted); 0 = off
    weight_power: float = 0.5          # gene weights ~ q^-power (0 = plain count sum, 1 = mean of per-gene ratios)
    floor: bool = True                 # group-specific non-CN expression-mismatch variance floor
    # segmentation / calling
    seg_penalty: float | None = None   # None -> calibrated on reference metacells
    penalty_grid: tuple = (6, 9, 12, 16, 22, 30, 40, 60)
    fp_target: float = 0.02            # max fraction of reference metacells with any call
    min_seg_bins: int = 2
    min_effect: float = 0.10           # |log2 ratio| of a segment to be called
    min_z: float = 4.0
    high_effect: float = 0.5           # flag high-level gain / deep loss
    recenter_iters: int = 2
    seed: int = 0


# --------------------------------------------------------------------------- genome helpers
def load_gene_table(path):
    t = pd.read_csv(path, dtype={'chr': str})
    t = t.drop_duplicates('symbol').copy()
    t['pos'] = ((t.start + t.end) // 2).astype(int)
    return t.set_index('symbol')


def load_centromeres(path):
    c = pd.read_csv(path, dtype={'chr': str}).set_index('chr')
    return {k: (r.p_end, r.q_start) for k, r in c.iterrows()}


def match_genes(names, table):
    """exact symbol match, then fall back to stripping a '.N' duplicate suffix."""
    out = []
    for n in names:
        if n in table.index: out.append(n)
        else:
            b = re.sub(r'\.\d+$', '', n)
            out.append(b if (b != n and b in table.index) else None)
    return out


# --------------------------------------------------------------------------- step 0/2: gene selection + bins
def select_genes_and_bins(X, gene_names, table, cen, cfg: Config, groups=None, ref_mask=None):
    """Return dict with kept genes (ordered along genome), bin table, gene->bin matrix M, gene weights w and
    per-bin Poisson constants a (relative variance of a weighted bin sum = a / total counts)."""
    n_cells = X.shape[0]
    det = np.asarray((X > 0).sum(0)).ravel() / n_cells
    mapped = match_genes(gene_names, table)
    bl = re.compile(cfg.blacklist_regex)
    rows = []
    chroms = [str(i) for i in range(1, 23)] + (['X'] if cfg.include_x else [])
    for j, (g, m) in enumerate(zip(gene_names, mapped)):
        if m is None or det[j] < cfg.min_cells_frac or bl.match(g): continue
        r = table.loc[m]
        if r.chr not in chroms: continue
        rows.append((j, g, r.chr, int(r.pos)))
    gdf = pd.DataFrame(rows, columns=['col', 'symbol', 'chr', 'pos'])
    gdf['chr_i'] = gdf.chr.map({c: i for i, c in enumerate(chroms)})
    gdf['arm'] = [('p' if p <= cen[c][0] else 'q') for c, p in zip(gdf.chr, gdf.pos)]
    gdf = gdf.sort_values(['chr_i', 'pos']).reset_index(drop=True)
    bins, gbin, b = [], np.full(len(gdf), -1), 0
    for (c, arm), sub in gdf.groupby(['chr_i', 'arm'], sort=False):
        if len(sub) < cfg.min_arm_genes: continue
        nb = max(1, int(round(len(sub) / cfg.genes_per_bin)))
        for part in np.array_split(sub.index.values, nb):
            gbin[part] = b
            bins.append(dict(bin=b, chr=sub.chr.iloc[0], arm=arm, start=int(gdf.pos[part[0]]), end=int(gdf.pos[part[-1]]),
                             n_genes=len(part), chr_i=int(sub.chr_i.iloc[0])))
            b += 1
    gdf['bin'] = gbin
    gdf = gdf[gdf.bin >= 0].reset_index(drop=True)
    bins = pd.DataFrame(bins); gb = gdf.bin.values
    M = sp.csr_matrix((np.ones(len(gdf)), (np.arange(len(gdf)), gb)), shape=(len(gdf), len(bins)))
    Xs = X[:, gdf.col.values]
    labs = pd.Series(groups).astype(str).values if groups is not None else np.zeros(n_cells, str)
    w = np.ones(len(gdf)); share = np.zeros(len(gdf)); qmix = np.zeros(len(gdf)); ug = np.unique(labs)
    for g in ug:
        P = np.asarray(Xs[labs == g].sum(0)).ravel().astype(float)
        T = M.T @ P
        share = np.maximum(share, P / np.maximum(T[gb], 1e-9))      # gene's share of its bin in the group where it is largest
        qmix += P / max(P.sum(), 1) / len(ug)                        # group-averaged expression proportion
    if cfg.max_gene_share < 1:                                       # cap dominant genes
        w = np.minimum(1.0, cfg.max_gene_share / np.maximum(share, 1e-9))
    nbg = np.asarray(M.sum(0)).ravel()
    if cfg.weight_power > 0:                                         # temper high-expression genes: w ~ q^-gamma
        Qb = M.T @ qmix
        rel = qmix / np.maximum((Qb / nbg)[gb], 1e-12)
        w = w * np.clip(rel, 0.05, 20.0) ** (-cfg.weight_power)
    qref = np.asarray(Xs[ref_mask].sum(0)).ravel().astype(float) if ref_mask is not None else qmix * 1.0
    qref = qref / qref.sum()
    Qb = M.T @ qref
    sq = M.T @ (w * qref)
    w = w * (Qb / np.maximum(sq, 1e-12))[gb]                         # keep counts of a bin on their original scale
    a = (M.T @ (w ** 2 * qref)) / np.maximum(M.T @ (w * qref), 1e-12) ** 2
    gdf['w'] = w; gdf['q_ref'] = qref
    return dict(genes=gdf, bins=bins, M=M, w=w, a=np.asarray(a).ravel())


# --------------------------------------------------------------------------- step 3: count splitting + metacells
def split_counts(X, frac, rng):
    """binomial thinning of a count matrix: A ~ Bin(x, frac), B = x - A (independent under Poisson)."""
    X = X.tocsr().astype(np.float32)
    a = rng.binomial(X.data.astype(np.int64), frac).astype(np.float32)
    A = X.copy(); A.data = a
    Bm = X.copy(); Bm.data = X.data - a
    A.eliminate_zeros(); Bm.eliminate_zeros()
    return A, Bm


def _balanced_assign(D, cap):
    """greedy capacity-constrained assignment of rows to columns of distance matrix D."""
    n, k = D.shape
    order = np.argsort(D, axis=None)
    assign = np.full(n, -1); load = np.zeros(k, int)
    for o in order:
        i, c = divmod(o, k)
        if assign[i] >= 0 or load[c] >= cap: continue
        assign[i] = c; load[c] += 1
        if (assign >= 0).all(): break
    return assign


def build_metacells(A_w, Bsig_w, M, groups, cfg: Config, rng):
    """A_w: structure fold (cells x genes, weighted); returns metacell id per cell (-1 = dropped)."""
    from sklearn.cluster import KMeans
    n = A_w.shape[0]
    meta = np.full(n, -1); mg = []
    groups = np.asarray(groups).astype(str)
    umi_B = np.asarray(Bsig_w.sum(1)).ravel()
    nxt = 0
    for g in pd.unique(groups):
        idx = np.where(groups == g)[0]
        if len(idx) < cfg.min_metacell_cells: continue
        size = cfg.metacell_cells
        if cfg.target_umi: size = int(np.ceil(cfg.target_umi / max(np.median(umi_B[idx]), 1)))
        if size > 1: size = max(size, cfg.min_metacell_cells)
        k = max(1, int(round(len(idx) / max(size, 1))))
        if size <= 1:                                   # single-cell mode
            for i in idx: meta[i] = nxt; mg.append(g); nxt += 1
            continue
        if k == 1:
            meta[idx] = nxt; mg.append(g); nxt += 1; continue
        Ab = (A_w[idx] @ M).toarray()
        tot = Ab.sum(1, keepdims=True)
        prop = (Ab + cfg.pseudocount) / (tot + cfg.pseudocount * Ab.shape[1])
        L = np.log2(prop / prop.mean(0, keepdims=True))   # log-ratio to group mean (structure fold only)
        L = (L - L.mean(0)) / (L.std(0) + 1e-6)
        nc = min(cfg.n_pcs, L.shape[1] - 1, L.shape[0] - 1)
        U, S, Vt = np.linalg.svd(L, full_matrices=False)
        P = U[:, :nc] * S[:nc]
        km = KMeans(k, n_init=3, random_state=cfg.seed).fit(P)
        D = ((P[:, None, :] - km.cluster_centers_[None]) ** 2).sum(-1)
        cap = int(np.ceil(len(idx) / k))
        asg = _balanced_assign(D, cap)
        for c in range(k):
            sel = idx[asg == c]
            if len(sel) < cfg.min_metacell_cells: continue
            meta[sel] = nxt; mg.append(g); nxt += 1
    return meta, np.array(mg)


def aggregate(Xw, meta, M):
    ok = meta >= 0
    n_meta = meta.max() + 1
    Mc = sp.csr_matrix((np.ones(ok.sum()), (meta[ok], np.where(ok)[0])), shape=(n_meta, len(meta)))
    return np.asarray(((Mc @ Xw) @ M).todense()), np.asarray(Mc.sum(1)).ravel()


# --------------------------------------------------------------------------- step 4/5: log-ratios, null, segmentation
def baseline_logp(Y, ref_sets, cfg, loo_row=None):
    """log2 baseline bin proportions from reference metacell sets. loo_row: reference metacell to leave out."""
    lps = []
    for rows in ref_sets:
        rows = [r for r in rows if r != loo_row]
        S = Y[rows].sum(0) + cfg.pseudocount
        lps.append(np.log2(S / S.sum()))
    if len(lps) == 1 or cfg.ref_combine == 'pooled':
        if len(lps) == 1: lp = lps[0]
        else:
            rows = [r for rr in ref_sets for r in rr if r != loo_row]
            S = Y[rows].sum(0) + cfg.pseudocount; lp = np.log2(S / S.sum())
    else:
        lp = np.median(np.vstack(lps), 0)
    p = 2 ** lp; p /= p.sum()
    return p


def log_ratio(Y, p, cfg):
    N = Y.sum(1, keepdims=True)
    mu = N * p[None]
    return np.log2((Y + cfg.pseudocount) / (mu + cfg.pseudocount)), N.ravel(), mu


@njit(cache=True)
def _opt_partition(S0, S1, S2, beta, min_len):
    n = len(S0) - 1
    F = np.full(n + 1, np.inf); F[0] = -beta
    last = np.zeros(n + 1, np.int64)
    for j in range(min_len, n + 1):
        best = np.inf; bi = -1
        for i in range(0, j - min_len + 1):
            if not np.isfinite(F[i]): continue
            w = S0[j] - S0[i]
            if w <= 0: continue
            c = (S2[j] - S2[i]) - (S1[j] - S1[i]) ** 2 / w
            v = F[i] + c + beta
            if v < best: best = v; bi = i
        F[j] = best; last[j] = bi
    cps = []
    j = n
    while j > 0:
        cps.append(j); j = last[j]
    cps.append(0)
    return np.array(cps[::-1])


def segment_profile(r, w, arm_id, beta, min_len):
    """piecewise-constant fit per arm. r, w: (n_bins,). returns breakpoints list [(a,b,mean,se)]."""
    segs = []
    for a in pd.unique(arm_id):
        ix = np.where(arm_id == a)[0]
        rr, ww = r[ix], w[ix]
        S0 = np.concatenate([[0], np.cumsum(ww)]); S1 = np.concatenate([[0], np.cumsum(ww * rr)])
        S2 = np.concatenate([[0], np.cumsum(ww * rr * rr)])
        ml = min(min_len, len(ix))
        cps = _opt_partition(S0, S1, S2, float(beta), int(ml))
        for s, e in zip(cps[:-1], cps[1:]):
            W = S0[e] - S0[s]
            segs.append((ix[s], ix[e - 1] + 1, (S1[e] - S1[s]) / W, 1 / np.sqrt(W)))
    return segs


def _fit(r, w, arm_id, beta, cfg):
    segs = segment_profile(r, w, arm_id, beta, cfg.min_seg_bins)
    fit = np.zeros_like(r)
    for a, b, m, se in segs: fit[a:b] = m
    return fit, segs


def _mode_level(segs, tol=0.08):
    """length-weighted modal segment level -> neutral baseline (handles aneuploid genomes)."""
    m = np.array([s[2] for s in segs]); L = np.array([s[1] - s[0] for s in segs], float)
    best, bl = -1, 0.0
    for v in m:
        sel = np.abs(m - v) <= tol
        if L[sel].sum() > best: best = L[sel].sum(); bl = np.average(m[sel], weights=L[sel])
    return bl


def call_segments(segs, fit_shift, cfg, lam=1.0, tau2=0.0):
    """z uses the calibrated segment variance  lam * se^2 + tau2  (white + arm-level correlated component)."""
    out = []
    for a, b, m, se in segs:
        mm = m - fit_shift
        z = mm / np.sqrt(lam * se ** 2 + tau2)
        call = 0
        if abs(mm) >= cfg.min_effect and abs(z) >= cfg.min_z and (b - a) >= cfg.min_seg_bins:
            call = int(np.sign(mm)) * (2 if abs(mm) >= cfg.high_effect else 1)
        out.append((a, b, mm, se, z, call))
    return out


@dataclass
class CNVResult:
    cfg: Config
    bins: pd.DataFrame
    genes: pd.DataFrame
    meta_table: pd.DataFrame
    cell_meta: np.ndarray
    Y: np.ndarray
    r: np.ndarray            # centred log2 ratio (metacells x bins)
    sd: np.ndarray           # null SD per metacell x bin
    z: np.ndarray
    seg: np.ndarray          # fitted piecewise-constant profile (log2)
    call: np.ndarray         # -2..2 per bin
    segments: pd.DataFrame
    group_seg: dict = field(default_factory=dict)
    group_r: dict = field(default_factory=dict)
    group_call: dict = field(default_factory=dict)
    group_segments: pd.DataFrame | None = None
    penalty: float | None = None
    calib: pd.DataFrame | None = None
    r_by_ref: dict = field(default_factory=dict)

    def pool(self, groups):
        """precision-weighted pooled profile over the metacells of the given groups, re-segmented.
        returns (r, seg, call) on the bin axis."""
        ms = np.where(self.meta_table.group.isin(groups).values)[0]
        W = 1.0 / self.sd[ms] ** 2
        rg = (self.r[ms] * W).sum(0) / W.sum(0); sdg = 1 / np.sqrt(W.sum(0))
        c, fit, sg = _center_and_segment(rg, sdg, self.arm_id, self.penalty, self.cfg, lam=self.lam, tau2=self.tau2 / max(len(ms), 1))
        call = np.zeros(len(rg), int)
        for a, b, mm, se, z, cl in sg: call[a:b] = cl
        return rg - c, fit, call

    def cell_matrix(self, which='seg'):
        M = getattr(self, which); out = np.full((len(self.cell_meta), M.shape[1]), np.nan)
        ok = self.cell_meta >= 0; out[ok] = M[self.cell_meta[ok]]; return out


def calibrate_segment_variance(R, SD, arm_id, Ls=(2, 3, 5, 8, 12, 20, 30)):
    """
    Fit  Var(window mean) = lam * Var_white + tau2  on reference metacells: window means over L consecutive bins
    inside an arm vs. the white-noise variance implied by SD. Robust (MAD) so rare real events do not inflate it.
    """
    from scipy.optimize import nnls
    obs, wh = [], []
    W = 1.0 / SD ** 2
    for L in Ls:
        means, whites = [], []
        for a in pd.unique(arm_id):
            ix = np.where(arm_id == a)[0]
            for s in range(0, len(ix) - L + 1, L):
                j = ix[s:s + L]
                ww = W[:, j]; means.append((R[:, j] * ww).sum(1) / ww.sum(1)); whites.append(1.0 / ww.sum(1))
        if not means: continue
        mm = np.concatenate(means); wv = np.concatenate(whites)
        obs.append((1.4826 * np.median(np.abs(mm - np.median(mm)))) ** 2); wh.append(np.median(wv))
    A = np.column_stack([wh, np.ones(len(wh))]); coef, _ = nnls(A / np.array(obs)[:, None], np.ones(len(obs)))
    return max(float(coef[0]), 1.0), float(coef[1])


def _ref_denoise(Rz, ref_idx, k):
    """remove top-k PCs of reference variation (cross-fitted for reference rows). Rz = R/SD."""
    out = Rz.copy()
    ref_idx = np.asarray(ref_idx)
    def proj(train_rows, rows):
        T = Rz[train_rows]; mu = T.mean(0)
        _, _, Vt = np.linalg.svd(T - mu, full_matrices=False); V = Vt[:k]
        d = Rz[rows] - mu
        return Rz[rows] - (d @ V.T) @ V
    nonref = np.setdiff1d(np.arange(len(Rz)), ref_idx)
    if len(nonref): out[nonref] = proj(ref_idx, nonref)
    for m in ref_idx: out[m] = proj(ref_idx[ref_idx != m], [m])[0]
    return out


def _center_and_segment(r, sd, arm_id, beta, cfg, c0=None, lam=1.0, tau2=0.0):
    """centre one metacell profile on its modal copy-number level and segment it."""
    w = 1.0 / sd ** 2
    c = np.median(r) if c0 is None else c0
    for _ in range(cfg.recenter_iters + 1):
        fit, segs = _fit(r - c, w, arm_id, beta, cfg)
        shift = _mode_level(segs)
        c = c + shift
    fit, segs = _fit(r - c, w, arm_id, beta, cfg)
    shift = _mode_level(segs)
    c = c + shift
    segs = [(a, b, m - shift, se) for a, b, m, se in segs]
    fit = fit - shift
    return c, fit, call_segments(segs, 0.0, cfg, lam, tau2)


def _arm_median_resid(r, arm_id):
    out = np.empty_like(r)
    for a in pd.unique(arm_id):
        ix = np.where(arm_id == a)[0]; out[ix] = r[ix] - np.median(r[ix])
    return out


def run(X, gene_names, table, cen, groups, ref_masks: dict, cfg: Config = None, lineage_offset=None,
        offset_var=None, verbose=True, sel=None):
    """
    X: cells x genes raw counts (CSR). groups: cluster label per cell (metacells never cross groups).
    ref_masks: {name: bool mask over cells} reference lineages (e.g. {'stroma':..., 'immune':...}); reference cells
               must sit in groups that contain only that reference set.
    lineage_offset / offset_var: optional per-bin (log2) offset between the epithelial test lineage and the reference
               lineage, and its variance (see lineage_offset_from_normals).
    """
    cfg = cfg or Config(); rng = np.random.default_rng(cfg.seed)
    groups = np.asarray(groups).astype(str)
    ref_union = np.any([np.asarray(m) for m in ref_masks.values()], 0)
    sel = sel if sel is not None else select_genes_and_bins(X, gene_names, table, cen, cfg, groups, ref_union)
    gdf, bins, M, w, a_b = sel['genes'], sel['bins'], sel['M'], sel['w'], sel['a']
    if verbose: print(f"[genes] {len(gdf)} kept -> {len(bins)} bins ({cfg.genes_per_bin} genes/bin, weight_power={cfg.weight_power})")
    Xs = X[:, gdf.col.values].tocsr().astype(np.float32)
    Wd = sp.diags(w.astype(np.float32))
    A, Bf = split_counts(Xs, cfg.split_frac, rng) if cfg.thin else (Xs, Xs)
    A_w, B_w = A @ Wd, Bf @ Wd
    cfg_m = cfg
    for attempt in range(4):          # shrink metacells if the reference ends up with too few of them
        meta, mgroup = build_metacells(A_w, B_w, M, groups, cfg_m, rng)
        ok = meta >= 0
        n_ref = min(len(np.unique(meta[np.asarray(mk) & ok])) for mk in ref_masks.values())      # smallest reference set
        small = (cfg_m.target_umi is None and cfg_m.metacell_cells <= cfg.min_metacell_cells) or cfg_m.metacell_cells <= 1
        if n_ref >= cfg.min_ref_metacells or small: break
        cfg_m = replace(cfg_m, target_umi=cfg_m.target_umi / 2) if cfg_m.target_umi else replace(cfg_m, metacell_cells=max(cfg.min_metacell_cells, cfg_m.metacell_cells // 2))
        if verbose: print(f"[metacells] only {n_ref} reference metacells -> shrinking metacell size ({'target_umi=%g' % cfg_m.target_umi if cfg_m.target_umi else 'cells=%d' % cfg_m.metacell_cells})")
    Y, ncell = aggregate(B_w, meta, M)
    Ncnt = np.bincount(meta[ok], weights=np.asarray(Bf.sum(1)).ravel()[ok], minlength=len(mgroup))
    if verbose: print(f"[metacells] {len(mgroup)} metacells, median {np.median(ncell):.0f} cells, {int((~ok).sum())} cells dropped")
    ref_sets, ref_names = [], []
    for name, mask in ref_masks.items():
        mask = np.asarray(mask); rows = np.unique(meta[mask & ok])
        for m in rows:
            assert mask[meta == m].mean() > 0.99, f"metacell {m} mixes reference '{name}' with other cells; refine groups"
        ref_sets.append(list(rows)); ref_names.append(name)
    all_ref = sorted({r for rs in ref_sets for r in rs}); is_ref = np.zeros(len(mgroup), bool); is_ref[all_ref] = True
    if verbose: print("[reference] " + ", ".join(f"{n}={len(r)} metacells" for n, r in zip(ref_names, ref_sets)))
    assert len(all_ref) >= 6, "need >= 6 reference metacells to estimate the null"
    p_full = baseline_logp(Y, ref_sets, cfg)
    R = np.zeros(Y.shape)
    for m in range(len(Y)):
        p = baseline_logp(Y, ref_sets, cfg, loo_row=m) if is_ref[m] else p_full
        R[m] = log_ratio(Y[m:m + 1], p, cfg)[0][0]
    # lineage offset: applied to QUERY metacells only (never to reference metacells, so the reference null / calibration stay offset-free).
    #   array -> every non-reference metacell; dict {group label: array} -> only those groups (others get 0)
    def _spec(spec, g):
        if spec is None: return np.zeros(len(bins))
        if isinstance(spec, dict): return np.asarray(spec[g], float) if g in spec else np.zeros(len(bins))
        return np.asarray(spec, float)
    off_mat = np.vstack([np.zeros(len(bins)) if is_ref[m] else _spec(lineage_offset, mgroup[m]) for m in range(len(Y))])
    ov_mat = np.vstack([np.zeros(len(bins)) if is_ref[m] else _spec(offset_var, mgroup[m]) for m in range(len(Y))])
    R_adj = R - off_mat
    pois = a_b[None] / (LN2 ** 2 * np.maximum(Ncnt[:, None], 1.0))

    def null_sd(Rm):
        """per-bin extra variance from reference LOO residuals (median-centred) + Poisson + offset variance."""
        Rr = Rm[all_ref] - np.median(Rm[all_ref], 1, keepdims=True)
        mad = 1.4826 * np.median(np.abs(Rr - np.median(Rr, 0)), 0)
        e_raw = np.maximum(mad ** 2 - pois[all_ref].mean(0), 0)
        e = np.maximum(cfg.var_shrink * np.median(e_raw) + (1 - cfg.var_shrink) * e_raw, 1e-4)
        return np.sqrt(e[None] + ov_mat + pois)

    SD = null_sd(R_adj)
    if cfg.n_ref_pcs > 0:
        Rm = R_adj - np.median(R_adj, 1, keepdims=True)
        R_adj = _ref_denoise(Rm / SD, all_ref, cfg.n_ref_pcs) * SD
        SD = null_sd(R_adj)
    arm_id = (bins.chr + bins.arm).values; n_arm = len(pd.unique(arm_id))
    lam, tau2 = calibrate_segment_variance(R_adj[all_ref] - np.median(R_adj[all_ref], 1, keepdims=True), SD[all_ref], arm_id)
    if verbose: print(f"[null] lambda={lam:.2f}, tau={np.sqrt(tau2):.3f}")
    # ---- calibration on reference metacells (their own LOO profiles)
    cfg0 = replace(cfg, min_z=0.0); calib = None; beta = cfg.seg_penalty
    if beta is None:
        rows = []
        for b in cfg.penalty_grid:
            nseg = np.mean([len(_center_and_segment(R_adj[m], SD[m], arm_id, b, cfg0, lam=lam, tau2=tau2)[2]) for m in all_ref])
            rows.append((b, nseg, nseg / n_arm))
        calib = pd.DataFrame(rows, columns=['penalty', 'ref_mean_segments', 'ratio_to_arms'])
        okc = calib[calib.ratio_to_arms <= 1.05]
        beta = float(okc.penalty.iloc[0]) if len(okc) else float(calib.penalty.iloc[-1])
    zmax = []
    for m in all_ref:
        sg = _center_and_segment(R_adj[m], SD[m], arm_id, beta, cfg0, lam=lam, tau2=tau2)[2]
        zz = [abs(s[4]) for s in sg if abs(s[2]) >= cfg.min_effect and (s[1] - s[0]) >= cfg.min_seg_bins]
        zmax.append(max(zz) if zz else 0.0)
    zmax = np.array(zmax); z_call = cfg.min_z
    for zc in np.arange(cfg.min_z, max(zmax.max(), cfg.min_z) + 0.5, 0.25):
        z_call = float(zc)
        if (zmax >= zc).mean() <= cfg.fp_target: break
    cfg = replace(cfg, min_z=z_call)
    if verbose: print(f"[calibration] penalty {beta}; call |z| >= {z_call:.2f} (reference max |z| median {np.median(zmax):.1f}, max {zmax.max():.1f})")

    def process(SDm):
        Rc = np.zeros_like(R_adj); FIT = np.zeros_like(R_adj); CALL = np.zeros(R_adj.shape, int); rows = []
        for m in range(len(Y)):
            c, fit, sg = _center_and_segment(R_adj[m], SDm[m], arm_id, beta, cfg, lam=lam, tau2=tau2)
            Rc[m] = R_adj[m] - c; FIT[m] = fit
            for a, b, mm, se, z, call in sg:
                CALL[m, a:b] = call
                rows.append((m, mgroup[m], bins.chr[a], bins.start[a], bins.end[b - 1], a, b, b - a, mm, se, z, call))
        return Rc, FIT, CALL, rows

    Rc, FIT, CALL, seg_rows = process(SD)
    # ---- group-specific floor: bin-level expression differences between a test group and the reference that are not
    # copy number (cell-type / state mismatch). Estimated robustly from the group profile around arm medians.
    floor = {}
    Wt = 1.0 / SD ** 2
    for g in pd.unique(mgroup):
        ms = np.where(mgroup == g)[0]
        if not cfg.floor or is_ref[ms].all(): floor[g] = 0.0; continue
        rg = (Rc[ms] * Wt[ms]).sum(0) / Wt[ms].sum(0); white = np.median(1.0 / Wt[ms].sum(0))
        res_ = _arm_median_resid(rg, arm_id)
        floor[g] = float(np.sqrt(max((1.4826 * np.median(np.abs(res_ - np.median(res_)))) ** 2 - white, 0.0)))
    if cfg.floor:
        SD2 = np.sqrt(SD ** 2 + np.array([floor[g] for g in mgroup])[:, None] ** 2)
        Rc, FIT, CALL, seg_rows = process(SD2)
    else:
        SD2 = SD
    segs = pd.DataFrame(seg_rows, columns=['metacell', 'group', 'chr', 'start', 'end', 'bin_start', 'bin_end', 'n_bins',
                                           'log2_ratio', 'se', 'z', 'call'])
    gseg, gr, gcall, grow = {}, {}, {}, []
    W = 1.0 / SD2 ** 2
    for g in pd.unique(mgroup):
        ms = np.where(mgroup == g)[0]
        rg = (Rc[ms] * W[ms]).sum(0) / W[ms].sum(0); sdg = 1 / np.sqrt(W[ms].sum(0))
        c, fit, sg = _center_and_segment(rg, sdg, arm_id, beta, cfg, lam=lam, tau2=tau2 / max(len(ms), 1))
        gr[g] = rg - c; gseg[g] = fit; cl = np.zeros(len(bins), int)
        for a, b, mm, se, z, call in sg:
            cl[a:b] = call; grow.append((g, bins.chr[a], bins.start[a], bins.end[b - 1], a, b, b - a, mm, se, z, call))
        gcall[g] = cl
    gsegs = pd.DataFrame(grow, columns=['group', 'chr', 'start', 'end', 'bin_start', 'bin_end', 'n_bins', 'log2_ratio', 'se', 'z', 'call'])
    mt = pd.DataFrame(dict(metacell=np.arange(len(mgroup)), group=mgroup, n_cells=ncell, umi=Ncnt, is_ref=is_ref,
                           floor=[floor[g] for g in mgroup]))
    out = CNVResult(cfg, bins, gdf, mt, meta, Y, Rc, SD2, Rc / SD2, FIT, CALL, segs, gseg, gr, gcall, gsegs, beta, calib)
    out.lam, out.tau2, out.z_call, out.floor, out.arm_id = lam, tau2, z_call, floor, arm_id
    return out



# --------------------------------------------------------------------------- count utilities
def stochastic_round(X, rng):
    """unbiased integer rounding of non-integer (e.g. ambient-corrected) counts: floor(x) + Bernoulli(frac)."""
    X = X.tocsr().copy(); x = X.data.astype(np.float64)
    fl = np.floor(x); X.data = (fl + (rng.random(len(x)) < (x - fl))).astype(np.float32)
    X.eliminate_zeros(); return X


def downsample_counts(X, frac, rng):
    """binomial thinning of every cell's UMIs to a fraction `frac` (simulates lower sequencing depth)."""
    if frac >= 1: return X
    X = X.tocsr().copy(); X.data = rng.binomial(X.data.astype(np.int64), frac).astype(np.float32)
    X.eliminate_zeros(); return X


def infercnv_like(X, sel, ref_mask, window=101, clip=3.0, chunk=4000):
    """
    Fast single-cell inferCNV-style baseline on the same genes/bins as the pipeline:
    log-CPM -> subtract reference mean per gene -> clip -> running mean over `window` genes within each chromosome ->
    per-cell median centring -> mean per bin. Returns cells x bins.
    """
    gdf, bins, M = sel['genes'], sel['bins'], sel['M']
    Xs = X[:, gdf.col.values].tocsr().astype(np.float32)
    tot = np.asarray(X.sum(1)).ravel().astype(np.float32)[:, None]
    refmean = None
    out = np.zeros((X.shape[0], len(bins)), np.float32)
    # reference means (log-CPM) on the full gene set
    ref_idx = np.where(ref_mask)[0]
    L = np.log1p(Xs[ref_idx].multiply(1e4 / np.maximum(tot[ref_idx], 1)).toarray()); refmean = L.mean(0); del L
    chr_i = gdf.chr_i.values; starts = np.r_[0, np.where(np.diff(chr_i) != 0)[0] + 1, len(gdf)]
    Mb = sp.csr_matrix(M); cnt = np.asarray(Mb.sum(0)).ravel()
    for s in range(0, X.shape[0], chunk):
        sl = slice(s, min(s + chunk, X.shape[0]))
        D = np.log1p(Xs[sl].multiply(1e4 / np.maximum(tot[sl], 1)).toarray()) - refmean
        D = np.clip(D, -clip, clip)
        S = np.empty_like(D)
        for a, b in zip(starts[:-1], starts[1:]):
            c = np.cumsum(np.pad(D[:, a:b], ((0, 0), (1, 0))), axis=1); n = b - a; h = window // 2
            idx = np.arange(n); lo = np.clip(idx - h, 0, n); hi = np.clip(idx + h + 1, 0, n)
            S[:, a:b] = (c[:, hi] - c[:, lo]) / (hi - lo)
        S -= np.median(S, axis=1, keepdims=True)
        out[sl] = (S @ Mb.toarray()) / np.maximum(cnt, 1)
    return out


# --------------------------------------------------------------------------- lineage utilities
def lineage_offset_from_normals(Y_epi, Y_ref, donor_epi, donor_ref, pseudocount=0.5):
    """
    Per-bin offset (log2) between a normal epithelial control and a reference lineage, estimated *within* each
    donor/dataset (so dataset effects cancel) and combined across donors.
    Y_*: (units x bins) aggregated bin counts; donor_*: label per unit.
    Returns offset (median over donors) and variance (between-donor, robust) for use as run(lineage_offset=, offset_var=).
    """
    ds = sorted(set(donor_epi) & set(donor_ref)); assert len(ds) >= 1
    D = []
    for d in ds:
        a = Y_epi[np.asarray(donor_epi) == d].sum(0) + pseudocount; b = Y_ref[np.asarray(donor_ref) == d].sum(0) + pseudocount
        D.append(np.log2(a / a.sum()) - np.log2(b / b.sum()))
    D = np.vstack(D)
    off = np.median(D, 0)
    var = (1.4826 * np.median(np.abs(D - off), 0)) ** 2 if len(ds) >= 3 else np.full(D.shape[1], np.var(D) if len(ds) > 1 else 0.0)
    return off, var, D


def bin_truth(res_bins, genes_df, truth: pd.Series):
    """aggregate gene-level truth (log2 scale) to the pipeline's bins (median over genes)."""
    t = truth.reindex(genes_df.symbol.values).values
    s = pd.Series(t).groupby(genes_df.bin.values).median()
    return s.reindex(np.arange(len(res_bins))).values


# --------------------------------------------------------------------------- step 1: clustering (optional)
def cluster_cells(X, n_hvg=2000, n_pcs=30, resolution=0.5, seed=0):
    """Leiden clustering on log-normalised HVGs (scanpy). Use when no curated clusters / annotation exist.
    Metacells are later built inside each cluster, so over-clustering is safer than under-clustering."""
    import scanpy as sc, anndata as ad
    a = ad.AnnData(sp.csr_matrix(X, dtype=np.float32))
    sc.pp.normalize_total(a, target_sum=1e4); sc.pp.log1p(a)
    sc.pp.highly_variable_genes(a, n_top_genes=n_hvg, flavor='seurat'); a = a[:, a.var.highly_variable].copy()
    sc.pp.scale(a, max_value=10); sc.tl.pca(a, n_comps=n_pcs, random_state=seed)
    sc.pp.neighbors(a, n_pcs=n_pcs, random_state=seed); sc.tl.leiden(a, resolution=resolution, random_state=seed, flavor='igraph', n_iterations=2, directed=False)
    return a.obs['leiden'].astype(str).values
