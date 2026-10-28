import numpy as np, pandas as pd
from sccnv import benchmark as bench, offsets as O
from sccnv.core import bin_truth


def test_truth_classes_and_mode_center():
    t = np.array([0.0, 0.02, 0.3, -0.3, 0.15, np.nan])
    c = bench.truth_classes(t)
    assert list(c[:4]) == [0, 0, 1, -1] and np.isnan(c[4]) and np.isnan(c[5])
    v = np.r_[np.full(60, 1.0), np.full(10, 1.6), np.full(10, 0.4)]
    assert abs(np.median(bench.mode_center(v))) < 1e-6


def test_eval_profiles_perfect_and_null_predictions():
    rng = np.random.default_rng(0); nb = 120
    truth = np.zeros(nb); truth[:20] = 0.5; truth[20:40] = -0.5
    n_t, n_n = 60, 60
    P = np.vstack([truth + rng.normal(0, 0.05, nb) for _ in range(n_t)] + [rng.normal(0, 0.05, nb) for _ in range(n_n)])
    key = np.array(['t'] * n_t + [''] * n_n); norm = key == ''
    d = bench.eval_profiles(P, key, {'t': truth}, norm)[0]
    assert d['sens_gain_fpr10'] > 0.9 and d['sens_loss_fpr10'] > 0.9 and d['auc_gain'] > 0.99 and d['pearson'] > 0.9
    Pn = rng.normal(0, 0.05, P.shape)
    d0 = bench.eval_profiles(Pn, key, {'t': truth}, norm)[0]
    assert d0['auc_gain'] < 0.6 and abs(d0['pearson']) < 0.3


def test_truth_threshold_arguments_are_honoured():
    truth = np.r_[np.full(30, 0.15), np.zeros(60), np.full(30, -0.15)]
    P = np.tile(truth, (40, 1)) + np.random.default_rng(0).normal(0, 0.02, (40, 120))
    key = np.array(['t'] * 40); norm = np.zeros(40, bool); norm[:5] = True; key[:5] = ''
    a = bench.eval_profiles(P, key, {'t': truth}, norm)[0]                                   # default +-0.25: no gain/loss bins
    b = bench.eval_profiles(P, key, {'t': truth}, norm, gain=0.12, loss=-0.12, neutral=0.06)[0]
    assert a['n_gain'] == 0 and b['n_gain'] == 30 and b['n_loss'] == 30


def test_segments_to_gene_log2_and_bin_truth():
    table = pd.DataFrame({'chr': ['1'] * 6, 'pos': [10, 20, 30, 40, 50, 60]}, index=[f'g{i}' for i in range(6)])
    seg = pd.DataFrame(dict(chr=['1', '1'], start=[0, 35], end=[34, 100], cn=[4.0, 2.0], status=['Gain', 'Neutral']))
    gl = bench.segments_to_gene_log2(seg, table, list(table.index))
    assert np.allclose(gl.values, [2, 2, 2, 1, 1, 1])
    bins = pd.DataFrame({'bin': [0, 1]}); genes = pd.DataFrame({'symbol': table.index, 'bin': [0, 0, 0, 1, 1, 1]})
    assert np.allclose(bin_truth(bins, genes, gl), [2, 1])


def test_offset_table_recovers_known_shift_and_cancels_donor_effects():
    rng = np.random.default_rng(0); nb = 60; shift = np.zeros(nb); shift[10:25] = 0.8
    Y, ct, unit = [], [], []
    for u in range(3):
        donor_fx = rng.normal(0, 0.3, nb)                              # donor-wide effect shared by epithelium and reference
        base = np.exp(rng.normal(0, 0.3, nb)) * 200
        for kind, fx in (('epi', shift), ('stroma', np.zeros(nb))):
            p = base * 2 ** (donor_fx + fx)
            Y.append(rng.poisson(p * 40, (80, nb)).astype(float) / 40); ct += [kind] * 80; unit += [f'D{u}'] * 80
    Y = np.vstack(Y)
    units, D = O.offset_table(Y, np.array(ct), np.array(unit), ['epi'], [['stroma']])
    off, var = O.combine_units(D, min_units=3)
    assert len(units) == 3
    assert np.corrcoef(off, shift)[0, 1] > 0.9 and off[10:25].mean() - off[:10].mean() > 0.5
    assert var.shape == off.shape and (var >= 0).all()
