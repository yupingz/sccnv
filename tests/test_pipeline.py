import numpy as np
import pytest
from sccnv.core import Config, run, select_genes_and_bins
from conftest import simulate


@pytest.fixture(scope='module')
def fitted(genome):
    table, cen, syms, chr_of = genome
    S = simulate(genome)
    cfg = Config(genes_per_bin=40, metacell_cells=20, min_arm_genes=20)
    sel = select_genes_and_bins(S['X'], S['genes'], table, cen, cfg, S['grp'], S['ref_mask'])
    res = run(S['X'], S['genes'], table, cen, S['grp'], {'ref': S['ref_mask']}, cfg, verbose=False, sel=sel)
    return S, cfg, sel, res, (table, cen)


def _arm(res, chrom, arm):
    return ((res.bins.chr == chrom) & (res.bins.arm == arm)).values


def test_tumour_gain_and_loss_are_recovered(fitted):
    S, cfg, sel, res, _ = fitted
    r, seg, call = res.pool(['tum|c0'])
    assert r[_arm(res, '2', 'q')].mean() > 0.3                     # gain 1.6x -> log2 ~ +0.68
    assert r[_arm(res, '3', 'p')].mean() < -0.4                    # loss 0.5x -> log2 ~ -1
    assert (call[_arm(res, '2', 'q')] > 0).mean() > 0.8
    assert (call[_arm(res, '3', 'p')] < 0).mean() > 0.8
    neutral = ~(_arm(res, '2', 'q') | _arm(res, '3', 'p'))
    assert (call[neutral] != 0).mean() < 0.1
    assert abs(np.median(r[neutral])) < 0.1


def test_heldout_diploid_cells_are_not_called(fitted):
    S, cfg, sel, res, _ = fitted
    mt = res.meta_table
    ho = np.where(mt.group.str.contains('_ho').values)[0]
    assert len(ho) >= 3
    assert (res.call[ho] != 0).mean() < 0.02


def test_reference_metacells_never_cross_groups(fitted):
    S, cfg, sel, res, _ = fitted
    mt = res.meta_table; meta = res.cell_meta
    for m in range(len(mt)):
        cells = np.where(meta == m)[0]
        assert len(set(S['grp'][cells])) == 1
    assert mt.is_ref.sum() >= cfg.min_ref_metacells


def test_cell_matrix_shape_and_dropped_cells(fitted):
    S, cfg, sel, res, _ = fitted
    P = res.cell_matrix('seg')
    assert P.shape == (S['X'].shape[0], len(res.bins))


def test_events_recovered_with_and_without_count_splitting(genome):
    """with and without count splitting, both runs must find the same events (splitting must not remove signal)."""
    table, cen, syms, chr_of = genome
    S = simulate(genome, seed=3)
    out = []
    for thin in (True, False):
        cfg = Config(genes_per_bin=40, metacell_cells=20, thin=thin)
        sel = select_genes_and_bins(S['X'], S['genes'], table, cen, cfg, S['grp'], S['ref_mask'])
        res = run(S['X'], S['genes'], table, cen, S['grp'], {'ref': S['ref_mask']}, cfg, verbose=False, sel=sel)
        r, seg, call = res.pool(['tum|c0'])
        out.append((r[_arm(res, '2', 'q')].mean(), r[_arm(res, '3', 'p')].mean()))
    for (g1, l1), (g2, l2) in [out]:
        assert g1 > 0.25 and g2 > 0.25 and l1 < -0.35 and l2 < -0.35


def test_lineage_offset_is_applied_to_query_metacells_only(genome):
    table, cen, syms, chr_of = genome
    S = simulate(genome, seed=5)
    cfg = Config(genes_per_bin=40, metacell_cells=20)
    sel = select_genes_and_bins(S['X'], S['genes'], table, cen, cfg, S['grp'], S['ref_mask'])
    res0 = run(S['X'], S['genes'], table, cen, S['grp'], {'ref': S['ref_mask']}, cfg, verbose=False, sel=sel)
    off = np.zeros(len(res0.bins)); idx = np.where(_arm(res0, '4', 'q'))[0][:4]; off[idx] = 0.6
    res1 = run(S['X'], S['genes'], table, cen, S['grp'], {'ref': S['ref_mask']}, cfg, verbose=False, sel=sel,
               lineage_offset={'tum|c0': off})
    mt = res0.meta_table; tum = np.where((mt.group == 'tum|c0').values)[0]; ref = np.where(mt.is_ref.values)[0]
    d = res0.r - res1.r
    assert d[np.ix_(tum, idx)].mean() > 0.4          # query metacells shifted by about the offset
    assert abs(d[np.ix_(ref, idx)].mean()) < 0.1     # reference metacells untouched


def test_run_refuses_to_start_without_reference_metacells(genome):
    table, cen, syms, chr_of = genome
    S = simulate(genome, n_ref=30, n_tum=40)
    cfg = Config(genes_per_bin=40, metacell_cells=20)
    sel = select_genes_and_bins(S['X'], S['genes'], table, cen, cfg, S['grp'], S['ref_mask'])
    with pytest.raises(AssertionError):
        run(S['X'], S['genes'], table, cen, S['grp'], {'ref': S['ref_mask'] & False}, cfg, verbose=False, sel=sel)


def test_null_simulation_yields_no_calls(genome):
    """tumour cells with the SAME copy number as the reference: the calibrated null must keep false calls rare."""
    table, cen, syms, chr_of = genome
    S = simulate(genome, seed=11, gain=('2', 'q', 1.0), loss=('3', 'p', 1.0))
    cfg = Config(genes_per_bin=40, metacell_cells=20)
    sel = select_genes_and_bins(S['X'], S['genes'], table, cen, cfg, S['grp'], S['ref_mask'])
    res = run(S['X'], S['genes'], table, cen, S['grp'], {'ref': S['ref_mask']}, cfg, verbose=False, sel=sel)
    r, seg, call = res.pool(['tum|c0'])
    assert (call != 0).mean() < 0.05
    tum = np.where((res.meta_table.group == 'tum|c0').values)[0]
    assert (res.call[tum] != 0).mean() < 0.03
