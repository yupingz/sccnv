import numpy as np
from sccnv.core import Config, select_genes_and_bins
from conftest import simulate


def test_bins_are_ordered_within_arms_and_sized_as_requested(genome):
    table, cen, syms, chr_of = genome
    S = simulate(genome, n_ref=120, n_tum=60)
    cfg = Config(genes_per_bin=50, min_arm_genes=20)
    sel = select_genes_and_bins(S['X'], S['genes'], table, cen, cfg, S['grp'], S['ref_mask'])
    b, g = sel['bins'], sel['genes']
    assert (b.groupby(['chr', 'arm']).size() >= 1).all() and len(b) >= 8 * 4 // 2
    assert (g.groupby('bin').size().between(30, 70)).all()
    assert g.bin.is_monotonic_increasing
    assert (b.start <= b.end).all()
    for _, grp in g.groupby(['chr', 'arm']):
        assert grp.pos.is_monotonic_increasing and grp.bin.nunique() >= 1
    assert sel['M'].shape == (len(g), len(b)) and np.all(sel['w'] > 0)


def test_blacklist_regex_removes_matching_genes(genome):
    table, cen, syms, chr_of = genome
    S = simulate(genome, n_ref=120, n_tum=60)
    cfg = Config(genes_per_bin=50, blacklist_regex=r'^G1_')
    sel = select_genes_and_bins(S['X'], S['genes'], table, cen, cfg, S['grp'], S['ref_mask'])
    assert not sel['genes'].symbol.str.startswith('G1_').any()
    assert set(sel['bins'].chr) == {'2', '3', '4'}
