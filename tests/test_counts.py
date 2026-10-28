import numpy as np, scipy.sparse as sp
from sccnv.core import stochastic_round, downsample_counts, split_counts


def test_stochastic_round_is_unbiased_and_integer():
    rng = np.random.default_rng(0)
    X = sp.random(2000, 300, density=0.2, random_state=1, format='csr') * 3.0
    Xi = stochastic_round(X, rng)
    assert np.all(Xi.data == np.round(Xi.data))
    assert abs(Xi.sum() / X.sum() - 1) < 0.01


def test_downsample_counts_fraction():
    rng = np.random.default_rng(0)
    X = sp.csr_matrix(rng.poisson(5, (500, 200)).astype(np.float32))
    Xd = downsample_counts(X, 0.1, rng)
    assert abs(Xd.sum() / X.sum() - 0.1) < 0.01
    assert (Xd.toarray() <= X.toarray()).all()
    assert downsample_counts(X, 1.0, rng) is X


def test_split_counts_partitions_every_umi():
    rng = np.random.default_rng(0)
    X = sp.csr_matrix(rng.poisson(4, (300, 150)).astype(np.float32))
    A, B = split_counts(X, 0.5, rng)
    assert (A + B - X).nnz == 0 or np.abs((A + B - X).toarray()).max() == 0
    assert 0.45 < A.sum() / X.sum() < 0.55
