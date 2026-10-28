import pytest
from sccnv.simulate import toy_genome, simulate as _simulate


@pytest.fixture(scope='session')
def genome():
    return toy_genome()


def simulate(genome, **kw):
    return _simulate(genome, **kw)
