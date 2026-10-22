import numpy as np
import networkx as nx
import janpa.clpo.matching as matching
from janpa.utils.simult_diag import simultaneous_diagonalization


def _random_graph(rng, n, p):
    return [(u, v, float(rng.uniform(0, 1) ** rng.integers(1, 6))) for u in range(n) for v in range(u) if rng.random() < p]


def _weight(edges, pairs):
    w = {frozenset((u, v)): x for u, v, x in edges}
    return sum(w[frozenset(e)] for e in pairs)


def test_max_weight_matching_like_networkx(monkeypatch):
    """The rustworkx matching (integer-scaled weights) and the networkx fallback find matchings as heavy as networkx's."""
    rng = np.random.default_rng(0)
    for n, p in [(6, 1.0), (30, 0.3), (80, 1.0)]:
        edges = _random_graph(rng, n, p)
        g = nx.Graph()
        g.add_weighted_edges_from(edges)
        expected = _weight(edges, nx.max_weight_matching(g))
        assert np.isclose(_weight(edges, matching.max_weight_matching_pairs(n, edges)), expected, rtol=1e-12)
        with monkeypatch.context() as m:
            m.setattr(matching, "rx", None)
            assert np.isclose(_weight(edges, matching.max_weight_matching_pairs(n, edges)), expected, rtol=1e-12)
    assert matching.max_weight_matching([(0, 1, 1.0), (1, 2, 3.0), (2, 3, 1.0)]) == {1: 2, 2: 1}


def test_simultaneous_diagonalization():
    """Matrices sharing a set of eigenvectors are diagonalized together, and the inputs are rotated in place."""
    rng = np.random.default_rng(0)
    q, _ = np.linalg.qr(rng.standard_normal((6, 6)))
    matrices = [q @ np.diag(rng.standard_normal(6)) @ q.T for _ in range(3)]
    originals = [m.copy() for m in matrices]
    u, n_iters = simultaneous_diagonalization(matrices, max_iters=100, threshold=1e-12)
    assert n_iters < 100
    assert np.allclose(u.T @ u, np.eye(6), atol=1e-12)
    for a, rotated in zip(originals, matrices):
        d = u.T @ a @ u
        assert np.abs(d - np.diag(np.diag(d))).max() < 1e-10
        assert np.allclose(rotated, d, atol=1e-10)
