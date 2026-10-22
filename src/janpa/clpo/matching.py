# janpa/clpo/matching.py
"""Maximum-weight matching for bonding graph construction.

Replaces the 600-line Blossom.java. Uses the compiled (Rust) blossom algorithm of rustworkx; networkx implements the
same algorithm in pure Python and is only used as a fallback, as it is about two orders of magnitude slower on the
complete hybrid graphs of the LPO step.
"""

from __future__ import annotations

try:
    import rustworkx as rx
except ImportError:  # pragma: no cover
    rx = None

# rustworkx needs integer weights: the largest weight is mapped to 2**50, which keeps ~15 significant digits (as many
# as a double) while leaving room in 64 bits for the sums of weights formed by the algorithm.
_WEIGHT_SCALE_BITS = 50


def max_weight_matching_pairs(n_nodes: int, edges: list[tuple[int, int, float]]) -> list[tuple[int, int]]:
    """Find maximum-weight (not maximum-cardinality) matching on a weighted graph.

    Parameters
    ----------
    n_nodes : number of vertices, labelled 0 .. n_nodes - 1
    edges : list of (u, v, weight) tuples

    Returns
    -------
    list of matched (u, v) pairs
    """
    edges = [(u, v, w) for u, v, w in edges if u != v]
    if not edges:
        return []
    if rx is None:
        import networkx as nx
        G = nx.Graph()
        G.add_nodes_from(range(n_nodes))
        G.add_weighted_edges_from(edges)
        return list(nx.max_weight_matching(G, maxcardinality=False, weight="weight"))

    max_w = max(abs(w) for _, _, w in edges)
    scale = 2.0 ** _WEIGHT_SCALE_BITS / max_w if max_w > 0 else 1.0
    G = rx.PyGraph(multigraph=False)
    G.add_nodes_from(range(n_nodes))
    G.add_edges_from([(u, v, round(w * scale)) for u, v, w in edges])
    return list(rx.max_weight_matching(G, max_cardinality=False, weight_fn=int))


def max_weight_matching(edges: list[tuple[int, int, float]]) -> dict[int, int]:
    """Find maximum-weight matching on a weighted graph.

    Parameters
    ----------
    edges : list of (u, v, weight) tuples

    Returns
    -------
    matching : dict mapping vertex → matched vertex (-1 if unmatched)
    """
    n_nodes = max((max(u, v) for u, v, _ in edges), default=-1) + 1
    matching: dict[int, int] = {}
    for u, v in max_weight_matching_pairs(n_nodes, edges):
        matching[u] = v
        matching[v] = u
    return matching
