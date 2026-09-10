# janpa/clpo/matching.py
"""Maximum-weight matching for bonding graph construction.

Replaces the 600-line Blossom.java with a single networkx call.
"""

from __future__ import annotations
import networkx as nx


def max_weight_matching(edges: list[tuple[int, int, float]]) -> dict[int, int]:
    """Find maximum-weight matching on a weighted graph.

    Parameters
    ----------
    edges : list of (u, v, weight) tuples

    Returns
    -------
    matching : dict mapping vertex → matched vertex (-1 if unmatched)
    """
    G = nx.Graph()
    for u, v, w in edges:
        if u != v:
            G.add_edge(u, v, weight=w)
    result = nx.max_weight_matching(G, maxcardinality=False, weight="weight")
    matching: dict[int, int] = {}
    for u, v in result:
        matching[u] = v
        matching[v] = u
    return matching