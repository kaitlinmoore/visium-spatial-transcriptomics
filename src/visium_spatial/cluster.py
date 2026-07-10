"""Leiden clustering as a compartment proxy.

Neighbors in *expression* space (``sc.pp.neighbors``) then Leiden
(``sc.tl.leiden``). These clusters stand in for tissue compartments and are the
grouping variable for the secondary neighborhood-enrichment read-out.

Important distinction the project keeps honest: this graph is the
**expression** kNN graph, NOT the spatial weights graph. The spatial graph
(build_graph.py) is a separate object built from coordinates. Conflating the two
is a classic spatial-omics error — clusters are transcriptomic, enrichment asks
whether those transcriptomic clusters are spatially adjacent.

Alternative considered: mclust / Louvain. Leiden is the current scanpy default
and gives well-defined, reproducible communities; resolution is exposed so the
compartment granularity is a visible knob.
"""

from __future__ import annotations


def leiden_clusters(adata, *, resolution: float = 1.0, n_neighbors: int = 15, key_added: str = "leiden"):
    """Run expression-space neighbors + Leiden; store labels in ``obs[key_added]``.

    Note: expression graph, not spatial. See module docstring.
    """
    raise NotImplementedError
