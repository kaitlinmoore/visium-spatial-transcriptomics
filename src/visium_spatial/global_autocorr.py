"""Rank spatially variable genes by global Moran's I (squidpy).

The cheap gate before the expensive per-spot pass: score every (QC-passing) gene
by global Moran's I with ``sq.gr.spatial_autocorr(mode="moran")`` on the shared
graph, rank, and hand only the top set to the local indicators (local_autocorr).

Why global first: global Moran's I is one number per gene and answers "is there
*any* spatial structure here?". LISA/Gi* are per-spot and far costlier; running
them genome-wide wastes compute on genes with no signal and worsens the
multiple-testing burden. Rank, then localize the hits.

The from-scratch implementation (moran_scratch.py) must agree with this squidpy
output on the same weights before the ranking is trusted — that agreement is the
defensibility check, run once, not per gene.

Alternatives considered:
- Sepal / other SVG methods (squidpy also ships ``mode="geary"``): different
  notions of spatial variability; Moran's I is the one that extends cleanly to
  the local indicators, so the global and local layers stay the same statistic.
"""

from __future__ import annotations


def rank_svgs(adata, *, genes=None, n_perms: int | None = None, seed: int = 0):
    """Score genes by global Moran's I and return them ranked.

    Wraps ``sq.gr.spatial_autocorr(adata, mode="moran", ...)`` on the graph built
    by build_graph.py. ``genes=None`` scores all; pass a list to restrict.

    Returns
    -------
    pandas.DataFrame
        Per-gene Moran's I, p-values, and rank, sorted by I descending.
    """
    raise NotImplementedError


def top_genes(ranked, *, n: int = 20, max_pval: float | None = 0.05):
    """Select the top-``n`` spatially variable genes to localize.

    Filters on the (corrected) global p-value first when ``max_pval`` is set,
    then takes the ``n`` highest-I genes.
    """
    raise NotImplementedError
