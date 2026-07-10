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

import numpy as np
import squidpy as sq

CONNECTIVITIES_KEY = "spatial_connectivities"


def _corrected_pval_col(df) -> str:
    """The multiple-testing-corrected p-value column squidpy wrote.

    squidpy names it ``pval_<norm|sim>_fdr_bh`` depending on ``n_perms``; fall
    back to the raw ``pval_*`` column if no corrected one is present.
    """
    for col in df.columns:
        if col.endswith("_fdr_bh"):
            return col
    for col in df.columns:
        if col.startswith("pval"):
            return col
    raise KeyError(f"no p-value column in {list(df.columns)}")


def rank_svgs(adata, *, genes=None, n_perms: int | None = None, seed: int = 0):
    """Score genes by global Moran's I and return them ranked.

    Wraps ``sq.gr.spatial_autocorr(adata, mode="moran", ...)`` on the graph built
    by build_graph.py. ``genes=None`` scores all; pass a list to restrict.
    ``n_perms=None`` uses squidpy's analytic (normal) p-values; an int uses the
    permutation null.

    Returns
    -------
    pandas.DataFrame
        Per-gene Moran's I, p-values (raw + BH-corrected across the tested genes),
        and an integer ``rank``, sorted by I descending.
    """
    if CONNECTIVITIES_KEY not in adata.obsp:
        raise ValueError(
            f'no spatial graph in obsp["{CONNECTIVITIES_KEY}"]; '
            "run build_graph.build_spatial_graph first"
        )

    genes = list(adata.var_names) if genes is None else list(genes)
    sq.gr.spatial_autocorr(adata, mode="moran", genes=genes, n_perms=n_perms, seed=seed)

    ranked = adata.uns["moranI"].sort_values("I", ascending=False).copy()
    ranked["rank"] = np.arange(1, len(ranked) + 1)
    return ranked


def top_genes(ranked, *, n: int = 20, max_pval: float | None = 0.05):
    """Select the top-``n`` spatially variable genes to localize.

    Filters on the BH-corrected global p-value first when ``max_pval`` is set (so
    the shortlist controls the false-discovery rate across genes), then takes the
    ``n`` highest-I genes.
    """
    df = ranked
    if max_pval is not None:
        df = df[df[_corrected_pval_col(df)] <= max_pval]
    return df.sort_values("I", ascending=False).head(n)
