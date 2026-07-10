"""Normalization and highly-variable-gene selection.

The conventional scanpy transform: ``normalize_total`` (library-size
normalization) then ``log1p``, then optional HVG selection.

Why this order and these choices:
- ``normalize_total`` removes per-spot sequencing-depth differences so that a
  spatial pattern in a gene reflects biology, not how deeply a spot was
  sequenced. Without it, total-count structure (often itself spatial) leaks into
  every gene's Moran's I.
- ``log1p`` stabilizes variance; Moran's I on raw counts is dominated by a few
  high-count spots.
- HVG selection is a candidate-narrowing convenience for clustering; the
  global→local autocorrelation ranking uses Moran's I directly, not HVG.

Alternative considered: SCTransform / analytic Pearson residuals. Deferred —
heavier, and the log-normalized path is the transparent baseline this learning
project wants to be able to explain end to end.
"""

from __future__ import annotations


def normalize(adata, *, target_sum: float | None = None):
    """Library-size normalize then ``log1p`` (in place). Stores raw counts first.

    ``target_sum=None`` uses the median library size (scanpy default). The
    pre-normalization counts are kept in ``layers["counts"]`` so the raw matrix
    (needed for e.g. Moran's I on counts, or re-normalization) is never lost.
    """
    import scanpy as sc

    adata.layers["counts"] = adata.X.copy()
    sc.pp.normalize_total(adata, target_sum=target_sum)
    sc.pp.log1p(adata)
    return adata


def select_hvg(adata, *, n_top_genes: int = 2000):
    """Flag highly variable genes for the clustering path (does not subset).

    ``n_top_genes`` is capped at the number of genes present so it is safe on
    small inputs. Sets ``var["highly_variable"]``.
    """
    import scanpy as sc

    n = min(n_top_genes, adata.n_vars)
    sc.pp.highly_variable_genes(adata, n_top_genes=n)
    return adata
