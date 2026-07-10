"""Tests for the global Moran's I ranking gate (squidpy) and top-gene selection.

The ranking is a thin squidpy wrapper, so the load-bearing test is defensibility:
the squidpy I that drives the ranking must agree with the from-scratch morans_i
on the SAME weights (CLAUDE.md), otherwise the whole "rank then localize" premise
is built on an unverified number.
"""

from __future__ import annotations

import numpy as np
import pytest


def test_rank_svgs_orders_by_moran(adata):
    from visium_spatial.build_graph import build_spatial_graph
    from visium_spatial.global_autocorr import rank_svgs

    build_spatial_graph(adata)
    ranked = rank_svgs(adata, seed=0)

    # sorted by I descending, with a 1..n rank column
    assert list(ranked["I"]) == sorted(ranked["I"], reverse=True)
    assert list(ranked["rank"]) == list(range(1, len(ranked) + 1))
    # the structured genes outrank the dispersed one
    assert ranked.index[0] in {"GENE_GRAD", "GENE_HH"}
    assert ranked.index[-1] == "GENE_DISP"


def test_rank_svgs_matches_scratch(adata):
    """Defensibility: squidpy's per-gene I (which drives the rank) equals the
    from-scratch morans_i on the row-standardized weights, and the two produce
    the same gene order."""
    from visium_spatial.build_graph import build_spatial_graph, to_libpysal_weights
    from visium_spatial.global_autocorr import rank_svgs
    from visium_spatial.moran_scratch import morans_i

    build_spatial_graph(adata)
    ranked = rank_svgs(adata, seed=0)

    W = to_libpysal_weights(adata, transform="r").full()[0]
    scratch = {
        g: morans_i(np.asarray(adata[:, g].X).ravel().astype(float), W)
        for g in adata.var_names
    }
    for g in adata.var_names:
        assert ranked.loc[g, "I"] == pytest.approx(scratch[g], abs=1e-6)

    scratch_order = sorted(scratch, key=scratch.get, reverse=True)
    assert list(ranked.index) == scratch_order


def test_top_genes_filters_and_limits(adata):
    from visium_spatial.build_graph import build_spatial_graph
    from visium_spatial.global_autocorr import rank_svgs, top_genes

    build_spatial_graph(adata)
    ranked = rank_svgs(adata, seed=0)

    # n caps the count and keeps the highest-I gene
    top1 = top_genes(ranked, n=1, max_pval=0.05)
    assert len(top1) == 1
    assert top1.index[0] == ranked.index[0]

    # the dispersed gene is not significant -> filtered out at max_pval=0.05
    kept = top_genes(ranked, n=10, max_pval=0.05)
    assert "GENE_DISP" not in kept.index

    # max_pval=None disables the filter -> all genes eligible
    assert len(top_genes(ranked, n=10, max_pval=None)) == len(ranked)
