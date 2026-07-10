"""Tests for LISA quadrant assignment and the local indicators.

The quadrant map is an owned seam (sign of value vs. sign of spatial lag ->
HH/LH/LL/HL). Skipped until local_autocorr is implemented.
"""

from __future__ import annotations

import pytest

pytestmark = pytest.mark.skip(reason="stub — implement src/visium_spatial/local_autocorr.py, then un-skip")


def test_planted_patch_recovers_high_high(adata, hh_indices):
    """The planted GENE_HH patch should come back as significant High-High
    (quadrant 1) — the whole compartment-recovery premise in miniature."""
    from visium_spatial.build_graph import build_spatial_graph, to_libpysal_weights
    from visium_spatial.local_autocorr import lisa_quadrants, local_moran

    build_spatial_graph(adata)
    W = to_libpysal_weights(adata, transform="r")
    lm = local_moran(adata[:, "GENE_HH"].X.ravel(), W, seed=0)
    labels = lisa_quadrants(lm, p_thresh=0.05)

    patch_labels = [labels[i] for i in hh_indices]
    assert patch_labels.count("HH") >= len(hh_indices) - 1  # allow one edge spot


def test_dispersed_gene_has_few_significant_clusters(adata):
    from visium_spatial.build_graph import build_spatial_graph, to_libpysal_weights
    from visium_spatial.local_autocorr import lisa_quadrants, local_moran

    build_spatial_graph(adata)
    W = to_libpysal_weights(adata, transform="r")
    lm = local_moran(adata[:, "GENE_DISP"].X.ravel(), W, seed=0)
    labels = lisa_quadrants(lm, p_thresh=0.05)
    n_sig = sum(1 for v in labels if v != "ns")
    assert n_sig <= 0.1 * adata.n_obs


def test_getis_gi_hotspot_over_patch(adata, hh_indices):
    from visium_spatial.build_graph import build_spatial_graph, to_libpysal_weights
    from visium_spatial.local_autocorr import getis_ord_gi

    build_spatial_graph(adata)
    W = to_libpysal_weights(adata, transform="r")
    gi = getis_ord_gi(adata[:, "GENE_HH"].X.ravel(), W, star=True, seed=0)
    # focal spots in the patch should be hot (positive z)
    assert (gi.Zs[hh_indices] > 0).mean() > 0.8
