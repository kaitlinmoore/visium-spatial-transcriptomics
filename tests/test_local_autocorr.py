"""Tests for LISA quadrant assignment and the local indicators.

The quadrant map is an owned seam (sign of value vs. sign of spatial lag ->
HH/LH/LL/HL). Skipped until local_autocorr is implemented.
"""

from __future__ import annotations

from types import SimpleNamespace

import numpy as np
import pytest


# --- owned seam: the pure quadrant/threshold map (no esda) --------------------

def test_lisa_quadrants_pure_mapping():
    """Significant spots get their quadrant label; the rest are 'ns'."""
    from visium_spatial.local_autocorr import lisa_quadrants

    result = SimpleNamespace(q=np.array([1, 2, 3, 4]), p_sim=np.array([0.01, 0.2, 0.01, 0.9]))
    labels = lisa_quadrants(result, p_thresh=0.05)
    assert list(labels) == ["HH", "ns", "LL", "ns"]


def test_lisa_quadrants_uses_corrected_pvals_when_given():
    """Passing FDR-corrected p-values overrides result.p_sim (the multitest plug
    point). Here the second spot is demoted to 'ns' by correction."""
    from visium_spatial.local_autocorr import lisa_quadrants

    result = SimpleNamespace(q=np.array([1, 1]), p_sim=np.array([0.01, 0.01]))
    labels = lisa_quadrants(result, p_thresh=0.05, pvals=np.array([0.01, 0.20]))
    assert list(labels) == ["HH", "ns"]


def test_lisa_quadrants_island_nan_is_ns():
    """A NaN p-value (island, undefined lag) must fall through to 'ns', not crash
    on int(q) when q is undefined."""
    from visium_spatial.local_autocorr import lisa_quadrants

    result = SimpleNamespace(q=np.array([1, 0]), p_sim=np.array([0.01, np.nan]))
    labels = lisa_quadrants(result, p_thresh=0.05)
    assert list(labels) == ["HH", "ns"]


# --- integration: esda on the synthetic fixture ------------------------------

def test_planted_patch_recovers_high_high(adata, hh_indices):
    """The planted GENE_HH patch should come back as significant High-High
    (quadrant 1) — the whole compartment-recovery premise in miniature.

    Gates on RAW p_sim by design: this tests signal DETECTION (does the local
    statistic localize the planted patch?), a different property from noise
    control. BH-FDR on this 49-spot fixture is so stringent it collapses the
    patch to 1/7 regardless of permutation count, so FDR belongs on the dispersed
    (noise-control) test, not here."""
    from visium_spatial.build_graph import build_spatial_graph, to_libpysal_weights
    from visium_spatial.local_autocorr import lisa_quadrants, local_moran

    build_spatial_graph(adata)
    W = to_libpysal_weights(adata, transform="r")
    lm = local_moran(np.asarray(adata[:, "GENE_HH"].X).ravel().astype(float), W, seed=0)
    labels = lisa_quadrants(lm, p_thresh=0.05)

    patch_labels = [labels[i] for i in hh_indices]
    assert patch_labels.count("HH") >= len(hh_indices) - 1  # allow one edge spot


def test_dispersed_gene_has_few_significant_clusters(adata):
    """Noise-control property: under BH-FDR a spatially random gene yields (near)
    no significant clusters. Gates on FDR-corrected p via pvals= — the counterpart
    to the patch test, which gates on raw p_sim for signal detection. On raw
    p_sim this gene admits ~6/49 false positives; FDR drives it to 0."""
    from visium_spatial.build_graph import build_spatial_graph, to_libpysal_weights
    from visium_spatial.local_autocorr import lisa_quadrants, local_moran
    from visium_spatial.multitest import fdr_within_gene

    build_spatial_graph(adata)
    W = to_libpysal_weights(adata, transform="r")
    lm = local_moran(np.asarray(adata[:, "GENE_DISP"].X).ravel().astype(float), W, seed=0)
    _, q = fdr_within_gene(lm.p_sim, alpha=0.05)
    labels = lisa_quadrants(lm, p_thresh=0.05, pvals=q)
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
