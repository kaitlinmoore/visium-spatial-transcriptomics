"""Tests for FDR correction and isolate handling (both owned seams).

Skipped until multitest is implemented.
"""

from __future__ import annotations

import numpy as np
import pytest


def test_fdr_is_monotone_and_bounded():
    from visium_spatial.multitest import fdr_within_gene

    p = np.array([0.001, 0.01, 0.02, 0.5, 0.9])
    rejected, q = fdr_within_gene(p, alpha=0.05)
    assert q.shape == p.shape
    assert np.all(q >= p)                 # BH adjusted p-values are never smaller
    assert np.all(np.diff(q[np.argsort(p)]) >= -1e-12)  # monotone in p-order
    assert rejected.dtype == bool


def test_all_null_rejects_nothing():
    from visium_spatial.multitest import fdr_within_gene

    p = np.full(100, 0.8)
    rejected, _ = fdr_within_gene(p, alpha=0.05)
    assert rejected.sum() == 0


def test_fdr_hand_computed():
    """Three p-values, n=3: q = sort(p)*n/rank then monotone. p=[.001,.02,.9]
    -> q=[.003, .03, .9]."""
    from visium_spatial.multitest import fdr_within_gene

    rejected, q = fdr_within_gene(np.array([0.001, 0.02, 0.9]), alpha=0.05)
    assert q == pytest.approx([0.003, 0.03, 0.9])
    assert list(rejected) == [True, True, False]


def test_fdr_ignores_nan_isolates():
    """A NaN p-value (isolate) is excluded from n, gets qval NaN, never rejected;
    the finite entries are corrected as if the NaN were absent (n=3)."""
    from visium_spatial.multitest import fdr_within_gene

    rejected, q = fdr_within_gene(np.array([0.001, 0.02, np.nan, 0.9]), alpha=0.05)
    assert np.isnan(q[2]) and not rejected[2]
    assert q[[0, 1, 3]] == pytest.approx([0.003, 0.03, 0.9])  # n=3, NaN excluded


def test_fdr_matches_statsmodels():
    """Owned BH must match statsmodels' fdr_bh exactly (cross-check, like the
    squidpy/esda parity for global Moran's)."""
    sm = pytest.importorskip("statsmodels.stats.multitest")
    from visium_spatial.multitest import fdr_within_gene

    p = np.array([0.001, 0.01, 0.02, 0.5, 0.9, 0.03, 0.2])
    _, q = fdr_within_gene(p)
    assert np.allclose(q, sm.multipletests(p, alpha=0.05, method="fdr_bh")[1])


def test_isolates_are_reported_not_dropped():
    """A spot with no neighbors must be counted and located, never silently
    removed — the tissue-boundary honesty rule."""
    from visium_spatial.build_graph import build_spatial_graph, to_libpysal_weights
    from visium_spatial.multitest import find_isolates, isolate_report

    adata = _adata_with_one_detached_spot()
    build_spatial_graph(adata)
    W = to_libpysal_weights(adata)

    iso = find_isolates(W)
    report = isolate_report(W)
    assert len(iso) == report["n_isolates"] >= 1
    # n_obs is unchanged: the spot is flagged, not deleted
    assert adata.n_obs == report["n_total"]


def test_isolate_masked_to_ns_in_lisa_pipeline():
    """End-to-end isolate honesty: esda hands a no-neighbor spot a spurious
    significant label; masking its p-value to NaN before FDR turns it into 'ns'
    rather than a fake hotspot."""
    from visium_spatial.build_graph import build_spatial_graph, to_libpysal_weights
    from visium_spatial.local_autocorr import lisa_quadrants, local_moran
    from visium_spatial.multitest import fdr_within_gene, find_isolates, mask_isolates

    adata = _adata_with_one_detached_spot()
    build_spatial_graph(adata)
    W = to_libpysal_weights(adata, transform="r")
    x = np.asarray(adata[:, "GENE_HH"].X).ravel().astype(float)
    lm = local_moran(x, W, seed=0)
    iso = find_isolates(W)

    # Without masking, the island is a spurious significant cluster...
    assert lisa_quadrants(lm, p_thresh=0.05)[iso[0]] != "ns"
    # ...masking to NaN before FDR makes it 'ns'.
    _, q = fdr_within_gene(mask_isolates(lm.p_sim, W), alpha=0.05)
    assert lisa_quadrants(lm, p_thresh=0.05, pvals=q)[iso[0]] == "ns"


def _adata_with_one_detached_spot():
    """Fixture helper: synthetic lattice plus a spot placed far away (isolate)."""
    from synthetic import make_synthetic_visium

    adata = make_synthetic_visium(seed=0)
    xy = adata.obsm["spatial"].copy()
    xy[-1] = xy[:-1].max(axis=0) + 10_000  # last spot: unreachable
    adata.obsm["spatial"] = xy
    return adata
