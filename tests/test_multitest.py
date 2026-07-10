"""Tests for FDR correction and isolate handling (both owned seams).

Skipped until multitest is implemented.
"""

from __future__ import annotations

import numpy as np
import pytest

pytestmark = pytest.mark.skip(reason="stub — implement src/visium_spatial/multitest.py, then un-skip")


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


def _adata_with_one_detached_spot():
    """Fixture helper: synthetic lattice plus a spot placed far away (isolate)."""
    from synthetic import make_synthetic_visium

    adata = make_synthetic_visium(seed=0)
    xy = adata.obsm["spatial"].copy()
    xy[-1] = xy[:-1].max(axis=0) + 10_000  # last spot: unreachable
    adata.obsm["spatial"] = xy
    return adata
