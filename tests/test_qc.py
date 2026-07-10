"""Tests for spot/gene QC filtering.

Each test builds a FRESH synthetic AnnData (not the shared session fixture),
because QC mutates obs/var and subsets spots — doing that to the shared fixture
would corrupt the raw-counts invariant the autocorrelation parity tests rely on.
"""

from __future__ import annotations

import numpy as np
import pytest

from synthetic import make_synthetic_visium


def test_compute_qc_metrics_adds_columns():
    from visium_spatial.qc import compute_qc_metrics

    a = make_synthetic_visium(seed=0)
    compute_qc_metrics(a)
    for col in ("total_counts", "n_genes_by_counts", "pct_counts_mt"):
        assert col in a.obs
    # synthetic genes are not mitochondrial -> zero mito fraction everywhere
    assert (a.obs["pct_counts_mt"] == 0).all()


def test_filter_spots_applies_thresholds_and_copies():
    from visium_spatial.qc import filter_spots

    a = make_synthetic_visium(seed=0)
    n0 = a.n_obs

    # permissive keeps all; impossibly strict drops all
    assert filter_spots(a, min_counts=1, min_genes=1, max_pct_mito=100.0).n_obs == n0
    assert filter_spots(a, min_counts=10_000, min_genes=1, max_pct_mito=100.0).n_obs == 0

    # a real threshold is actually enforced on the survivors
    mid = filter_spots(a, min_counts=20, min_genes=1, max_pct_mito=100.0)
    assert (mid.obs["total_counts"] >= 20).all()

    assert a.n_obs == n0  # original untouched (a copy is returned)


def test_filter_genes_drops_rare():
    from visium_spatial.qc import filter_genes

    a = make_synthetic_visium(seed=0)
    X = a.X.copy()
    X[:, 1] = 0
    X[0, 1] = 5  # GENE_DISP now detected in a single spot
    a.X = X

    filter_genes(a, min_spots=10)
    assert "GENE_DISP" not in a.var_names
    assert a.n_vars == 2
