"""Tests for normalization and HVG flagging (fresh AnnData per test)."""

from __future__ import annotations

import numpy as np
import pytest

from synthetic import make_synthetic_visium


def test_normalize_preserves_raw_counts_and_logs():
    from visium_spatial.preprocess import normalize

    a = make_synthetic_visium(seed=0)
    raw_max = float(np.asarray(a.X).max())

    normalize(a)
    assert "counts" in a.layers
    assert float(np.asarray(a.layers["counts"]).max()) == pytest.approx(raw_max)  # raw kept
    assert float(np.asarray(a.X).max()) < raw_max                                 # log1p shrank it
    assert "log1p" in a.uns                                                        # scanpy records it


def test_select_hvg_flags_without_subsetting():
    from visium_spatial.preprocess import normalize, select_hvg

    a = make_synthetic_visium(seed=0)
    normalize(a)
    n_vars = a.n_vars

    select_hvg(a, n_top_genes=2000)  # capped to n_vars internally, safe on 3 genes
    assert "highly_variable" in a.var
    assert a.n_vars == n_vars                       # flags, does not subset
    assert 0 < a.var["highly_variable"].sum() <= n_vars
