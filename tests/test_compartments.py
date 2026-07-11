"""Tests for the compartment-recovery concordance seam.

The concordance metric (Jaccard of HH spot sets) is owned logic, tested as a pure
function; hh_spot_set is checked end-to-end on the planted synthetic patch.
"""

from __future__ import annotations

import numpy as np
import pytest


def test_jaccard_values():
    from visium_spatial.compartments import jaccard

    assert jaccard({1, 2, 3}, {2, 3, 4}) == pytest.approx(2 / 4)  # {2,3} / {1,2,3,4}
    assert jaccard({1, 2}, {3, 4}) == 0.0
    assert jaccard({1, 2}, {1, 2}) == 1.0
    assert np.isnan(jaccard(set(), set()))  # both empty -> undefined


def test_concordance_matrix_block_structure():
    """Two same-compartment sets overlap; a disjoint set does not; diagonal is 1."""
    from visium_spatial.compartments import concordance_matrix

    sets = {"A1": {1, 2, 3}, "A2": {2, 3, 4}, "B1": {10, 11}}
    names, M = concordance_matrix(sets)
    ia1, ia2, ib1 = names.index("A1"), names.index("A2"), names.index("B1")

    assert M[ia1, ia1] == 1.0
    assert M[ia1, ia2] == pytest.approx(0.5) == M[ia2, ia1]  # within-compartment
    assert M[ia1, ib1] == 0.0                                # between-compartment
    assert M.shape == (3, 3)


def test_hh_spot_set_recovers_patch(adata, hh_indices):
    """On the planted GENE_HH patch, the significant HH spots are a subset of the
    planted patch (no HH called outside it)."""
    from visium_spatial.build_graph import build_spatial_graph, to_libpysal_weights
    from visium_spatial.compartments import hh_spot_set
    from visium_spatial.local_autocorr import local_moran

    build_spatial_graph(adata)
    W = to_libpysal_weights(adata, transform="r")
    lm = local_moran(np.asarray(adata[:, "GENE_HH"].X).ravel().astype(float), W, seed=0)

    hh = hh_spot_set(lm, W)
    patch = {int(i) for i in hh_indices}
    assert len(hh) >= 1
    assert hh <= patch  # every HH spot lies within the planted patch
