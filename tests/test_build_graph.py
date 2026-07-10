"""Tests for the hex graph build and the squidpy->libpysal weights bridge.

The bridge is an owned seam: the SAME neighbor set must reach both squidpy and
esda. Skipped until build_graph is implemented.
"""

from __future__ import annotations

import anndata as ad
import numpy as np
import pytest
import scipy.sparse as sp


def test_interior_spots_have_six_neighbors(adata):
    from visium_spatial.build_graph import assert_hex_neighbors, build_spatial_graph

    build_spatial_graph(adata, n_neighs=6, coord_type="grid")
    report = assert_hex_neighbors(adata, expected=6)
    assert report["n_isolates"] == 0
    assert report["mode_degree"] == 6


def test_bridge_preserves_squidpy_neighbors(adata):
    """libpysal W must encode exactly obsp['spatial_connectivities']."""
    from visium_spatial.build_graph import build_spatial_graph, to_libpysal_weights

    build_spatial_graph(adata)
    W = to_libpysal_weights(adata, transform="O")  # untransformed for adjacency check
    conn = adata.obsp["spatial_connectivities"]
    for i, name in enumerate(adata.obs_names):
        squidpy_nbrs = set(np.flatnonzero(conn[i].toarray().ravel()))
        pysal_nbrs = {list(adata.obs_names).index(j) for j in W.neighbors[name]}
        assert squidpy_nbrs == pysal_nbrs


def test_row_standardized_rows_sum_to_one(adata):
    from visium_spatial.build_graph import build_spatial_graph, to_libpysal_weights

    build_spatial_graph(adata)
    W = to_libpysal_weights(adata, transform="r")
    for name in adata.obs_names:
        w = W.weights[name]
        if w:  # non-isolate
            assert sum(w) == pytest.approx(1.0)


def _adata_from_adjacency(A: np.ndarray):
    """Minimal AnnData carrying a hand-built connectivity graph in obsp."""
    a = ad.AnnData(np.zeros((A.shape[0], 1), dtype=float))
    a.obsp["spatial_connectivities"] = sp.csr_matrix(A.astype(float))
    return a


def test_isolates_are_preserved_not_dropped():
    """A disconnected spot must survive the bridge as an island, not vanish."""
    from visium_spatial.build_graph import to_libpysal_weights

    # spots 0-1 connected, spot 2 isolated
    a = _adata_from_adjacency(np.array([[0, 1, 0], [1, 0, 0], [0, 0, 0]]))
    W = to_libpysal_weights(a, transform="r")
    names = list(a.obs_names)
    assert names[2] in W.islands
    assert sum(W.weights[names[0]]) == pytest.approx(1.0)  # others still standardized


def test_assert_hex_neighbors_raises_on_wrong_graph():
    """A ring (every spot degree 2) is not the hex neighborhood -> must raise."""
    from visium_spatial.build_graph import assert_hex_neighbors

    n = 6
    A = np.zeros((n, n))
    for i in range(n):
        A[i, (i + 1) % n] = A[(i + 1) % n, i] = 1
    a = _adata_from_adjacency(A)
    with pytest.raises(ValueError, match="modal degree"):
        assert_hex_neighbors(a, expected=6)


def test_scratch_matches_squidpy_and_esda_global_moran(adata):
    """The project's core defensibility check: from-scratch Moran's I reproduces
    BOTH squidpy and esda on the SAME row-standardized weights. All three
    row-standardize, so parity holds only with transform='r' (verified: binary
    weights give a different I)."""
    import squidpy as sq
    from esda.moran import Moran

    from visium_spatial.build_graph import build_spatial_graph, to_libpysal_weights
    from visium_spatial.moran_scratch import morans_i

    build_spatial_graph(adata)
    W = to_libpysal_weights(adata, transform="r")
    W_dense = W.full()[0]
    x = np.asarray(adata[:, "GENE_GRAD"].X).ravel().astype(float)

    sq.gr.spatial_autocorr(adata, mode="moran", genes=["GENE_GRAD"], n_perms=None, seed=0)
    squidpy_I = float(adata.uns["moranI"].loc["GENE_GRAD", "I"])
    esda_I = float(Moran(x, W, permutations=0).I)
    scratch_I = morans_i(x, W_dense)

    assert scratch_I == pytest.approx(squidpy_I, abs=1e-6)
    assert scratch_I == pytest.approx(esda_I, abs=1e-6)
