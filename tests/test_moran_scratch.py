"""Tests for the from-scratch global Moran's I (the owned defensibility check).

The ``morans_i`` tests below are self-contained: each expected value is computed
by hand (shown in the test's comment) or against numpy's closed form, never
against esda/squidpy. That is deliberate — this function EXISTS to cross-check
the libraries, so it must be pinned to first principles, not to them.

Scaffolds for functions that aren't implemented yet (``row_standardize``, the
``adata``/``build_graph`` recovery tests, permutation inference) are kept but
skipped individually; un-skip each as its function lands.
"""

from __future__ import annotations

import numpy as np
import pytest
import scipy.sparse as sp

from visium_spatial.moran_scratch import morans_i


# --- hand-computed point values ---------------------------------------------

def test_morans_i_triangle_hand_computed():
    """Triangle (all spots mutually adjacent), binary weights, x = [1, 2, 3].

    z = [-1, 0, 1]. The only neighbor pair with both z != 0 is spots 0 and 2,
    adjacent in both directions: cross = 2 * (-1)(1) = -2. S0 = 6 (six directed
    edges), denom = sum z^2 = 2, n = 3.
        I = (3/6) * (-2/2) = -0.5
    Negative is correct: the extreme values (1 and 3) are neighbors, so like
    sits next to unlike.
    """
    W = np.array([[0, 1, 1], [1, 0, 1], [1, 1, 0]], dtype=float)
    x = np.array([1.0, 2.0, 3.0])
    assert morans_i(x, W) == pytest.approx(-0.5)


def test_morans_i_line_positive_hand_computed():
    """4-spot line 0-1-2-3, binary weights, x = [1, 1, 2, 2].

    z = [-.5, -.5, .5, .5]. Directed edges contribute (each counted both ways):
    (0,1)=+.25, (1,2)=-.25, (2,3)=+.25 -> cross = 2*(.25-.25+.25) = 0.5.
    S0 = 6, denom = sum z^2 = 1, n = 4.
        I = (4/6) * (0.5/1) = 1/3
    Positive is correct: low values cluster left, high values right.
    """
    W = np.array(
        [[0, 1, 0, 0], [1, 0, 1, 0], [0, 1, 0, 1], [0, 0, 1, 0]], dtype=float
    )
    x = np.array([1.0, 1.0, 2.0, 2.0])
    assert morans_i(x, W) == pytest.approx(1.0 / 3.0)


# --- differential check against numpy's closed form -------------------------

def test_morans_i_matches_quadratic_form():
    """The nested loop must equal the matrix expression (n/S0)*(z'Wz)/(z'z) on a
    non-trivial symmetric, zero-diagonal weight matrix. Guards against index
    transposition bugs that the two point-values above could miss (both are
    symmetric in a way that hides i<->j swaps)."""
    W = np.array(
        [
            [0.0, 2.0, 0.5, 0.0],
            [2.0, 0.0, 1.0, 3.0],
            [0.5, 1.0, 0.0, 0.0],
            [0.0, 3.0, 0.0, 0.0],
        ]
    )
    x = np.array([4.0, 1.0, 7.0, 2.0])
    z = x - x.mean()
    expected = (len(x) / W.sum()) * (z @ W @ z) / (z @ z)
    assert morans_i(x, W) == pytest.approx(expected)


def test_morans_i_invariant_to_weight_scaling():
    """Moran's I is unchanged if every weight is scaled by a constant, because
    S0 scales with the numerator. Here the row-standardized triangle (every
    weight 1/2) gives the same -0.5 as the binary triangle. Built with literal
    0.5s so this does not depend on row_standardize (still a stub)."""
    W_binary = np.array([[0, 1, 1], [1, 0, 1], [1, 1, 0]], dtype=float)
    W_rowstd = W_binary / 2.0  # each row already sums to 2 -> becomes 1
    x = np.array([1.0, 2.0, 3.0])
    assert morans_i(x, W_rowstd) == pytest.approx(morans_i(x, W_binary))
    assert morans_i(x, W_rowstd) == pytest.approx(-0.5)


# --- input handling ----------------------------------------------------------

def test_morans_i_accepts_sparse_weights():
    """Real weights arrive sparse from squidpy's obsp; result must match dense."""
    W = np.array([[0, 1, 1], [1, 0, 1], [1, 1, 0]], dtype=float)
    x = np.array([1.0, 2.0, 3.0])
    assert morans_i(x, sp.csr_matrix(W)) == pytest.approx(morans_i(x, W))


# --- error paths (degenerate inputs are metrics, not silent NaNs) -----------

def test_rejects_shape_mismatch():
    x = np.array([1.0, 2.0, 3.0])
    W = np.zeros((2, 2))
    with pytest.raises(ValueError, match="shape"):
        morans_i(x, W)


def test_rejects_nonzero_diagonal():
    """A self-weight (w_ii != 0) signals a graph-construction bug upstream."""
    x = np.array([1.0, 2.0, 3.0])
    W = np.array([[1, 1, 0], [1, 0, 1], [0, 1, 0]], dtype=float)
    with pytest.raises(ValueError, match="diagonal"):
        morans_i(x, W)


def test_rejects_empty_graph():
    """No edges -> S0 = 0 -> I undefined; must raise, not divide by zero."""
    x = np.array([1.0, 2.0, 3.0])
    W = np.zeros((3, 3))
    with pytest.raises(ValueError, match="zero"):
        morans_i(x, W)


def test_rejects_constant_gene():
    """Zero variance -> denom = 0 -> I undefined."""
    x = np.array([2.0, 2.0, 2.0])
    W = np.array([[0, 1, 1], [1, 0, 1], [1, 1, 0]], dtype=float)
    with pytest.raises(ValueError, match="constant"):
        morans_i(x, W)


# --- permutation inference ---------------------------------------------------

def _chain_weights(n: int) -> np.ndarray:
    """Binary adjacency of a 1-D chain 0-1-...-(n-1). A monotone x on this graph
    has strong positive autocorrelation, which makes a clean signal fixture."""
    W = np.zeros((n, n))
    for i in range(n - 1):
        W[i, i + 1] = W[i + 1, i] = 1.0
    return W


def test_permutation_is_deterministic_under_seed():
    from visium_spatial.moran_scratch import morans_i_permutation

    W = _chain_weights(10)
    x = np.arange(10, dtype=float)
    r1 = morans_i_permutation(x, W, n_perm=200, seed=7)
    r2 = morans_i_permutation(x, W, n_perm=200, seed=7)
    assert r1.p_sim == pytest.approx(r2.p_sim)
    assert np.array_equal(r1.sim, r2.sim)


def test_permutation_observed_matches_morans_i():
    from visium_spatial.moran_scratch import morans_i_permutation

    W = _chain_weights(8)
    x = np.array([3.0, 1.0, 4.0, 1.0, 5.0, 9.0, 2.0, 6.0])
    r = morans_i_permutation(x, W, n_perm=50, seed=0)
    assert r.I == pytest.approx(morans_i(x, W))
    assert r.sim.shape == (50,)


def test_permutation_detects_strong_autocorrelation():
    """A smooth gradient on a chain is near-maximally clustered, so almost no
    shuffle beats it -> p_sim sits at the 1/(n_perm+1) floor."""
    from visium_spatial.moran_scratch import morans_i_permutation

    W = _chain_weights(15)
    x = np.arange(15, dtype=float)
    r = morans_i_permutation(x, W, n_perm=999, seed=0)
    assert r.I > 0.5
    assert r.p_sim <= 0.01


def test_permutation_null_mean_is_minus_one_over_n_minus_one():
    """The permutation null of Moran's I has exact mean E[I] = -1/(n-1),
    independent of the data. With many shuffles the sample mean recovers it."""
    from visium_spatial.moran_scratch import morans_i_permutation

    n = 12
    W = _chain_weights(n)
    x = np.arange(n, dtype=float)
    r = morans_i_permutation(x, W, n_perm=9999, seed=0)
    assert r.sim.mean() == pytest.approx(-1.0 / (n - 1), abs=0.03)


def test_permutation_pseudo_p_within_bounds():
    from visium_spatial.moran_scratch import morans_i_permutation

    W = _chain_weights(10)
    x = np.arange(10, dtype=float)
    r = morans_i_permutation(x, W, n_perm=999, seed=1)
    assert 1.0 / (999 + 1) <= r.p_sim <= 1.0


def test_permutation_accepts_sparse_weights():
    from visium_spatial.moran_scratch import morans_i_permutation

    W = _chain_weights(10)
    x = np.arange(10, dtype=float)
    dense = morans_i_permutation(x, W, n_perm=200, seed=3)
    sparse = morans_i_permutation(x, sp.csr_matrix(W), n_perm=200, seed=3)
    assert sparse.p_sim == pytest.approx(dense.p_sim)
    assert np.allclose(sparse.sim, dense.sim)


# --- scaffolds for not-yet-implemented functions (un-skip as they land) ------

def test_row_standardize_rows_sum_to_one():
    from visium_spatial.moran_scratch import row_standardize

    W = np.array([[0, 1, 1], [1, 0, 1], [1, 1, 0]], dtype=float)
    Wr = row_standardize(W)
    assert np.allclose(Wr.sum(axis=1), 1.0)


def test_row_standardize_values():
    """Row [0, 1, 1] (two equal neighbors) -> [0, .5, .5]."""
    from visium_spatial.moran_scratch import row_standardize

    W = np.array([[0, 1, 1], [3, 0, 1], [1, 1, 0]], dtype=float)
    Wr = row_standardize(W)
    assert Wr[0] == pytest.approx([0.0, 0.5, 0.5])
    assert Wr[1] == pytest.approx([0.75, 0.0, 0.25])  # 3 and 1 -> 3/4 and 1/4


def test_row_standardize_leaves_isolates_as_zero():
    """An all-zero row (isolate) is left untouched, not divided by zero. Its row
    stays 0 while the other rows still sum to 1."""
    from visium_spatial.moran_scratch import row_standardize

    W = np.array([[0, 1, 0], [1, 0, 0], [0, 0, 0]], dtype=float)  # spot 2 isolated
    Wr = row_standardize(W)
    assert np.all(Wr[2] == 0.0)
    assert Wr[:2].sum(axis=1) == pytest.approx([1.0, 1.0])
    assert np.isfinite(Wr).all()  # no nan/inf from a 0/0


def test_row_standardize_accepts_sparse():
    from visium_spatial.moran_scratch import row_standardize

    W = np.array([[0, 1, 1], [1, 0, 1], [1, 1, 0]], dtype=float)
    Wr = row_standardize(sp.csr_matrix(W))
    assert Wr == pytest.approx(row_standardize(W))


def test_row_standardize_feeds_morans_i_consistently():
    """morans_i computes S0 itself, so hand-standardizing the weights first must
    give the identical statistic — the two isolate/standardization stories agree."""
    from visium_spatial.moran_scratch import row_standardize

    W = np.array([[0, 1, 1], [1, 0, 1], [1, 1, 0]], dtype=float)
    x = np.array([1.0, 2.0, 3.0])
    assert morans_i(x, row_standardize(W)) == pytest.approx(morans_i(x, W))


def test_planted_cluster_is_positive(adata):
    from visium_spatial.build_graph import build_spatial_graph, to_libpysal_weights
    from visium_spatial.moran_scratch import morans_i as _morans_i

    build_spatial_graph(adata)
    W = to_libpysal_weights(adata, transform="r").full()[0]
    x = np.asarray(adata[:, "GENE_HH"].X).ravel().astype(float)
    assert _morans_i(x, W) > 0


def test_dispersed_gene_near_zero(adata):
    from visium_spatial.build_graph import build_spatial_graph, to_libpysal_weights
    from visium_spatial.moran_scratch import morans_i as _morans_i

    build_spatial_graph(adata)
    W = to_libpysal_weights(adata, transform="r").full()[0]
    x = np.asarray(adata[:, "GENE_DISP"].X).ravel().astype(float)
    assert abs(_morans_i(x, W)) < 0.2
