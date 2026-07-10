"""From-scratch global Moran's I

This is the statistic the project actually implements, as opposed to calls into
squidpy/esda. Its job is defensibility: reproduce squidpy's and esda's global
Moran's I on the *same* weights, within tolerance, before the local layer
(LISA / Gi*) is trusted.

Moran's I measures spatial autocorrelation. Do nearby spots carry similar
values? For a variable ``x`` on ``n`` spots with spatial weights ``W = [w_ij]``:

    z_i = x_i - mean(x)
    I   = (n / S0) * (sum_ij w_ij z_i z_j) / (sum_i z_i^2)

where ``S0 = sum_ij w_ij``. With row-standardized weights (each row sums to 1),
``S0 = n`` and the leading factor collapses to 1, so I is just the ratio of the
spatially-lagged cross-product to the total variance. That row-standardization
is why build_graph.py uses ``transform="r"``. It also gives the analytic
p-values a defined null.

Inference: expected value under no spatial structure is ``E[I] = -1/(n-1)``
(slightly negative, not zero). I test observed I against a **permutation null**
— shuffle ``x`` across spots many times, recompute I, and place the observed
value in that distribution. Permutation is assumption-light and matches how esda
reports LISA significance, keeping global and local inference consistent.

Alternatives considered:
- Analytic normality/randomization z-scores (esda offers these): fast, but lean
  on distributional assumptions that Visium counts violate; kept only as a
  cross-check.
- Geary's C: sensitive to local differences rather than the global covariance;
  Moran's I is the standard first move and the one the local indicators extend.
"""

from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np
import scipy.sparse as sp

from visium_spatial.types import SpatialWeights


def row_standardize(W: SpatialWeights) -> np.ndarray:
    """Return ``W`` with each non-empty row scaled to sum to 1.

    Turns the spatial lag from a sum over neighbors into a weighted average, so
    every spot's neighborhood contributes with equal total weight regardless of
    degree, and (for non-isolate rows) S0 = n. This is the from-scratch mirror of
    libpysal's ``transform="r"`` used on the real pipeline.

    Rows for isolates (no neighbors) sum to 0 and are left as zeros — not rescued
    here. That also keeps the isolate signal legible downstream: after this,
    "row sums to 1" means has-neighbors and "row sums to 0" means isolate.

    Accepts a dense array or scipy sparse matrix; always returns a dense float
    array (this is the small-input validation path, matching ``morans_i``).
    """
    w = (W.toarray() if sp.issparse(W) else np.asarray(W)).astype(float)
    if w.ndim != 2 or w.shape[0] != w.shape[1]:
        raise ValueError(f"weights must be a square 2-D matrix, got shape {w.shape}")

    row_sums = w.sum(axis=1)
    out = np.zeros_like(w)
    has_neighbors = row_sums != 0.0  # isolates stay all-zero, never divided
    out[has_neighbors] = w[has_neighbors] / row_sums[has_neighbors, None]
    return out


def morans_i(x: np.ndarray, W: SpatialWeights) -> float:
    """Global Moran's I of vector ``x`` under weights ``W``. Textbook definition.
    Validation reference, and so coded literally. O(N^2) on purpose.

    Parameters
    ----------
    x:
        Length-``n`` values (one per spot).
    W:
        ``n x n`` spatial weights, dense array or scipy sparse (converted to
        dense internally); row-standardized upstream for analytic parity.

    Returns
    -------
    float
        The global Moran's I statistic.
    """

    x = np.asarray(x, dtype=float)
    n = x.shape[0]
    z = x - x.mean()

    # Accept dense or sparse. Work in dense (reference, tiny inputs).
    w = W.toarray() if sp.issparse(W) else np.asarray(W, dtype=float)

    # Confirm matrix w is correct shape and matches len(x).
    if w.shape != (n, n):
        raise ValueError(f"weights shape {w.shape} does not match x length {n}")

    # Moran's I excludes self-pairs (a spot is not its own neighbor). Why reject
    # rather than tolerate: squidpy's hex graph never self-loops, so a nonzero
    # diagonal signals a construction bug upstream — better to fail loudly here
    # than to silently inflate I with spurious w_ii·z_i² terms.
    if not np.allclose(np.diag(w), 0.0):
        raise ValueError("Moran's I is undefined for nonzero self-weights (w_ii != 0); diagonal must be zero.")

    cross = 0.0     # ∑_i ∑_j w_ij z_i z_j
    w_sum = 0.0     # W = ∑_i ∑_j w_ij
    for i in range(n):
        for j in range(n):
            cross += w[i, j] * z[i] * z[j]
            w_sum += w[i, j]

    denom = np.sum(z ** 2)   # ∑_i z_i²

    # Avoid denominators that result in dividing by 0.
    if w_sum == 0:
        raise ValueError("weights sum to zero (empty graph); Moran's I is undefined")
    if denom == 0:
        raise ValueError("gene is constant across spots (zero variance); Moran's I is undefined")
    
    # Why n and w_sum can disagree: n counts every spot, but w_sum (S0) omits
    # isolates' all-zero rows, while denom still includes their z_i². That is the
    # literal formula, not a bug — isolates carry variance but no spatial
    # covariance, so they enter the denominator but not the weighted cross-sum.
    return (n / w_sum) * (cross / denom)


@dataclass(frozen=True)
class MoranPermutationResult:
    """Outcome of permutation inference for global Moran's I.

    I:      observed Moran's I.
    p_sim:  pseudo p-value, esda's folded two-sided convention (floor
            ``1 / (n_perm + 1)``).
    sim:    the ``n_perm`` simulated I values under the shuffled null.
    """

    I: float
    p_sim: float
    sim: np.ndarray = field(repr=False)  # kept out of repr; it's a long array


def morans_i_permutation(
    x: np.ndarray, W: SpatialWeights, *, n_perm: int = 999, seed: int = 0
) -> MoranPermutationResult:
    """Permutation inference for global Moran's I.

    Holds the graph ``W`` fixed and shuffles ``x`` across spots ``n_perm`` times,
    recomputing I each shuffle to build the null distribution of "location does
    not matter". Returns the observed I, the pseudo p-value, and the null array.

    Parameters
    ----------
    x:
        Length-``n`` values (one per spot).
    W:
        ``n x n`` spatial weights (dense or sparse); the SAME graph is reused for
        every shuffle.
    n_perm:
        Number of shuffles. The smallest reportable p is ``1 / (n_perm + 1)``.
    seed:
        Seed for the shuffling RNG; identical ``seed`` gives identical output.
    """
    x = np.asarray(x, dtype=float)
    n = x.shape[0]

    # Observed statistic via the literal reference. This also validates W (shape,
    # zero diagonal, non-empty graph) and x (non-constant) once, up front.
    observed = float(morans_i(x, W))

    # A shuffle only reorders values, so the mean is invariant -> z is merely
    # permuted -> Sum z^2 (denom) and S0 are constant across shuffles. Only the
    # cross-product z'Wz moves. Precompute the constant leading factor, then each
    # draw is one matvec W@z (kept in W's native form so sparse stays a matvec,
    # not a Python double loop).
    z = x - x.mean()
    factor = n / float(W.sum()) / float(z @ z)

    rng = np.random.default_rng(seed)
    sim = np.empty(n_perm, dtype=float)
    for k in range(n_perm):
        zp = rng.permutation(z)
        sim[k] = factor * float(zp @ (W @ zp))

    # Pseudo p-value, esda's convention: count sims at least as extreme as the
    # observed, fold to the smaller tail (two-sided, so strong negative
    # autocorrelation is caught too), and add 1 to both numerator and denominator
    # so the observed arrangement counts as one realization (p is never 0).
    above = int(np.sum(sim >= observed))
    larger = above if (n_perm - above) >= above else n_perm - above
    p_sim = (larger + 1.0) / (n_perm + 1.0)

    return MoranPermutationResult(I=observed, p_sim=p_sim, sim=sim)
