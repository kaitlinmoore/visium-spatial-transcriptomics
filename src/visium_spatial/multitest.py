"""Multiple testing and isolate handling — both first-class, neither silent.

Two failure modes that a naive LISA/Gi* pass hides:

1. **Many spots x many genes.** Each gene's local pass produces one pseudo
   p-value per spot; calling every small p a "hotspot" guarantees false
   positives. We FDR-correct (Benjamini-Hochberg) the per-spot p-values
   *within a gene*. Correcting within a gene keeps the correction interpretable;
   comparing significance *across* genes is treated cautiously and flagged,
   because the genes have different marginal distributions and the LISA null is
   conditional (CLAUDE.md).

2. **Isolates / islands.** Spots on the tissue boundary can end up with no
   neighbors. A no-neighbor spot has an undefined spatial lag; left alone, esda
   still returns a result for it — verified: a disconnected spot comes back as a
   *significant* Low-Low cluster (q=3, p_sim=0.001), a pure artifact. Passing
   esda's ``island_weight`` (threaded through local_autocorr's wrappers) is the
   visible policy choice, but it does NOT suppress that spurious significance on
   its own. The decisive step is ``mask_isolates``: set isolate p-values to NaN
   *before* FDR, so they are excluded from the BH denominator and label as "ns"
   rather than a fake hotspot. Isolates are counted and located (``find_isolates``
   / ``isolate_report``), never deleted.

Alternatives considered:
- Bonferroni instead of BH-FDR: far too conservative for thousands of spots;
  BH controls the false-discovery rate, which is the right target for
  hotspot discovery.
- Dropping isolates: convenient but dishonest — it silently changes ``n`` and
  the tissue footprint, and hides how much boundary there is.
"""

from __future__ import annotations

import numpy as np


def fdr_within_gene(pvals: np.ndarray, *, alpha: float = 0.05):
    """Benjamini-Hochberg FDR correction over one gene's per-spot p-values.

    NaN p-values (isolates with no defined lag) are excluded from the correction:
    they do not count toward ``n``, get ``qval = NaN``, and are never rejected.
    Counting them would silently weaken every other spot's significance.

    Returns
    -------
    (rejected, qvals): boolean mask of significant spots and BH-adjusted
    p-values (monotone in p-order, ``>= pvals``, clipped to 1; NaN preserved).
    """
    p = np.asarray(pvals, dtype=float)
    qvals = np.full(p.shape, np.nan, dtype=float)
    rejected = np.zeros(p.shape, dtype=bool)

    idx = np.flatnonzero(np.isfinite(p))  # isolates (NaN) are not tests
    n = idx.size
    if n == 0:
        return rejected, qvals

    order = idx[np.argsort(p[idx])]       # positions sorted by ascending p
    ranked = p[order]
    q_raw = ranked * n / np.arange(1, n + 1)
    # enforce monotonicity: adjusted p is the running min from the largest p down
    q_mono = np.clip(np.minimum.accumulate(q_raw[::-1])[::-1], 0.0, 1.0)

    qvals[order] = q_mono
    rejected[order] = q_mono <= alpha
    return rejected, qvals


def find_isolates(W) -> np.ndarray:
    """Integer indices of spots with no neighbors in ``W`` (islands).

    A reporting primitive: the count and location of isolates is a data-quality
    metric surfaced in the results, not a filter applied behind the scenes.
    Indices are positions in ``W.id_order`` (== adata row order via build_graph).
    """
    pos = {id_: k for k, id_ in enumerate(W.id_order)}
    return np.array(sorted(pos[i] for i in W.islands), dtype=int)


def isolate_report(W) -> dict:
    """Summarize isolates: count, fraction, total, and their spot indices."""
    iso = find_isolates(W)
    n_total = int(W.n)
    return {
        "n_isolates": int(iso.size),
        "n_total": n_total,
        "fraction": float(iso.size / n_total) if n_total else 0.0,
        "isolate_indices": iso,
    }


def mask_isolates(pvals: np.ndarray, W) -> np.ndarray:
    """Return a copy of ``pvals`` with isolate positions set to NaN.

    This is what makes isolate handling honest end to end: NaN'd p-values are
    excluded from :func:`fdr_within_gene` and fall through to "ns" in
    ``local_autocorr.lisa_quadrants``, so a no-neighbor spot cannot surface as a
    spurious hotspot. The original array is not mutated.
    """
    p = np.array(pvals, dtype=float, copy=True)
    p[find_isolates(W)] = np.nan
    return p
