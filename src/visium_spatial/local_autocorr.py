"""Local spatial autocorrelation: LISA (local Moran's I) and Getis-Ord Gi*.

The deliverable the global gate leads to, and the part squidpy does not provide
out of the box — it comes from PySAL (``esda``). Both run on the SAME
``libpysal.weights.W`` produced by build_graph.py.

Local Moran's I (LISA, ``esda.Moran_Local``) decomposes the global statistic
into a per-spot contribution and classifies each spot into one of four
quadrants by the sign of its own value vs. the sign of its spatial lag:

    quadrant 1  High-High  (HH)  hot core        <- the compartment signal
    quadrant 2  Low-High   (LH)  spatial outlier
    quadrant 3  Low-Low    (LL)  cold core
    quadrant 4  High-Low   (HL)  spatial outlier

High-High clusters for compartment markers are what should recover known tissue
architecture (follicle markers over follicles, T-zone over paracortex).

Getis-Ord Gi* (``esda.G_Local`` with ``star=True``) is a complementary hotspot
statistic: it asks whether a spot *and its neighbors* form a high-value cluster
(the star includes the focal spot). LISA finds both clustering and outliers and
gives the HH/LL/HL/LH typology; Gi* gives a clean hot/cold z-score. Reporting
both, on one graph, cross-checks the hotspots.

Significance is quadrant/label-gated on the per-spot pseudo p-values, which are
FDR-corrected within a gene by multitest.py before anything is called a hotspot.

Alternative considered: analytic (normal-approximation) p-values. esda's
conditional-permutation p-values are assumption-lighter and consistent with the
global permutation null; kept as the default.
"""

from __future__ import annotations

import numpy as np
from esda.getisord import G_Local
from esda.moran import Moran_Local

# LISA quadrant codes as returned by esda.Moran_Local.q
QUADRANT_LABELS = {1: "HH", 2: "LH", 3: "LL", 4: "HL"}


def local_moran(x, W, *, n_perm: int = 999, seed: int = 0, island_weight: float = 0.0):
    """Local Moran's I / LISA for one gene's values ``x`` on weights ``W``.

    Thin wrapper over ``esda.Moran_Local`` (the statistic is PySAL's, not ours).
    Returns the fitted object exposing per-spot ``Is``, quadrant codes ``q``, and
    pseudo p-values ``p_sim``. ``seed`` makes the conditional-permutation
    p-values reproducible.

    ``island_weight`` is passed to esda as the explicit isolate policy (a visible
    choice, per CLAUDE.md). Note it does NOT prevent isolates from getting a
    spurious significant label — that is handled by masking their p-values via
    ``multitest.mask_isolates`` before FDR.
    """
    x = np.asarray(x, dtype=float)
    return Moran_Local(x, W, permutations=n_perm, seed=seed, island_weight=island_weight)


def lisa_quadrants(local_moran_result, *, p_thresh: float = 0.05, pvals=None):
    """Map an esda LISA result to per-spot quadrant labels. The owned seam.

    A pure sign/threshold map: a spot significant at ``p_thresh`` gets its
    HH/LH/LL/HL label from :data:`QUADRANT_LABELS`; everything else is ``"ns"``
    (not significant). Keeping it a pure function of (quadrant, p-value) is what
    makes it unit-testable without running esda.

    By default it reads ``result.p_sim``. Pass FDR-corrected p-values via
    ``pvals`` to gate on the multiple-testing-corrected significance instead —
    this is the plug point for ``multitest.fdr_within_gene``. NaN p-values
    (islands with no defined lag) fall through to ``"ns"``.
    """
    q = np.asarray(local_moran_result.q)
    p = np.asarray(local_moran_result.p_sim if pvals is None else pvals, dtype=float)

    labels = np.empty(q.shape[0], dtype=object)
    for i in range(q.shape[0]):
        significant = p[i] <= p_thresh  # NaN -> False -> "ns"
        labels[i] = QUADRANT_LABELS.get(int(q[i]), "ns") if significant else "ns"
    return labels


def getis_ord_gi(
    x, W, *, star: bool = True, n_perm: int = 999, seed: int = 0, island_weight: float = 0.0
):
    """Getis-Ord Gi* hotspot statistic for ``x`` on ``W`` (``esda.G_Local``).

    Thin wrapper. ``star=True`` includes the focal spot (Gi* vs Gi). Returns the
    fitted object with per-spot z-scores ``Zs`` and pseudo p-values ``p_sim``.
    ``island_weight`` is the explicit isolate policy passed to esda (see
    :func:`local_moran`).
    """
    x = np.asarray(x, dtype=float)
    return G_Local(x, W, star=star, permutations=n_perm, seed=seed, island_weight=island_weight)
