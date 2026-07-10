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

# LISA quadrant codes as returned by esda.Moran_Local.q
QUADRANT_LABELS = {1: "HH", 2: "LH", 3: "LL", 4: "HL"}


def local_moran(x, W, *, n_perm: int = 999, seed: int = 0):
    """Local Moran's I / LISA for one gene's values ``x`` on weights ``W``.

    Thin wrapper over ``esda.Moran_Local``. Returns the fitted object exposing
    per-spot ``Is``, quadrant codes ``q``, and pseudo p-values ``p_sim``.
    """
    raise NotImplementedError


def lisa_quadrants(local_moran_result, *, p_thresh: float = 0.05, labels=None):
    """Map an esda LISA result to per-spot quadrant labels.

    Spots whose (already FDR-corrected) p-value exceeds ``p_thresh`` are labeled
    ``"ns"`` (not significant); the rest get HH/LH/LL/HL from
    :data:`QUADRANT_LABELS`. Pass corrected p-values via ``labels``/upstream so
    this stays a pure sign/threshold map that is easy to unit-test.
    """
    raise NotImplementedError


def getis_ord_gi(x, W, *, star: bool = True, n_perm: int = 999, seed: int = 0):
    """Getis-Ord Gi* hotspot statistic for ``x`` on ``W`` (``esda.G_Local``).

    ``star=True`` includes the focal spot (Gi*). Returns the fitted object with
    per-spot z-scores and pseudo p-values.
    """
    raise NotImplementedError
