"""Quantified compartment-recovery validation.

The section ships with **no ground-truth compartment annotation**, so "recovery"
cannot be a supervised accuracy. It is assessed as *concordance* of High-High
LISA hotspots:

- markers of the **same** compartment should share HH spots (they are co-expressed
  in the same tissue region);
- markers of **different** compartments (B-cell follicle vs T-cell paracortex)
  should be spatially **segregated** (low HH overlap).

A block structure — high within-compartment, low between-compartment overlap —
is the evidence that the local statistics recover known lymph-node architecture.
The overlap is a number (Jaccard of HH spot sets), not an eyeballed overlay.

Honest caveats (see docs/methodology.md §8): within-compartment overlap is partly
trivial (co-expression), so the *stronger* signal is between-compartment
segregation; and each spot is multiple cells, so "compartment" here is a
dominant-signal notion, not a cell-level label.
"""

from __future__ import annotations

import numpy as np

# Published lymph-node marker -> tissue-domain panel, used instead of an
# ad-hoc marker list so the segregation/nesting argument rests on a citable
# source (see docs/methodology.md §8 for the reference). Genes verified present
# and spatially autocorrelated in V1_Human_Lymph_Node. The four *compact*
# compartments (follicle, germinal_center, T-zone, medulla) are where HH-hotspot
# concordance is meaningful; the thin/linear domains (lymphatic, blood_vessel)
# and the transitional B-T interface are included for completeness but do not
# form clean HH blocks.
LYMPH_NODE_MARKERS = {
    "T-zone": ["TRBC1", "TRAC"],
    "follicle": ["FDCSP", "CR2"],
    "germinal_center": ["BCL6", "MYBL1"],
    "B-T_interface": ["THY1"],
    "medulla": ["IGHG1", "IGHG2"],
    "lymphatic": ["LYVE1", "PROX1"],
    "blood_vessel": ["VWF", "PECAM1"],
}

# Compartments compact enough for the HH-hotspot concordance test to be meaningful.
CORE_COMPARTMENTS = ("follicle", "germinal_center", "T-zone", "medulla")


def hh_spot_set(local_moran_result, W, *, p_thresh: float = 0.05) -> set[int]:
    """Indices of significant High-High spots for one gene.

    Composes the honest hotspot pipeline: isolate p-values masked to NaN, FDR
    within the gene, then the LISA quadrant map — and keeps only spots labelled
    ``"HH"`` (high value surrounded by high values, significant after correction).
    """
    from visium_spatial.local_autocorr import lisa_quadrants
    from visium_spatial.multitest import fdr_within_gene, mask_isolates

    _, qvals = fdr_within_gene(mask_isolates(local_moran_result.p_sim, W))
    labels = lisa_quadrants(local_moran_result, p_thresh=p_thresh, pvals=qvals)
    return set(np.flatnonzero(labels == "HH").tolist())


def jaccard(a, b) -> float:
    """Jaccard overlap |a∩b| / |a∪b| of two spot-index sets.

    Returns NaN when both sets are empty (overlap undefined), so an all-``ns``
    marker does not masquerade as either concordant or segregated.
    """
    a, b = set(a), set(b)
    union = a | b
    if not union:
        return float("nan")
    return len(a & b) / len(union)


def concordance_matrix(hh_sets: dict[str, set[int]]):
    """Pairwise Jaccard overlap of HH spot sets.

    Parameters
    ----------
    hh_sets:
        Mapping of marker name -> its HH spot-index set (from :func:`hh_spot_set`).

    Returns
    -------
    (names, matrix):
        ``names`` in input order and the symmetric ``len(names) x len(names)``
        Jaccard matrix (diagonal 1 for non-empty sets, NaN for empty).
    """
    names = list(hh_sets)
    n = len(names)
    matrix = np.full((n, n), np.nan)
    for i, ni in enumerate(names):
        for j, nj in enumerate(names):
            matrix[i, j] = jaccard(hh_sets[ni], hh_sets[nj])
    return names, matrix
