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
   neighbors. A no-neighbor spot has an undefined spatial lag; libpysal will warn
   and, left alone, can poison the statistic or get dropped silently. We handle
   them explicitly — assign an ``island_weight``, count and report them as a
   metric, and mark their local results as undefined rather than deleting the
   spots. Degenerate spots are a number to report, not an exception to swallow.

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

    Returns
    -------
    (rejected, qvals): boolean mask of significant spots and adjusted p-values.
    """
    raise NotImplementedError


def find_isolates(W) -> np.ndarray:
    """Indices of spots with no neighbors in ``W`` (islands).

    A reporting primitive: the count and location of isolates is a data-quality
    metric surfaced in the results, not a filter applied behind the scenes.
    """
    raise NotImplementedError


def isolate_report(W) -> dict:
    """Summarize isolates: count, fraction, and their spot indices."""
    raise NotImplementedError


def apply_island_weight(W, *, island_weight: str | float = "island_weight"):
    """Set the self/neighbor weight libpysal assigns to islands, explicitly.

    Makes the isolate policy a visible choice in the weights object rather than a
    library default that silently decides the answer for boundary spots.
    """
    raise NotImplementedError
