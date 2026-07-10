"""Build the ONE spatial weights graph and bridge it into libpysal.

This module owns the single most load-bearing choice in the project: the spatial
graph *is* the spatial weights matrix, and the same matrix must feed both the
global statistic (squidpy) and the local ones (esda). Two different graphs would
silently corrupt every Moran's I, LISA, Gi*, and enrichment result.

What it does:
1. ``build_spatial_graph`` — ``sq.gr.spatial_neighbors`` with
   ``coord_type="grid"``, ``n_neighs=6``. Visium spots sit on a hexagonal grid;
   6 is the hex neighborhood. This writes ``obsp["spatial_connectivities"]`` and
   ``obsp["spatial_distances"]``.
2. ``assert_hex_neighbors`` — verify the neighbor structure instead of assuming
   it: interior spots should have 6 neighbors. The wrong graph is a silent bug,
   so we assert (CLAUDE.md).
3. ``to_libpysal_weights`` — bridge ``obsp["spatial_connectivities"]`` into a
   ``libpysal.weights.W`` so global and local share one matrix, and
   row-standardize (``transform="r"``) for analytic p-values.

Why bridge rather than build two graphs: esda/libpysal expect a ``W`` object;
squidpy produces a sparse adjacency. Converting the squidpy adjacency guarantees
byte-for-byte the same neighbor set, which is the whole point.

Alternative considered: ``coord_type="generic"`` with Delaunay or radius-based
neighbors. Rejected for Visium — it ignores the known hex geometry and can
connect across tissue gaps, giving a graph that looks plausible but is wrong.
"""

from __future__ import annotations

import warnings

import numpy as np
import squidpy as sq
from libpysal.weights import WSP

CONNECTIVITIES_KEY = "spatial_connectivities"
DISTANCES_KEY = "spatial_distances"


def build_spatial_graph(adata, *, n_neighs: int = 6, coord_type: str = "grid"):
    """Build the hex spatial-neighbors graph in place (squidpy).

    Populates ``adata.obsp["spatial_connectivities"]`` /
    ``["spatial_distances"]``. Defaults encode the Visium hex neighborhood; the
    structure they produce is *verified* by :func:`assert_hex_neighbors`, which is
    the real guard against a wrong graph (see CLAUDE.md).
    """
    if "spatial" not in adata.obsm:
        raise ValueError('adata.obsm["spatial"] is required to build the graph')

    sq.gr.spatial_neighbors(adata, coord_type=coord_type, n_neighs=n_neighs)
    return adata


def assert_hex_neighbors(adata, *, expected: int = 6) -> dict:
    """Verify the neighbor structure; return a report, raise on a wrong graph.

    The Visium hex neighborhood means the *modal* spot degree should be
    ``expected`` (6). Interior spots have 6 neighbors; boundary spots legitimately
    have fewer, so we check the mode, not every spot. A modal degree that is not
    ``expected`` signals the wrong ``coord_type``/``n_neighs`` — the silent-bug
    class — and raises.

    Isolates (degree-0 spots) are *reported, not raised on*: a real tissue
    boundary produces some, so they are a metric, not an exception (CLAUDE.md).

    Returns
    -------
    dict
        ``mode_degree``, ``n_isolates``, ``degree_hist`` (degree -> count, the
        evidence for the hex assumption), and ``n_obs``.
    """
    if CONNECTIVITIES_KEY not in adata.obsp:
        raise ValueError(
            f'no spatial graph in obsp["{CONNECTIVITIES_KEY}"]; '
            "run build_spatial_graph first"
        )

    conn = adata.obsp[CONNECTIVITIES_KEY]
    degree = np.asarray((conn > 0).sum(axis=1)).ravel()
    vals, counts = np.unique(degree, return_counts=True)
    degree_hist = {int(v): int(c) for v, c in zip(vals, counts)}
    mode_degree = int(vals[np.argmax(counts)])
    n_isolates = int(degree_hist.get(0, 0))

    report = {
        "mode_degree": mode_degree,
        "n_isolates": n_isolates,
        "degree_hist": degree_hist,
        "n_obs": int(adata.n_obs),
    }

    if mode_degree != expected:
        raise ValueError(
            f"expected modal degree {expected} (hex neighborhood) but got "
            f"{mode_degree}; degree histogram {degree_hist}. Wrong "
            "coord_type/n_neighs?"
        )
    return report


def to_libpysal_weights(adata, *, transform: str = "r"):
    """Bridge ``obsp["spatial_connectivities"]`` into a ``libpysal.weights.W``.

    Converts the squidpy sparse adjacency so the SAME neighbor set feeds the
    global (squidpy) and local (esda) statistics. The ``W`` is keyed by
    ``adata.obs_names`` so spots stay identifiable end to end.

    ``transform="r"`` row-standardizes rows to sum to 1 — required for parity with
    squidpy's and esda's global Moran's I, both of which row-standardize
    internally (verified: on this graph binary weights give a *different* I).
    Pass ``"O"`` for the untransformed adjacency (e.g. to check neighbor sets).

    Isolates are preserved as empty rows (they appear in ``W.islands``), not
    dropped; libpysal's island warning is silenced here because they are handled
    explicitly and reported via :func:`assert_hex_neighbors`.

    Returns
    -------
    libpysal.weights.W
    """
    if CONNECTIVITIES_KEY not in adata.obsp:
        raise ValueError(
            f'no spatial graph in obsp["{CONNECTIVITIES_KEY}"]; '
            "run build_spatial_graph first"
        )

    conn = adata.obsp[CONNECTIVITIES_KEY]
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")  # island warning handled explicitly
        W = WSP(conn).to_W()

    W.remap_ids(list(adata.obs_names))  # integer ids -> spot names, order preserved
    W.transform = transform
    return W
