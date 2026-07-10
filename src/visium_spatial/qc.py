"""Spot- and gene-level quality control and filtering.

Standard scanpy QC, kept separate from loading and from normalization so the
thresholds are visible and auditable rather than buried in a mega-function.

- Spot QC: total counts, number of genes per spot, mitochondrial fraction.
  Low-count / high-mito spots are tissue-edge or damaged and distort the
  spatial signal.
- Gene QC: drop genes detected in too few spots (they carry no spatial
  information and only inflate the multiple-testing burden downstream).

Choice: filter genes *before* ranking by global Moran's I, so the FDR
denominator reflects genes that could plausibly show structure. Alternative
(rank everything, filter after) inflates the correction and wastes the per-spot
local pass on noise.
"""

from __future__ import annotations


def compute_qc_metrics(adata):
    """Annotate ``adata`` with per-spot and per-gene QC metrics (in place).

    Wraps ``scanpy.pp.calculate_qc_metrics`` plus a mitochondrial-fraction column.
    """
    raise NotImplementedError


def filter_spots(adata, *, min_counts: int = 500, min_genes: int = 250, max_pct_mito: float = 30.0):
    """Return a view/copy of ``adata`` with low-quality spots removed.

    Thresholds are explicit arguments so they land in the notebook record, not a
    hidden default.
    """
    raise NotImplementedError


def filter_genes(adata, *, min_spots: int = 10):
    """Drop genes detected in fewer than ``min_spots`` spots."""
    raise NotImplementedError
