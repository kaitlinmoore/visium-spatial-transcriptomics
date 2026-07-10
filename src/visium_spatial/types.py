"""Shared type aliases for the package.

``SpatialWeights`` is the contract every weights-consuming function speaks. The
union exists because spatial weights reach us two ways: dense ``np.ndarray`` from
hand-built synthetic fixtures (tests, moran_scratch validation), and sparse
(``spmatrix``/``sparray``) from squidpy's ``obsp["spatial_connectivities"]`` on a
real section. Centralizing the alias keeps that dual-source reality in one place
rather than re-declared at every call site.
"""

import numpy as np
import scipy.sparse as sp
from typing import TypeAlias

SpatialWeights: TypeAlias = np.ndarray | sp.spmatrix | sp.sparray