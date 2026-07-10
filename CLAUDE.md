# CLAUDE.md — visium-spatial-statistics

## What this is
Spatial-statistics analysis of a human lymph node 10x Visium section. Primary:
global Moran's I extended to local spatial autocorrelation (local Moran's I /
LISA and Getis-Ord Gi*) via esda, on squidpy's spatial weights graph, with
hotspots overlaid on the aligned H&E and validated against known immune
compartments. Secondary: neighborhood enrichment across clusters.

## Standing principles (do not violate)
- This is a learning and exploratory project. AI-assisted code is fine, but
  every non-trivial piece comes with a plain-language explanation in comments or
  docs/methodology.md.
- Explain as you build. For each transform, statistic, and graph choice, record
  the alternatives considered and why this one won.
- Honest framing. Never claim more than the code supports. Separate what I
  implemented from what squidpy/esda do, state limitations plainly.
- Plans before execution. Propose a plan and wait for approval before large
  changes, restructures, or deleting work. I decide.
- Tests where logic is testable. The Moran's I, the graph/weights
  bridge, the LISA quadrant assignment, and the FDR/isolate logic get pytest
  tests against committed synthetic fixtures. Test the seams I own, not esda or
  squidpy internals.

## Repo conventions
- Keep raw and derived data out of git. Commit code, configs, synthetic
  fixtures, and docs only.
- README carries the exact dataset ID, source URL, Space Ranger version, the 10x
  reference transcriptome, and access date. There is no genome build here.
- docs/methodology.md explains approach, choices, and evaluation.

## Spatial-coordinate & scalefactor hygiene
- No CHROM:POS join and no reference build, but the coordinate frame is a real
  silent-bug class.
- The spatial graph IS the spatial weights matrix. Visium spots sit on a
  hexagonal grid; sq.gr.spatial_neighbors must use coord_type="grid", n_neighs=6
  (the hex neighborhood). The wrong graph silently corrupts every global Moran's,
  LISA, Gi*, and enrichment result. Assert the neighbor structure; do not assume.
- Coordinate frame / orientation. obsm["spatial"] is in full-resolution image
  pixel coordinates. Array (row/col) vs pixel (x/y) and y-axis orientation differ
  between image and plot conventions. Record and check which frame is in use.
- Scalefactors. To overlay statistics on tissue_hires_image.png, multiply
  obsm["spatial"] by tissue_hires_scalef from scalefactors_json.json (lowres has
  its own factor). Mixing full-res coordinates with a scaled image is a silent
  misalignment — the tissue-image analog of a CRS/projection mismatch. Apply the
  factor explicitly in overlay.py; do not eyeball it.
- Gene-ID/annotation hygiene from Flagship A still applies if joining across
  sections or to external annotation (harmonize to stable Ensembl IDs; report the
  drop rate). Not needed for a single section.

## Environment
Run under WSL2 (Ubuntu). uv-managed env (pyproject.toml + uv.lock), Python 3.11.
Python deps: squidpy, scanpy, anndata, esda, libpysal, leidenalg, igraph,
numpy, pandas, matplotlib, pytest. All ship Linux wheels on PyPI; no conda.

## Project-specific instructions
- Graph: build once with sq.gr.spatial_neighbors (hex), and use the SAME graph
  for global (squidpy) and local (esda). Bridge squidpy's
  obsp["spatial_connectivities"] into a libpysal weights object so the global and
  local statistics share one matrix. Row-standardize (transform="r") for analytic
  p-values.
- Validate: moran_scratch.py must reproduce squidpy's and esda's global Moran's I
  on the same weights, within tolerance, before the local layer is trusted.
- Rank spatially variable genes by global Moran's I first, then run LISA/Gi* only on
  the top set.
- Multiple testing is real: many spots x many genes. Correct per-spot p-values
  (FDR) within a gene; be explicit and cautious about comparing across genes.
- Isolates are first-class, not silent drops. Spots with no neighbors (edge/tissue
  boundary) get handled explicitly (esda island_weight), reported, not dropped
  quietly. Treat unparseable/degenerate spots as a metric, not an exception.
- Image is an aligned backdrop only. Overlay hotspots on the H&E via the
  scalefactor; do NOT extract image features (that is a deferred stretch goal).
- Check whether High-High LISA clusters for known markers recover known compartments
  (e.g. B-cell follicle markers -> follicles, T-zone markers -> paracortex). Cite
  marker sources; don't assert.