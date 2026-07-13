# Methodology

> Records *why* each transform, statistic, and graph choice was made — the
> alternatives considered and why this one won — per the "explain as you build"
> principle. The pipeline is implemented and validated end to end on a committed
> **synthetic** fixture (deterministic, code-defined; see
> `tests/fixtures/synthetic.py`). Claims below are marked **[verified on
> synthetic]** where a committed test backs them, and **[awaits real section]**
> where they can only be assessed on the actual lymph-node data.

## 1. Data and provenance
- Dataset, source, Space Ranger version, reference transcriptome, access date:
  see the provenance table in the README. Single section; no genome build, no
  CHROM:POS join.
- The synthetic fixture is a real hexagonal lattice (7×7) in full-resolution
  pixel coordinates with three planted genes — a High-High patch (`GENE_HH`), a
  dispersed noise gene (`GENE_DISP`), and a smooth gradient (`GENE_GRAD`) — so
  every owned seam has a case with a known right answer.

## 2. Coordinate frame & scalefactor hygiene
- `obsm["spatial"]` is full-resolution image pixel coordinates (x, y). The frame
  is recorded in `uns["spatial_frame"]` and checked by `assert_spatial_frame`,
  not assumed.
- Array (row/col) vs pixel (x/y) and y-axis orientation differ between image and
  plot conventions. `overlay.plot_hotspots` draws spots in the same hires-pixel
  frame as the H&E: with a background image, `imshow` sets the origin (y down)
  and no extra flip is applied; without one, the y-axis is inverted to keep the
  pixel convention.
- The scalefactor is applied in exactly one place — `overlay.to_hires_coords`
  (`hires_xy = spatial_xy * tissue_hires_scalef`). Mixing full-res coordinates
  with the scaled image is a silent misalignment (the tissue-image analog of a
  CRS/projection mismatch). **[verified on synthetic]** a test asserts the
  plotted scatter offsets equal the *scaled* coordinates, so plotting raw
  full-res coords would fail loudly. Overlays now render on the *real* hires H&E
  in the notebook; only the **visual** orientation-vs-histology check is pending
  (the automated test uses a dummy image array, so it verifies the scalefactor
  math, not visual registration).

## 3. The spatial weights graph (the load-bearing choice)
- One graph, built once with `sq.gr.spatial_neighbors(coord_type="grid",
  n_neighs=6)` — the Visium hex neighborhood — and shared by the global (squidpy)
  and local (esda) statistics.
- **[verified on synthetic]** on the 7×7 lattice this yields exactly the 25
  interior spots at degree 6, 0 isolates, symmetric; `assert_hex_neighbors`
  raises if the modal degree is not 6 (the wrong-graph silent-bug guard).
- The bridge into esda: `libpysal.weights.WSP(conn).to_W()` then
  `remap_ids(obs_names)`, so the libpysal `W` encodes byte-for-byte the same
  neighbor set squidpy produced, keyed by spot name. Row-standardized
  (`transform="r"`).
- _Alternatives:_ `coord_type="generic"` (Delaunay / radius) — rejected, ignores
  known hex geometry and can bridge tissue gaps.
- _Note:_ squidpy 1.8.2 deprecates `spatial_neighbors` (removed in 1.9, use
  `spatial_neighbors_grid`). The pinned version works; the call site is the one
  place to update on a squidpy bump.

## 4. Global Moran's I — owned implementation and validation
- `moran_scratch.morans_i` computes the statistic literally from the definition
  (`I = (n/S0) · Σ w_ij z_i z_j / Σ z_i²`), coded as an O(N²) reference, with
  `row_standardize` and a permutation null (`morans_i_permutation`).
- The permutation null shuffles values across a fixed graph; only the
  cross-product `z'Wz` varies (the mean, and hence `Σz²` and `S0`, are invariant
  under a shuffle). Its expected value is exactly `E[I] = -1/(n-1)`
  **[verified on synthetic]** the sample mean of 9,999 shuffles recovers it. The
  pseudo p-value uses esda's folded two-sided convention with the `+1`
  correction (floor `1/(n_perm+1)`), so the observed arrangement counts as one
  realization and p is never 0.
- **The central defensibility result [verified on synthetic]:** `morans_i`
  reproduces *both* `sq.gr.spatial_autocorr(mode="moran")` and `esda.Moran` on
  the same weights — but **only when the weights are row-standardized**. On
  `GENE_GRAD`, all three give `0.637977`; binary weights give `0.588864`. Both
  squidpy and esda row-standardize internally by default, so parity holds only
  with `transform="r"`. This is exactly the kind of silent mismatch the project
  exists to catch, and it is pinned by a committed test at both the single-gene
  and the full-ranking level.
- _Alternatives:_ analytic randomization z-scores (kept as cross-check only);
  Geary's C (measures local differences, not the global covariance).

## 5. Global → local
- Rank genes by global Moran's I first (`global_autocorr.rank_svgs`, a cheap gate
  — one number per gene), then run LISA / Gi* only on the top set
  (`top_genes`, the expensive per-spot pass). Running the local layer genome-wide
  wastes compute on genes with no signal and worsens the multiple-testing burden.
- **[verified on synthetic]** the squidpy ranking order equals the order by
  from-scratch `morans_i` on the same weights.

## 6. Local indicators — LISA and Getis-Ord Gi*
- `esda.Moran_Local` (LISA) gives the HH/LH/LL/HL quadrant typology;
  `local_autocorr.lisa_quadrants` is the owned sign/threshold map (quadrant code
  + p-value → label or `"ns"`), unit-tested independently of esda. High-High
  clusters for compartment markers are the compartment signal.
- `esda.G_Local(star=True)` (Gi*) gives a complementary hot/cold z-score. Both
  run on the same shared weights.
- **[verified on synthetic]** the planted `GENE_HH` patch returns as significant
  High-High (7/7 on raw p_sim); Gi* z is positive over 100% of patch spots.
- _Gotcha:_ on row-standardized weights with no diagonal, esda's Gi* *assumes*
  the star self-weight equals the maximum row weight (it emits a warning). This
  is a modeling default we currently accept rather than an explicit choice; worth
  pinning (`star=<value>` or `fill_diagonal`) when Gi* is hardened.

## 7. Multiple testing & isolates
- Per-spot p-values are FDR-corrected (Benjamini-Hochberg) **within a gene**
  (`multitest.fdr_within_gene`). The BH is hand-rolled (an owned seam) and
  **[verified on synthetic]** matches `statsmodels` `fdr_bh` exactly. Across-gene
  comparison is treated cautiously and flagged (different marginals, conditional
  LISA null).
- **Isolates are handled at the p-value layer, not just via esda.** Passing
  esda's `island_weight` (threaded through the wrappers as the explicit policy)
  is *not sufficient*: **[verified on synthetic]** a no-neighbor spot still comes
  back from esda as a *significant* Low-Low cluster (`q=3`, `p_sim=0.001`) — a
  pure artifact. The decisive fix is `multitest.mask_isolates`, which NaNs
  isolate p-values *before* FDR so they are excluded from the BH denominator and
  fall through to `"ns"` rather than a fake hotspot. `find_isolates` /
  `isolate_report` count and locate them; the spot is flagged, never dropped.
- **Detection vs noise-control use different gates by design.** On the 49-spot
  synthetic fixture, BH-FDR is so stringent it collapses the planted patch to
  1/7 significant — **and this does not improve with more permutations**
  (verified at 999 / 4,999 / 9,999; the patch p-values are genuinely moderate,
  not floored by resolution). So the *detection* property (planted patch → HH) is
  evaluated on raw `p_sim`, and the *noise-control* property (dispersed gene →
  ~no clusters) on FDR-corrected p. This is a small-n artifact of the toy
  fixture; on a real section (thousands of spots, larger compartments) FDR
  behaves far less harshly relative to signal.
- _Alternatives:_ Bonferroni (too conservative); dropping isolates (dishonest —
  silently changes `n` and the tissue footprint).

## 8. Validation against tissue architecture — **[verified on real section]**
- No compartment annotation ships in the raw section download, so "recovery"
  is assessed as **concordance of High-High LISA hotspots** (`compartments.py`):
  for each marker, the significant HH spot set (isolate-masked, FDR-corrected);
  then the pairwise Jaccard overlap of those sets. High *within*-compartment and
  low *between*-compartment overlap — a block structure — is the evidence.
- **Marker panel** — from Grasso et al. 2025 (*Eur. J. Immunol.*,
  10.1002/eji.202451218; PMC12154172), a single published human-lymph-node
  marker→domain reference rather than an assembled list. Domains and markers (all
  verified present and spatially autocorrelated here): T-zone `TRBC1`/`TRAC`;
  follicle `FDCSP`/`CR2`; germinal center `BCL6`/`MYBL1`; medulla `IGHG1`/`IGHG2`;
  B–T interface `THY1`; lymphatic `LYVE1`/`PROX1`; blood vessel `VWF`/`PECAM1`.
- **Result (V1_Human_Lymph_Node, seed 0, FDR α=0.05):** the four *compact*
  compartments form a clean block structure — mean Jaccard **0.49 within vs 0.02
  between** (follicle, germinal center, T-zone, medulla), each block completely
  segregated from the others (0.00). Same-compartment markers concordant (T-zone
  TRBC1↔TRAC 0.60; medulla IGHG1↔IGHG2 0.72; follicle FDCSP↔CR2 0.44). GC markers
  correlate with the follicle block (CR2↔MYBL1 0.24) and are 0.00 with T-zone and
  medulla — GCs nested inside follicles. The **thin/linear domains do not form HH
  hotspots**: blood vessels give 0 HH spots, lymphatics are incoherent (LYVE1 132 /
  PROX1 2 spots, mutual 0.01), the B–T interface (`THY1`) yields no block — an
  honest limit of a compartment-*hotspot* method (it captures compact regions, not
  thin structures), not a failure.
- **Honest caveats.** This concordance is markers-vs-markers: same-compartment
  overlap is partly trivial (co-expression), so the load-bearing signal is
  between-compartment *segregation*. Each spot is multiple cells, so "compartment"
  is a dominant-signal notion, not a cell-level label; and it is one section,
  seed-fixed. That is why the external benchmark below matters.
- **External benchmark — [verified, non-circular].** Our GC-marker LISA hotspots
  are scored against `germ_center`, a **morphology-drawn germinal-centre
  annotation** (368 spots) from the cell2location/Chrysalis processing of this
  exact section (Zenodo 8247779). It is non-circular w.r.t. our method: spatially
  contiguous (all 368 GC spots adjoin another) and *not* a threshold of GC
  cell-type abundance (abundance predicts it at AUC 0.97 but with real overlap —
  36 non-GC spots outrank the median GC spot). Barcodes join exactly (all 3,985
  annotated spots are in our data). Result (`annotation_auc` /
  `hotspot_enrichment`): combined GC-score **AUC 0.925** (RGS13 0.92, MYBL1 0.85,
  BCL6 0.82, AICDA 0.67); the BCL6/MYBL1 HH set is **92% precise** (154 of 168
  spots inside annotated GCs; odds ratio 185) at **0.42 recall** — FDR-gating
  captures GC *cores*, so precision is high and recall moderate (AUC, being
  threshold-free, sidesteps this). This is the project's strongest, non-circular
  validation. Overlay orientation against the real H&E remains a notebook check.

## 9. Secondary read-out: neighborhood enrichment — **[verified on real section]**
- `cluster.leiden_clusters` runs Leiden on the **expression** PCA/kNN graph
  (`sc.pp.neighbors`) — a separate object from the spatial weights graph;
  conflating the two is a classic spatial-omics error. `sq.gr.nhood_enrichment`
  across those clusters then asks whether transcriptomic clusters are spatially
  adjacent (`n_jobs=1` sidesteps a Windows multiprocessing-spawn issue).
- **Result (V1_Human_Lymph_Node, res=1.0, seed 0):** 11 Leiden clusters. Annotated
  by mean marker expression, the T-like cluster (max `CD3D`, also high `CCL21`) and
  the B-like cluster (max `MS4A1`, high `CXCL13`) are each strongly
  **self-enriched** (spatial-neighbourhood z ≈ +70 and +29) and **mutually
  depleted** (cross z ≈ **−11.6**, i.e. less adjacent than chance). An independent
  method — expression clustering + spatial adjacency, not marker LISA — recovers
  the same follicle-vs-paracortex segregation as §8, which is the point of running
  it as a corroborating read-out.
- **Caveat:** Visium spots are multi-cell mixtures, so cluster→compartment
  annotation is approximate — the extreme clusters annotate cleanly, the middle
  ones (mixed B/T signal) less so; the segregation of the *extremes* is the
  robust part.

## 10. Limitations
- Single section; spot-level (multi-cell) resolution; no image-feature analysis
  in v1 (the H&E is a registered backdrop only).
- The statistics run on the real section and compartment recovery is validated
  there (§8, including the external germinal-centre benchmark). What remains
  **uncalibrated / not eyeballed on real data**: QC thresholds and HVG/cluster
  granularity use defaults; and overlay orientation against the real H&E renders
  but has not been visually confirmed against the histology.
- BH-FDR on small spot counts is very stringent (§7); interpret local
  significance on the real section with the spot count in mind.
- squidpy `spatial_neighbors` is deprecated in the pinned version (§3).
