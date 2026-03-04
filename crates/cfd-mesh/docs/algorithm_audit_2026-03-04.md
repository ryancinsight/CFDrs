# cfd-mesh Algorithm Audit (2026-03-04)

## Scope audited

- `src/application/csg/clip/polygon2d/cdt.rs`
- `src/application/csg/arrangement/planar.rs`
- `src/domain/topology/orientation.rs`
- `src/application/csg/arrangement/mod.rs` (integration path review)

## Baseline limitations observed

1. **PSLG shattering in CDT clipping performed repeated full-point scans**
   - `collect_points_on_segment_interior` was called per polygon edge and scanned the full merged point set each time.
   - Effective complexity: `O(E * P)` where `E` is polygon edges and `P` is merged vertices/intersections.
2. **Orientation repair used hash-set visitation state**
   - Correct, but with unnecessary hashing overhead and higher memory traffic in hot repair paths.
3. **Sparse-grid pathological spans were not bounded in indexed candidate queries**
   - Tiny cell widths and long segment AABBs can produce very large integer cell ranges.

## Modern methods reviewed (primary references)

1. **Arrangement/corefinement-based robust booleans**
   - CGAL corefinement docs: https://doc.cgal.org/latest/Polygon_mesh_processing/group__PMP__corefinement__grp.html
2. **Exact-predicate / exact-construction CSG direction**
   - Levy et al. (2024): https://arxiv.org/abs/2405.12949
3. **Robust exact orientation predicates**
   - Shewchuk robust predicates page: https://www.cs.cmu.edu/~quake/robust.html
4. **Generalized winding-number based robust mesh booleans (arrangement workflow)**
   - libigl tutorial (Boolean chapter): https://libigl.github.io/tutorial/
5. **Robust 2D polygon clipping with integer/exact strategy tradeoffs**
   - Clipper2 overview: https://angusj.com/clipper2/Docs/Overview.htm

## Enhancements implemented in this change

### 1) Indexed segment-interior collection for CDT shattering

- Added `PlanarPointGridIndex` in `arrangement/planar.rs`.
- Added `collect_points_on_segment_interior_indexed(...)`.
- Updated `cdt_clip` to build one point grid and reuse it for all edge shattering calls.

#### Theorem documented

- **AABB-filter completeness** for indexed candidate retrieval:
  - if a point passes `dist(point, segment) < tol`, it must lie in the segment AABB expanded by `tol`; therefore cell queries over that AABB cannot miss valid candidates.

#### Performance/memory impact

- Eliminates repeated global scans over `unique` for each segment.
- Reuses a single scratch candidate vector per polygon pass.
- Adds sparse fallback in `collect_aabb_candidates` to avoid pathological dense-cell iteration.

### 2) Dense visited-state orientation repair

- Updated `fix_orientation` to use `Vec<bool>` visited flags instead of `HashSet<FaceId>`.
- Preserved component-complete BFS traversal over disconnected face-adjacency components.

#### Theorem documented

- **Component-complete consistency repair**:
  - BFS parity propagation on each connected component repairs local opposite-edge-direction constraints; iterating all unvisited seeds covers all disconnected components.

#### Performance/memory impact

- Lower per-face visitation overhead in orientation repair.
- Reduced memory and hash churn in repeated watertight repair loops.

## New tests added

1. `arrangement/planar.rs`
   - `indexed_segment_collection_finds_all_collinear_points`
   - `indexed_segment_collection_matches_bruteforce` (proptest)
2. `domain/topology/orientation.rs`
   - `fix_orientation_repairs_disconnected_component`
   - `fix_orientation_restores_consistency_under_random_flips` (proptest)

## Verification executed

- `cargo test -p cfd-mesh --lib --no-default-features`
- `cargo test -p cfd-mesh --no-default-features`

Both pass in this workspace after the upgrades.

---

## Update B — BVH Build Path Audit (2026-03-04)

### Scope audited

- `src/infrastructure/spatial/bvh/build.rs`
- `src/infrastructure/spatial/bvh/mod.rs`
- `src/application/csg/broad_phase.rs` (integration path)

### Additional baseline limitations observed

1. **Per-node temporary allocation in partition step**
   - SAH recursion previously allocated fresh `left_buf` and `right_buf` vectors in every partition call.
   - This produced avoidable allocator traffic and transient memory churn in large CSG broad-phase builds.
2. **Degenerate SAH split fallback did not spatially partition primitive order**
   - When all primitives landed on one side of the SAH split, fallback used an index midpoint without centroid-based reordering.
   - This guaranteed recursion progress but could weaken spatial locality in pathological centroid distributions.

### Modern methods reviewed (primary references)

1. **Bucketed SAH + nth-element style split refinement**
   - PBRT v4 BVH chapter: https://www.pbr-book.org/4ed/Primitives_and_Intersection_Acceleration/Bounding_Volume_Hierarchies
2. **In-place BVH construction and memory-conscious build pipelines**
   - Karras (NVIDIA): https://research.nvidia.com/publication/2012-06_maximizing-parallelism-construction-bvhs-octrees-and-k-d-trees
3. **Recent practical acceleration in clipping/corefinement stacks**
   - CGAL clip/split update: https://www.cgal.org/2025/06/06/new_clip/

### Enhancements implemented in this update

1. **In-place two-pointer partition (allocation-free hot path)**
   - Replaced per-call temporary vectors with in-place swap partition in `bvh/build.rs::partition`.
   - Reduced auxiliary memory from per-recursion dynamic allocations to constant additional state.
2. **Deterministic median partition fallback for degenerate SAH splits**
   - Added `partition_median(...)` using `select_nth_unstable_by` on centroid axis values.
   - Ensures balanced recursion (size difference at most 1) when SAH split is degenerate.

### New theorem documentation added

1. **Partition soundness/completeness** (`partition`)
   - Left side strictly satisfies `< split`; right side satisfies `>= split`; output is a permutation.
2. **Balanced fallback partition** (`partition_median`)
   - Nth-selection guarantees split cardinalities `k` and `n-k` with `|k-(n-k)| <= 1`.

### New tests added

1. `infrastructure/spatial/bvh/build.rs`
   - `partition_predicate_and_permutation_hold` (proptest)
   - `median_partition_balances_degenerate_centroids`
   - `median_partition_keeps_all_indices`

### Verification executed for this update

- `cargo test -p cfd-mesh --lib --no-default-features` (pass)

---

## Update C — Centroid-Cached BVH + Query-Stream AABB Reuse (2026-03-04)

### Scope audited

- `src/infrastructure/spatial/bvh/build.rs`
- `src/infrastructure/spatial/bvh/mod.rs`
- `src/application/csg/broad_phase.rs`

### Remaining limitations identified

1. **Repeated centroid recomputation in SAH build**
   - `Aabb::center()` was repeatedly evaluated inside `sah_split`, `partition`, and `partition_median`.
   - This increased arithmetic overhead and cache pressure for deep trees.
2. **Repeated query-side triangle AABB computation in broad-phase**
   - Query loops recomputed `triangle_aabb` per face while traversing BVH.
   - This duplicates work when the same query side is traversed once per broad-phase call.

### Modern methods reviewed (primary references)

1. **PBRT v4 BVH chapter (SAH + HLBVH practical build variants)**
   - https://www.pbr-book.org/4ed/Primitives_and_Intersection_Acceleration/Bounding_Volume_Hierarchies
2. **Exact-kernel and autorefine guidance for robust Boolean pipelines (CGAL)**
   - https://doc.cgal.org/latest/Polygon_mesh_processing/index.html
   - https://www.cgal.org/2024/10/22/cgal601/
3. **Exact Weiler-model CSG direction**
   - https://arxiv.org/abs/2405.12949
4. **Robust predicates foundation**
   - https://www.cs.cmu.edu/~quake/robust.html

### Enhancements implemented

1. **Centroid cache for BVH build path**
   - Added `build_centroids(aabbs)` and passed cached centroids to:
     - `build_recursive`
     - `sah_split`
     - `partition`
     - `partition_median`
   - Removes repeated center reconstruction from AABB min/max in inner loops.
2. **Broad-phase AABB stream precompute on both sides**
   - `broad_phase_pairs` now computes `aabbs_a` and `aabbs_b` once and reuses them for BVH build/query loops.
   - Keeps pair semantics unchanged while reducing repeated geometry extraction.

### New theorem documentation added

1. **Centroid cache equivalence theorem** (`build_centroids`)
   - Cached centroid lookups and on-demand center computation are decision-equivalent for all centroid-axis comparisons.
2. **AABB precompute equivalence theorem** (`broad_phase_pairs`)
   - Querying with precomputed per-face AABBs is equivalent to on-demand AABB construction at the same face index.

### New tests added

1. `infrastructure/spatial/bvh/build.rs`
   - `centroid_cache_matches_direct_centers`
2. `application/csg/broad_phase.rs`
   - `broad_phase_matches_bruteforce_aabb_overlap` (proptest)

### Verification executed for this update

- `cargo test -p cfd-mesh infrastructure::spatial::bvh::build::tests --no-default-features` (pass)
- `cargo test -p cfd-mesh application::csg::broad_phase::tests --no-default-features` (pass)

---

## Update D — Deterministic Broad-Phase Emission Without Global Pair Sort (2026-03-04)

### Scope audited

- `src/application/csg/broad_phase.rs`

### Limitation identified

1. **Global candidate-pair sort in hot path**
   - `pairs.sort_unstable_by_key(|p| (p.face_a, p.face_b))` imposed full `O(k log k)` ordering over all discovered overlap pairs `k`.
   - For high-overlap scenes, sorting cost can dominate broad-phase time.

### Modern methods reviewed (primary references)

1. **Traversal-order aware BVH output strategies**
   - PBRT v4 BVH chapter: https://www.pbr-book.org/4ed/Primitives_and_Intersection_Acceleration/Bounding_Volume_Hierarchies
2. **High-performance BVH traversal/batching practice (Embree)**
   - https://www.embree.org/
3. **Robust CSG exactness context (unchanged narrow-phase correctness model)**
   - https://arxiv.org/abs/2405.12949

### Enhancement implemented

1. **Removed global pair sort and replaced it with deterministic streaming/bucketed emission**
   - Branch `build(B), query(A)`:
     - sort per-query hit list (`hits.sort_unstable()`) and emit directly in ascending `(face_a, face_b)` order.
   - Branch `build(A), query(B)`:
     - append each `face_b` hit into `face_a`-indexed buckets while iterating `face_b` in ascending order.
     - flatten buckets in ascending `face_a` order.
   - Result: deterministic lexicographic output without global `k log k` sort.

### New theorem documentation added

1. **Deterministic bucketed emission equivalence theorem** (`broad_phase_pairs`)
   - Bucket/stream emission preserves exact candidate set and reproduces lexicographic ordering equivalent to global sort.

### New tests added

1. `application/csg/broad_phase.rs`
   - `broad_phase_output_is_lexicographically_sorted` (proptest)
   - Existing brute-force equivalence property retained:
     - `broad_phase_matches_bruteforce_aabb_overlap` (proptest)

### Verification executed for this update

- `cargo test -p cfd-mesh application::csg::broad_phase::tests --no-default-features` (pass)
- `cargo test -p cfd-mesh --lib --no-default-features` (pass)

---

## Update E — Adaptive DDA Segment-Corridor Indexing for Planar CDT Shattering (2026-03-04)

### Scope audited

- `src/application/csg/arrangement/planar.rs`
- `src/application/csg/clip/polygon2d/cdt.rs` (call-site behavior check)

### Limitation identified

1. **Area-based candidate retrieval in segment shattering path**
   - Even after point indexing, candidate enumeration depended on expanded-AABB cell sweeps.
   - For long/thin segments this remains area-driven instead of segment-length-driven, increasing unnecessary bucket touch count.

### Modern methods reviewed (primary references)

1. **Fast voxel traversal (DDA stepping)**
   - Amanatides & Woo (1987): https://www.researchgate.net/publication/2611491_A_Fast_Voxel_Traversal_Algorithm_for_Ray_Tracing
2. **Production DDA traversal patterns**
   - PBRT v4 (DDA grid stepping in rendering integrators): https://www.pbr-book.org/4ed/Volume_Scattering/Media
3. **Robust exact-predicate Boolean context (unchanged correctness model)**
   - Levy et al. (2024): https://arxiv.org/abs/2405.12949

### Enhancement implemented

1. **Adaptive segment-corridor candidate collection**
   - Added DDA cell traversal over the segment with Chebyshev neighborhood expansion radius:
     - `r = ceil(tol / cell_size)`.
   - Added adaptive dispatch guard:
     - use DDA when estimated visits are proportional to occupied bins,
     - otherwise fallback to sparse occupied-bin AABB filtering.
2. **No change to geometric acceptance rules**
   - Final point acceptance remains exact projection + distance checks.
   - Only candidate generation strategy changed.

### New theorem documentation added

1. **Adaptive corridor completeness theorem** (`collect_segment_corridor_candidates`)
   - DDA branch: corridor-neighborhood completeness.
   - Sparse fallback branch: AABB-filter completeness.
   - Dispatch over complete supersets preserves completeness.
2. **Indexed collection correctness theorem** (`collect_points_on_segment_interior_indexed`)
   - Index stage affects performance only; exact acceptance predicates preserve geometric output equivalence.

### New tests added

1. `application/csg/arrangement/planar.rs`
   - `segment_corridor_candidates_cover_exact_accepts` (proptest)
   - Existing regression properties retained:
     - `indexed_segment_collection_matches_bruteforce` (proptest)
     - `indexed_segment_collection_finds_all_collinear_points`

### Verification executed for this update

- `cargo check -p cfd-mesh --no-default-features` (pass)
- `cargo test -p cfd-mesh application::csg::arrangement::planar::tests --no-default-features` (pass)
- `cargo test -p cfd-mesh --lib --no-default-features` (pass)

---

## Update F — Dense-Index + Position-Cache Vertex Consolidation (2026-03-04)

### Scope audited

- `src/application/csg/arrangement/phase3.rs`

### Limitation identified

1. **Hash-heavy inner loop in cross-mesh consolidation**
   - Previous implementation performed repeated:
     - `pool.position(...)` lookups for candidate/probe vertices,
     - `HashMap<VertexId, usize>` index resolution per accepted pair.
   - In dense seam regions this increased hash traffic and reduced cache locality.

### Modern methods reviewed (primary references)

1. **Data-oriented / cache-friendly dense indexing**
   - Fabian Giesen (DOD overview): https://fgiesen.wordpress.com/2013/07/09/a-trip-through-the-graphics-pipeline-2011-part-1/
2. **Union-find implementation practice (flat parent arrays + path compression)**
   - Tarjan (1975) + standard DSU references:
     https://en.wikipedia.org/wiki/Disjoint-set_data_structure
3. **Robust Boolean pipeline context (unchanged exact predicates upstream)**
   - Levy et al. (2024): https://arxiv.org/abs/2405.12949

### Enhancement implemented

1. **Dense local indexing for DSU unions**
   - Replaced per-pair `VertexId -> index` hash lookups with deterministic local
     indices over `[A..., B...]`.
2. **One-time position cache**
   - Copied participating `A`/`B` positions into contiguous vectors once, then
     reused in neighborhood probes.
3. **Behavior preserved**
   - Same grid traversal order and same `B -> A` parent direction, preserving
     canonical A-root behavior.

### New theorem documentation added

1. **Position-cache equivalence theorem** (`build_cross_mesh_merge_map`)
   - Cached positions are value-identical to pool lookups, so distance
     predicates and union decisions are unchanged.
2. Existing theorem set retained:
   - **Indexed-DSU equivalence**
   - **A-root invariant**

### New tests added

1. `application/csg/arrangement/phase3.rs`
   - `dense_cache_merge_map_matches_reference` (proptest):
     - New dense/cache implementation is checked against a preserved reference
       implementation of the prior hash-heavy algorithm on randomized inputs.

### Verification executed for this update

- `cargo check -p cfd-mesh --no-default-features` (pass)
- `cargo test -p cfd-mesh application::csg::arrangement::phase3::tests --no-default-features` (pass)
- `cargo test -p cfd-mesh --lib --no-default-features` (pass)
