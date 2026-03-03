# cfd-mesh Algorithm Audit (2026-03-03)

## Audited areas

- `application/csg/clip/polygon2d/cdt.rs`
- `application/csg/arrangement/seam.rs`
- `application/csg/broad_phase.rs`
- `application/csg/arrangement/mod.rs`

## Baseline limitations observed

1. `cdt_clip` previously used quadratic point dedup (`O(P^2)`) and dense edge-pair scans.
2. Seam-repair fallback used brute-force greedy nearest search (`O(n^2)`), which can dominate on dense boundary loops.
3. Exact-predicate decisions were strong, but some fallback paths still relied on costly global scans that are unnecessary with locality-aware indexing.

## Modern methods researched

1. CGAL 6.1 plane clip/split reimplementation reports large speedups by exploiting clipper structure and avoiding generic boolean overhead:
   - https://www.cgal.org/2025/06/06/new_clip/
   - https://www.cgal.org/2025/10/01/cgal61/
2. CGAL PMP `autorefine_triangle_soup()` + exact-kernel guidance for robust consecutive booleans:
   - https://www.cgal.org/2024/10/22/cgal601/
   - https://doc.cgal.org/5.1/Polygon_mesh_processing/index.html
3. Exact Weiler-model CSG with exact predicates/constructions (Levy 2024/2025):
   - https://arxiv.org/abs/2405.12949
   - https://www.bibsonomy.org/bibtex/9a203311a87a6f2b7fcce191d97a0e73
4. Robust manifold CSG design tradeoffs from ManifoldCAD:
   - https://github.com/elalish/manifold/wiki/Manifold-Library
   - https://manifoldcad.org/docs/html/classmanifold_1_1_manifold.html

## Enhancements implemented in this change

1. `cdt.rs`: sweep-index broad-phase for edge candidate pruning + spatial-hash welding for near-linear expected behavior.
2. `seam.rs`: replaced brute-force fallback nearest matching with grid-accelerated nearest search while preserving merge semantics.
3. Maintained deterministic ordering in candidate processing to avoid seam-regression instability.

## Formal documentation and tests added

1. `cdt.rs`: theorem on 3x3 neighborhood completeness for welding cells.
2. `cdt.rs`: theorem on sweep-query overlap completeness.
3. `cdt.rs` property tests:
   - sweep pairs match brute-force pairs.
   - weld idempotency and near-point merge behavior.
4. `seam.rs`: theorem on 27-cell completeness for bounded nearest fallback.
5. `seam.rs` tests:
   - deterministic small-case equivalence vs brute-force reference.
   - property-based equivalence vs brute-force reference on randomized point sets.

