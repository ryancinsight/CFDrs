# cfd-mesh clipping audit (2026-03-03)

## Scope audited

- `src/application/csg/clip/polygon2d/cdt.rs`
- `src/application/csg/clip/polygon2d/sutherland_hodgman.rs`
- `src/application/csg/clip/halfspace.rs`
- `src/application/csg/coplanar/operations.rs`

## Baseline limitations observed

1. `cdt_clip` used all-pairs edge testing for subject-vs-clip edges (`O(n*m)` candidate generation).
2. `cdt_clip` used quadratic point deduplication over all collected points (`O(P^2)`).
3. Coplanar emit path sent triangles through generic fan triangulation, adding avoidable allocation/work.

## Modern methods reviewed

1. Exact arrangement/corefinement kernels (CGAL PMP corefinement): robust boolean reliability via exact predicates/constructions.
2. Generalized winding-number arrangement pipelines (libigl mesh booleans): robust classification after arrangement extraction.
3. Sweep-line segment intersection in production geometry stacks (geo crate `sweep_points_to` docs): practical broad-phase reduction for segment events.
4. Integer-robust overlay/clipping practice (`Clipper2`): practical robustness strategy for 2D polygon booleans.

References:

- CGAL Polygon Mesh Processing corefinement docs:
  https://doc.cgal.org/latest/Polygon_mesh_processing/group__PMP__corefinement__grp.html
- libigl tutorial chapter 6:
  https://libigl.github.io/tutorial/
- geo crate boolean ops internals:
  https://docs.rs/geo/latest/geo/algorithm/bool_ops/fn.sweep_points_to.html
- Clipper2 overview:
  https://angusj.com/clipper2/Docs/Overview.htm

## Enhancements implemented in this change

1. Added sweep broad-phase index for clip edges (`ClipEdgeSweepIndex2d`) and per-subject overlap query:
   - preserves deterministic original clip-edge processing order,
   - prunes non-overlapping edge pairs before exact segment intersection.
2. Replaced quadratic dedup with uniform-grid spatial welding (`SpatialHashWeld2d`):
   - expected near-linear insertion behavior,
   - lower transient memory than canonical/remap/seen quadratic pass.
3. Added coplanar emit triangle fast-path:
   - direct lift-and-emit for 3-vertex polygons,
   - bypasses generic triangulation path for the common CDT output case.

## Formal documentation and tests added

1. Theorem doc in `cdt.rs`: 3x3-neighborhood sufficiency for spatial hash welding.
2. Theorem doc in `cdt.rs`: sweep query completeness for AABB overlap discovery.
3. Property tests in `cdt.rs`:
   - sweep overlap pairs equal brute-force overlap pairs.
   - spatial hash weld idempotency and close-point merge behavior.

