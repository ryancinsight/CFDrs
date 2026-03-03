# Checklist

- [x] **Phase 1: Diagnosis & Analysis**
    - [x] Run `examples/csg_cube_sphere.rs` and capture output/error.
    - [x] Audit `src/application/csg/bsp.rs` for correct plane classification logic.
    - [x] Audit `src/application/csg/split.rs` for triangle splitting numerical stability.
- [x] Check `src/application/csg/boolean` for correct set operation logic (Union, Intersection, Difference) and decompose flat file.

- [x] **Phase 2: Core Boolean Logic Fixes**
    - [ ] Robustify plane/triangle intersection tests (e.g., epsilon handling).
    - [ ] Fix re-triangulation of split polygons (ensure CCW winding/normals).
    - [x] Handle "on-plane" cases correctly (coplanar surfaces) and decompose logic.
    - [x] Coplanar cube–cylinder CSG (caps flush with walls): union, intersection, difference — all watertight with <1% volume error.
    - [x] Fix O(n²) coplanar Difference bottleneck: AABB per-fragment pre-screening → 64-seg completes in 8ms.

- [ ] **Phase 3: Tooling & Verification**
    - [x] Create `inspect_stl` binary — watertightness, volume, normals, aspect ratio, bounds.
    - [ ] Add "watertight" check (edge manifoldness). *(done inside inspect_stl)*
    - [x] Add normal orientation check (outward facing).
    - [x] Add degenerate triangle check (zero area).
    - [x] `csg_cube_cylinder_coplanar` example: union, intersection, difference with coplanar end caps (64 segments, 8ms Difference).


- [x] **Phase 4: TPMS Primitives Expansion**
    - [x] Neovius surface — `tpms/neovius.rs` + `primitives/neovius_sphere.rs` + example
    - [x] Lidinoid — `tpms/lidinoid.rs` + `primitives/lidinoid_sphere.rs` + example
    - [x] I-WP (Schoen) — `tpms/iwp.rs` + `primitives/iwp_sphere.rs` + example
    - [x] Split P — `tpms/split_p.rs` + `primitives/split_p_sphere.rs` + example
    - [x] FRD (Schoen F-RD) — `tpms/frd.rs` + `primitives/frd_sphere.rs` + example
    - [x] Fischer-Koch C(Y) — `tpms/fischer_koch_cy.rs` + `primitives/fischer_koch_cy_sphere.rs` + example
    - [x] 398 tests pass; examples compile clean

- [x] **Phase 5: Performance, Memory & Correctness Audit (Sprint 2)**
    - [x] `reconstruct.rs`: replace `Vec<Option<VertexId>>` id_map with `HashMap` — O(face_count) vs O(pool_size)
    - [x] `reconstruct.rs`: add `region_invalid_uses_untagged_path` + `large_pool_small_face_set_no_over_alloc` tests
    - [x] `classify.rs`: fix all garbled UTF-8 mojibake section separators
    - [x] `classify.rs`: add module doc with GWN theorem, `classify_fragment` decision flowchart, complexity table, references
    - [x] `classify.rs`: optimize `gwn()` — use `norm_squared()` for near-vertex guard (avoids 3 sqrt calls per face)
    - [x] `classify.rs`: add 8 unit tests: gwn interior/exterior/clamped/empty, classify inside/outside/coplanar, centroid/tri_normal helpers
    - [x] `vertex_pool.rs`: add snap-rounding architecture diagram + 27-neighbour completeness theorem
    - [x] 419 tests pass (up from 408); 2 pre-existing v_shape failures unchanged

- [x] **Phase 6: Retain Largest Component Precision Preservation (Sprint 3)**
    - [x] Audit `cfd-mesh` hierarchy, finding Delaunay tools are optimally implemented `O(N log N)`.
    - [x] Spot precision-loss bug in `IndexedMesh::retain_largest_component` fallback to default.
    - [x] Extract `empty_clone()` to `VertexPool` and `IndexedMesh` to precisely preserve config.
    - [x] Repair `retain_largest_component` to use `empty_clone()`.
    - [x] Write two new tests: `empty_clone_preserves_custom_tolerance` and `retain_largest_component_preserves_tolerance`.
    - [x] Test suite 421 passes (2 pre-existing v-shape failures persist unaffected).

