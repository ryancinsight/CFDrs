# Checklist

- [ ] **Phase 1: Diagnosis & Analysis**
    - [ ] Run `examples/boolean_sphere_cube.rs` and capture output/error.
    - [ ] Audit `src/csg/bsp.rs` for correct plane classification logic.
    - [ ] Audit `src/csg/split.rs` for triangle splitting numerical stability.
    - [ ] Check `src/csg/boolean.rs` for correct set operation logic (Union, Intersection, Difference).

- [ ] **Phase 2: Core Boolean Logic Fixes**
    - [ ] Robustify plane/triangle intersection tests (e.g., epsilon handling).
    - [ ] Fix re-triangulation of split polygons (ensure CCW winding/normals).
    - [ ] Handle "on-plane" cases correctly (coplanar surfaces).

- [ ] **Phase 3: Tooling & Verification**
    - [ ] Create `inspect_stl` binary.
    - [ ] Add "watertight" check (edge manifoldness).
    - [ ] Add normal orientation check (outward facing).
    - [ ] Add degenerate triangle check (zero area).
