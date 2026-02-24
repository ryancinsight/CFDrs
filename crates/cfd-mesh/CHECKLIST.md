# Checklist

- [x] **Phase 1: Diagnosis & Analysis**
    - [x] Run `examples/csg_cube_sphere.rs` and capture output/error.
    - [x] Audit `src/application/csg/bsp.rs` for correct plane classification logic.
    - [x] Audit `src/application/csg/split.rs` for triangle splitting numerical stability.
- [x] Check `src/application/csg/boolean` for correct set operation logic (Union, Intersection, Difference) and decompose flat file.

- [ ] **Phase 2: Core Boolean Logic Fixes**
    - [ ] Robustify plane/triangle intersection tests (e.g., epsilon handling).
    - [ ] Fix re-triangulation of split polygons (ensure CCW winding/normals).
    - [x] Handle "on-plane" cases correctly (coplanar surfaces) and decompose logic.

- [ ] **Phase 3: Tooling & Verification**
    - [ ] Create `inspect_stl` binary.
    - [ ] Add "watertight" check (edge manifoldness).
    - [ ] Add normal orientation check (outward facing).
    - [ ] Add degenerate triangle check (zero area).
