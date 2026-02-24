# Gap Audit: cfd-mesh

## 1. Architectural Violations
- **Flat Hierarchy**: The `src/` directory currently houses too many root-level modules (`csg`, `delaunay`, `geometry`, `welding`, `watertight`, etc.) that violate the strict Domain -> Application -> Infrastructure -> Presentation clean architecture layering.
- **Dependency Inversion**: Custom structures should inverted; core domains should depend on traits, not specific third-party data structures.

## 2. Inexact Implementations / Dependencies
- **[RESOLVED]** The `Cargo.toml` specifies third party dependencies (`kiddo`, `rstar`, `bvh`, `spade`) which might introduce non-exact floating point approximations.
- A mathematically pure system requires operations based on robust geometric predicates and exact arithmetic, which we must implement directly.

## 3. Missing Infrastructure
- Complete lack of Sparse Voxel Directed Acyclic Graph (SSVDAG) / Sparse Voxel Octree (SVO) implementations. The user has specified a transition toward these state-of-the-art exact structures over k-d trees (`kiddo`).
