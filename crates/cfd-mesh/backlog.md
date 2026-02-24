# Backlog: cfd-mesh

## High Priority
- [x] Remove `kiddo`, `rstar`, `bvh`, `spade` dependencies from `Cargo.toml`.
- [x] Restructure source tree into `domain`, `application`, `infrastructure`, `presentation` according to persona constraints.
- [x] Implement mathematically pure `BVH` based on exact arithmetic predicates.
- [ ] Implement `SSVDAG / SVO` data structures for exact spatial queries.
- [ ] Verify CSG operations remain exact and watertight post-refactor.

## Medium Priority
- [ ] Consolidate duplicated mathematical utilities across `csg` and `geometry` modules.
- [ ] Increase Proptest coverage over the new exact bounds structures.
