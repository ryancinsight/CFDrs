# Sprint Task List - CFD Blood Flow Implementation

## Status: February 13, 2026

### Completed Validations (Previously)
- [x] 1D Bifurcation - Analytical validation (0.00% error)
- [x] 2D Poiseuille - Analytical validation (0.72% error)
- [x] Cross-package validation - Python_CFD (11/11 tests)
- [x] Serpentine physics fixes - All 6 tests passing
- [x] Blood rheology - Rust/Python 0.00% error
- [x] Placeholder audit - 1 critical fix applied

### Completed This Sprint
- [x] Fixed cfd-3d trifurcation validation imports (NonNewtonianFluid)
- [x] Fixed cfd-3d serpentine validation imports (NonNewtonianFluid)
- [x] Fixed cfd-3d venturi validation imports (NonNewtonianFluid)
- [x] Fixed fem_boundary_conditions.rs test (n_corner_nodes argument)
- [x] Added SPRINT_TASK_LIST.md for tracking
- [x] Committed and pushed changes to git

### Pending Validations
- [ ] 3D FEM solver convergence - needs block preconditioner for saddle-point systems
- [ ] fluidsim cross-package comparison
- [ ] Complete trifurcation 3D validation
- [ ] Complete bifurcation 3D validation
- [ ] Complete serpentine 3D validation

### External Validation Sources (Available)
- external/Python_CFD/ - DrZGan Python CFD notebooks
- external/pmocz_cfd/ - Princeton CFD comparison
- fluidsim (to install) - Full CFD suite

### Next Sprint Goals
1. Investigate 3D FEM matrix conditioning
2. Implement specialized preconditioner for saddle-point systems
3. Complete validation against fluidsim
4. Run full test suite
