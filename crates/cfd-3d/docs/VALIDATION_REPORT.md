# cfd-3d Validation Report

## Models
- VOF/Level-Set multiphase solvers; IBM for immersed boundaries; FEM fluid modules with stabilization; spectral methods (Fourier/Chebyshev) for Poisson and related.

## Tests
- Spectral basis and Poisson validations; Chebyshev tests present.
- Planned MMS for level-set/VOF advection and IBM interpolation; FEM element consistency and stabilization checks.

## Invariants
- Conservation for VOF/Level-Set mass; IBM coupling consistency; symmetry/definiteness in elliptic solves.

## Units
- SI units across modules; dimensional checks to be integrated in extended suites.

## Assumptions
- Solver-specific constraints documented in module docs; mesh/geometry compatibility for FEM.

## References
- Multiphase and IBM literature; spectral method standards.
