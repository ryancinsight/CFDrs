# cfd-2d Validation Report

## Models and Schemes
- Momentum, pressure–velocity coupling (SIMPLE/PISO), discretizations (central/upwind/WENO/TVD), turbulence models (k–ε, k–ω SST, SA, DES, LES variants).

## Tests
- Literature validations for turbulence modules; algorithm validation for momentum schemes; time scheme tests for explicit/implicit/multistep and adaptive.
- Stability via CFL and scheme-specific constraints; manufactured solutions for FDM Poisson/diffusion/advection–diffusion.

## Invariants
- Conservation (mass/momentum/energy) where solver applies; symmetry properties in elliptic operators; limiter monotonicity.

## Units
- Consistent SI units across fields and constants; dimensional checks embedded in tests.

## Assumptions
- Grid structures (structured/unstructured) documented; boundary condition semantics per module; limiter and scheme ranges.

## References
- Standard CFD literature for turbulence and pressure–velocity coupling; module docs reference specific models.
