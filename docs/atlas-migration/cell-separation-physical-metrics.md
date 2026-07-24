# Cell-separation physical metrics

## Decision

The public cell-separation contracts preserve physical quantities with Aequitas
types. `CellProperties` uses `Length` for cell diameter and `MassDensity` for
cell density. Kappa-aware cascade stages use `Length` for treatment and
recovery hydraulic diameters and `Velocity` for parent inflow. The public
Zweifach–Fung functions also accept typed channel diameter. Optimization
summaries expose treatment and total branch widths as `Length`.

Dimensionless flow fractions, stiffness exponents, confinement ratios, and
population fractions remain scalar because their contracts are dimensionless.
Formula kernels extract base scalars only at the numerical boundary. In-tree
constructors and tests use typed values directly; no scalar compatibility
facade is retained.

## Verification contract

The focused cell-separation unit, integration, and validation tests preserve
the existing analytical routing, confinement, and enrichment assertions.
`cfd-1d` Nextest (728/728, 3 skipped), `cfd-optim` Nextest (137/137), and the
focused validation binary (16/16) pass. The full `cfd-validation` package gate
remains separately tracked because two pre-existing Venturi tests exceed the
committed 30-second budget; that exact residual is recorded in `gap_audit.md`.
