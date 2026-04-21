# CFDrs Backlog

## Structural Improvements
- [x] `cfd-schematics`: Consolidate shared components and implement unified access.
- [x] `cfd-schematics`: Enforce Clean Architecture, Dependency Inversion Principle, Single Responsibility Principle, and Single Source of Truth.
- [x] `cfd-schematics`: Replace silent no-op plotters drawer methods with real rendering and output-verified tests.
- [x] `cfd-1d`: Preserve transient composition sampling error context in the time-config projection path.
- [x] `cfd-1d`: Preallocate merged timepoint schedules in the transient composition time-config path.
- [x] `cfd-1d`: Remove the silent 1% Quemada viscosity floor and fall back to Secomb when the Quemada domain is invalid.
- [x] `cfd-1d`: Remove the short-channel Reynolds floor from the Durst entrance correction and use the published low-Re formula directly.
- [x] `cfd-2d`: Reuse a persistent pressure-velocity state workspace and validate repeated SIMPLE iterations.
- [x] `cfd-2d`: Remove the transient PIMPLE outer residual snapshot allocation by reusing the corrected velocity workspace.
- [x] `cfd-2d`: Reset Rhie-Chow coefficient caches on every update and remove the dead SIMPLEC diagonal workspace.
- [x] `cfd-2d`: Wire the pressure-velocity solver boundary-condition and viscosity inputs into the reused state workspace.
- [x] `cfd-2d`: Reject mismatched initial state layouts and non-physical fluid inputs before pressure-velocity stepping.
- [x] `cfd-math`: Preserve direct sparse solver conversion failures and reject non-finite fallback output.
- [x] `cfd-optim`: Remove the sonosensitizer activation floor from SDT report metrics and enforce zero-dose propagation at zero cavitation.
- [x] `cfd-optim`: Remove synthetic 1% synergy floors from objective and pool scoring paths.
- [x] `cfd-optim`: Remove synthetic 1% floors from venturi selectivity metric aggregation.
- [x] `cfd-optim`: Make hard-constraint candidate scoring return exact zero for infeasible designs.
- [x] `cfd-optim`: Remove the synthetic non-cavitating floor from the cavitation score.
- [x] `cfd-3d`: Fall back to first-order upwind when the WENO5-Z stencil is unavailable on small grids.
- [x] `cfd-3d`: Scale PLIC plane-bisection tolerance with the cell dimensions and reuse precomputed normals for curvature axis selection.
- [x] `cfd-3d`: Remove the 1% void-fraction damage cutoff and accumulate cavitation damage from any nonzero void fraction.
- [x] `cfd-3d`: Classify VOF mixed cells by the exact 0 < α < 1 criterion and drive compression strength from the configured coefficient.
- [ ] `cfd-schematics`: Optimize memory layout (e.g. flat vectors, zero-copy mapping).

## Rigor & Correctness
- [ ] Review all numerical bounds and geometry assumptions in `cfd-schematics`.
