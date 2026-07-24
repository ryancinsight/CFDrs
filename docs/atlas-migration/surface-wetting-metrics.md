# Surface and wetting physical metrics

## Context

CFDrs exposed surface roughness, contact angle, surface energy, and
fluid-solid surface tension as bare scalars at public Rust boundaries. The
documentation identified metres, radians, joules per square metre, and
newtons per metre, but those units were not represented in the types. This
allowed values from incompatible physical dimensions to be assigned without a
boundary conversion.

## Decision

Use Aequitas quantities for the public surface contracts:

- `cfd-1d::SurfaceProperties` stores roughness as `Length`, contact angle as
  `Angle`, and surface energy as `EnergyPerArea`.
- `cfd-core::WettingProperties` and `FluidSolidInterface` store contact angles
  as `Angle` and surface tension as `SurfaceTension`.
- The `InterfaceProperties` adhesion default returns `EnergyPerArea`.

The resistance model converts roughness to its canonical scalar at the Darcy
formula boundary. The adhesion law converts the typed tension and angle to
base scalars for the cosine evaluation, then returns a typed energy-per-area
quantity. No scalar forwarding fields or compatibility constructors remain.

## Rejected alternative

Keeping the public fields scalar and documenting their units would preserve the
same assignment ambiguity. A CFDrs-local wrapper would duplicate Aequitas's
quantity and serde contracts, so it is not retained.

## Verification contract

The surface-property regression asserts that typed roughness, angle, and
surface energy preserve their SI values. The material-interface regression
asserts typed water-air tension and advancing/receding angles plus the
adhesion-energy result. In-tree channel, validation, and material constructors
must compile with the typed fields; residue scans must find no raw surface
field in these public contracts. Locked `cfd-core` and `cfd-1d` checks and
Nextest, plus focused validation binaries, are the acceptance gates. Any
dependency-graph or runtime-budget blocker remains recorded in `gap_audit.md`.
