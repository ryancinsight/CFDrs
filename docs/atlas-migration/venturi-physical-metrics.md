# Venturi and selective-cavitation physical metrics

## Context

The selective-cavitation model and its Venturi consumer already computed
pressure, density, velocity, length, viscosity, radius, and surface-tension
values with SI meanings. Their public Rust contracts still exposed those
values as `f64`, so callers could exchange pressure, velocity, and geometry
without a type-level unit boundary. `cfd-optim` then converted the typed
`cfd-1d` result back to scalars in `VenturiPlacementMetrics` and
`BlueprintVenturiMetrics`.

## Decision

Use Aequitas quantities at the public physical-metric boundaries:

- `cfd-core` selective cavitation uses `Pressure`, `MassDensity`, `Length`,
  and `SurfaceTension` for input, population, and result fields.
- `cfd-1d` Venturi screening uses `Pressure`, `MassDensity`, `Velocity`,
  `Length`, and `DynamicViscosity` for public input and physical output fields.
- `cfd-optim` placement and blueprint Venturi metrics preserve those quantity
  types through the optimization boundary.
- Formula kernels convert once to the scalar representation required by their
  numerical implementation. Serialized report DTOs remain scalar display-unit
  contracts and perform the conversion at that explicit reporting boundary.

The existing field names and scalar SI wire representation remain stable
through Aequitas's canonical scalar serde implementation. No scalar forwarding
facade is retained.

## Rejected alternative

Keeping raw `f64` fields in the optimization metrics would preserve source
compatibility but would reintroduce the exact unit-loss boundary this change
removes. A local wrapper would duplicate Aequitas's quantity contract and
create a second owner.

## Verification

The acceptance checks are source residue scans for the migrated public fields,
touched-file rustfmt and diff checks, value-semantic cfd-core/cfd-1d/cfd-optim
tests, and warning-denied package gates. The current CFDrs package gates are
blocked before source compilation by the peer root dependency transition that
resolves duplicate Aequitas/Eunomia/Proteus source identities; workspace
metadata is also blocked by the missing `D:\tmp\cutile-rs\cutile\Cargo.toml`.
