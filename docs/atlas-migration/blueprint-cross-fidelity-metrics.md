# Blueprint Cross-Fidelity Metrics

## Contract

`cfd-3d::blueprint_integration` exposes the reference configuration and the
channel/node traces produced while comparing schematic, cfd-1d, cfd-2d, and
cfd-3d paths. Physical values remain typed at that public boundary:

- `MassDensity` and `DynamicViscosity` configure the reference fluid.
- `VolumetricFlowRate` configures flow and carries channel/node flow traces.
- `Volume` carries schematic and meshed channel volumes.
- `Pressure` carries pressure-drop, wall-shear, and nodal-pressure traces.
- `Velocity` carries reference mean velocity traces.

Percent errors, normalized pressure-drop coefficients, solver tolerances, and
mesh/count fields remain dimensionless, numerical, or structural. The
cfd-1d/cfd-2d solver APIs currently accept base scalars, so extraction occurs
only at those adapters. Millimetre values from the mesh provider are converted
once into the Aequitas `Volume` contract at the trace boundary.

## Migration

The unit-suffixed public field names were removed. Callers construct typed
values with `from_base` and extract with `into_base` only for formula assertions
or provider adapters. No forwarding fields or compatibility aliases remain.

## Verification

The `blueprint_integration` test target covers bifurcation, serpentine,
trifurcation, and Venturi cross-fidelity traces. The acceptance gate is a
source scan with no raw dimensional public fields in the migrated structs,
typed construction in every in-tree caller, and value-semantic Nextest
coverage for the full target.
