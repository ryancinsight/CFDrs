# Transient composition control metrics

## Contract

The public transient-composition control boundary carries physical and
dimensionless metrics as Aequitas quantities:

- event activation times and `SimulationTimeConfig` intervals use `Time`;
- inlet hematocrit uses `Dimensionless`;
- edge-flow overrides and composition snapshots use
  `VolumetricFlowRate`;
- pressure-boundary overrides use `Pressure`;
- the segmented transport Courant limit uses `Dimensionless`.

The simulator converts quantities to its scalar representation only when it
sorts or compares timepoints, applies pressure and flow to numerical network
solvers, or evaluates the CFL substep formula. Mixture fraction maps remain
scalar in this increment because their dimensionless storage contract is a
separate public representation change; solver residuals remain scalar because
their units depend on the assembled and scaled equation.

## Rejected alternative

Keeping event and snapshot values as `T` and relying on field documentation was
rejected: a timestamp, pressure, flow rate, and hematocrit are all
interchangeable to the compiler despite requiring different network behavior.
A compatibility facade was rejected because it would preserve the unsafe raw
construction path instead of moving callers to the Aequitas contract.

## Verification

- `cargo check -p cfd-1d --tests --offline` passes after the public contract
  migration.
- Existing transient droplet and literature lanes remain value-semantic
  regression coverage for the shared event and snapshot boundary.
- A dedicated typed-control regression covers time, flow, pressure,
  hematocrit, and CFL value preservation.
