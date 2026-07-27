# Transient composition control metrics

## Contract

The public transient-composition control boundary carries physical and
dimensionless metrics as Aequitas quantities:

- event activation times and `SimulationTimeConfig` intervals use `Time`;
- requested, calculated, and returned timepoint vectors use `Vec<Time>`;
- inlet hematocrit uses `Dimensionless`;
- edge-flow overrides and composition snapshots use
  `VolumetricFlowRate`;
- pressure-boundary overrides use `Pressure`;
- the segmented transport Courant limit uses `Dimensionless`.
- mixture fraction storage, blood hematocrit accessors, weighted blend inputs,
  comparison tolerances, and node/edge concentration queries use
  `Dimensionless`.

The simulator converts quantities to its scalar representation only when it
sorts or compares timepoints, applies pressure and flow to numerical network
solvers, or evaluates the CFL substep formula. Private transport kernels may
therefore retain scalar time vectors, but no public timepoint contract does.
Mixture arithmetic extracts base scalars only for normalization and numerical
transport. Solver residuals remain scalar because their units depend on the
assembled and scaled equation.

## Rejected alternative

Keeping event and snapshot values as `T` and relying on field documentation was
rejected: a timestamp, pressure, flow rate, and hematocrit are all
interchangeable to the compiler despite requiring different network behavior.
A compatibility facade was rejected because it would preserve the unsafe raw
construction path instead of moving callers to the Aequitas contract.

## Verification

- `cargo check -p cfd-1d --tests --offline` passes after the public contract
  migration, including the typed timepoint vector boundary.
- Existing transient droplet and literature lanes remain value-semantic
  regression coverage for the shared event and snapshot boundary.
- A dedicated typed-control regression covers time, flow, pressure,
  hematocrit, and CFL value preservation.
- Composition parity, droplet parity, and literature validation retain
  value-semantic fraction, hematocrit, and concentration assertions after the
  public representation change.
