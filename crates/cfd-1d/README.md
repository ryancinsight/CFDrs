# cfd-1d

1D network solvers for microfluidic and millifluidic devices, part of the
[CFDrs](https://github.com/ryancinsight/CFDrs) simulation suite.

Each channel is reduced to a lumped hydraulic resistance and the network is solved as
a graph (the electrical circuit analogy). A whole device solves in microseconds, which
makes this the right fidelity for layout exploration and flow distribution; use
`cfd-2d` when you need resolved velocity and pressure fields in a single channel.

## Modules

- **`domain::network`** — nodes, edges, and their properties; the solved graph.
- **`domain::channel`** — geometry, cross sections, surface properties, flow state.
- **`domain::components`** — pumps, valves, flow sensors, micromixers, porous
  membranes, and organ compartments for organ-on-chip networks.
- **`domain::junctions::branching`** — two- and three-way branch junction physics,
  solvers, and validators, kept in one place rather than spread across call sites.
- **`physics::resistance`** — resistance models (Darcy-Weisbach and others) selected
  per channel geometry and flow conditions.
- **`physics::{hemolysis, cell_separation, vascular}`** — blood-specific models.
- **`solver::core`** — steady-state flow solvers plus the transient composition and
  droplet pipelines.
- **`solver::analysis`** — network analysis, performance metrics, and shear-limit
  checks against blood damage thresholds.

`prelude` re-exports the commonly used subset.

## Transient pipeline

1. Solve the steady-state network flow field.
2. Run `TransientCompositionSimulator` for time-varying inlet fluid or mixture
   schedules. `simulate_with_flow_events` and `simulate_with_pressure_events` accept
   time-scheduled edge flow-rate and boundary-pressure changes, the latter
   re-solving the hydraulics each step.
3. Run `TransientDropletSimulator` over the composition states for droplet injection,
   occupancy spans, and lifecycle tracking. It also has direct
   `simulate_with_flow_events` / `simulate_with_pressure_events` entry points, each
   with a `_and_policy` variant.

Droplets carry finite length, occupy channel spans, transition to sink/trapped states,
and split or merge at junctions conserving total volume. Split behavior is selected by
`DropletSplitPolicy`: `AutoFlowWeighted` (default — split only when the secondary
branch flow fraction and minimum child volume thresholds are met), `AlwaysSplit`, or
`NeverSplit`.

## Validation

`tests/transient_literature_validation.rs` anchors the transient pipeline to known
results: laminar pressure-flow scaling (`Q ∝ ΔP`, Hagen-Poiseuille regime),
mass-conservative flow-weighted junction mixing, droplet advection against
`dx = (Q/A) dt`, and Reynolds-range guards keeping scenarios inside the laminar
applicability bounds.

## Documentation

- [CFDrs book](https://ryancinsight.github.io/CFDrs/)
- [Workspace README](https://github.com/ryancinsight/CFDrs#readme)

## License

MIT OR Apache-2.0
