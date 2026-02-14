# cfd-1d

1D CFD simulations for microfluidic and millifluidic networks.

## Features

- Channel-based flow networks
- Electrical circuit analogy solvers
- Pressure-driven flow simulations
- Component-based microfluidic devices (pumps, valves, sensors)
- Porous membrane and organ-compartment modeling for organ-on-chip networks
- Transient inlet composition scheduling and network-wide mixture propagation
- Finite-length droplet tracking with channel occupancy spans and boundary points
- Junction droplet split/merge behavior (flow-weighted, volume-conserving)
- Integration with the `scheme` library for 2D schematic visualization

## Scheme Integration

This crate can integrate with the [scheme](https://github.com/ryancinsight/scheme) library to provide 2D schematic visualization of 1D microfluidic networks. This is similar to how electronic circuit design tools visualize circuit layouts.

### Enabling Scheme Integration

Add the feature flag to your `Cargo.toml`:

```toml
[dependencies]
cfd-1d = { version = "0.1", features = ["scheme-integration"] }
```

### Usage Example

```rust
use cfd_1d::prelude::*;

#[cfg(feature = "scheme-integration")]
fn visualize_network() -> Result<(), Box<dyn std::error::Error>> {
    // Create a bifurcation network schematic
    let system = helpers::create_bifurcation_schematic(200.0, 100.0, 2)?;
    
    // Export to PNG
    helpers::export_to_png(&system, "network_schematic.png")?;
    
    Ok(())
}
```

### Current Limitations

**Note**: As of January 2025, the `scheme` library uses unstable Rust features (`f64::midpoint`) that require a nightly Rust compiler. To use scheme integration:

1. Install nightly Rust:

   ```bash
   rustup install nightly
   ```

2. Build with nightly:

   ```bash
   cargo +nightly build --features scheme-integration
   ```

Alternatively, you can wait for the `scheme` library to be updated to use stable Rust features, or use it without the integration by designing schematics separately and importing/exporting via JSON.

## Architecture

The crate is organized into the following modules:

- `network`: Core network representation and graph algorithms
- `channel`: Channel geometry and properties
- `components`: Microfluidic components (pumps, valves, etc.)
- `junctions/branching`: Two-way/three-way branch junction physics, solver, validation
- `solver`: Flow solvers plus transient composition and droplet tracking pipelines
- `resistance`: Resistance models for different channel geometries
- `scheme_integration`: Optional integration with the scheme library

### Naming and hierarchy conventions

- **Deep vertical hierarchy by domain**: e.g. `junctions/branching/{physics,solver,validation}`.
- **SSOT (single source of truth)**: branching implementation logic is centralized in `junctions/branching`.
- **SoC/SRP**: physics equations, solver orchestration, and validation are separated at module boundaries.
- **No compatibility shims**: branching APIs are exposed through canonical `junctions/branching` terminology.

## Native Transient Pipeline

`cfd-1d` now provides a native transient microfluidics pipeline using existing network and solver abstractions:

1. Solve the network flow field using existing steady-state solvers.
2. Run `TransientCompositionSimulator` for time-varying inlet fluid/mixture schedules.
3. Run `TransientDropletSimulator` over composition states for droplet injection, occupancy, and lifecycle tracking.

For native MMFT-style transient controls, composition also supports:

- `simulate_with_flow_events(...)` for time-scheduled edge flow-rate changes.
- `simulate_with_pressure_events(...)` for time-scheduled pressure boundary changes with per-step hydraulic re-solve.

Droplet tracking can be run directly on pressure-driven transients via:

- `TransientDropletSimulator::simulate_with_pressure_events(...)`
- `TransientDropletSimulator::simulate_with_pressure_events_and_policy(...)`

Droplet tracking can also be run directly on flow-event-driven transients via:

- `TransientDropletSimulator::simulate_with_flow_events(...)`
- `TransientDropletSimulator::simulate_with_flow_events_and_policy(...)`

Literature-anchored transient validation coverage is provided in
`tests/transient_literature_validation.rs`, including:

- Laminar pressure-flow scaling checks (`Q ∝ ΔP`, Hagen-Poiseuille regime).
- Flow-weighted junction mixing checks (mass-conservative instantaneous mixing).
- Droplet advection checks against `dx = (Q/A) dt` kinematics.
- Reynolds-range guard checks to ensure transient validation scenarios stay
   within laminar microfluidic applicability bounds.

The droplet stage supports finite-length boundaries, occupancy spans, sink/trapped transitions, and flow-weighted split/merge behavior at junctions while conserving total droplet volume.

Split behavior is policy-driven via `DropletSplitPolicy`:

- `AutoFlowWeighted` (default): split only when secondary branch flow fraction and minimum child volume thresholds are satisfied.
- `AlwaysSplit`: force flow-weighted splitting across outgoing branches.
- `NeverSplit`: route to dominant outgoing branch only.

## Design Principles

This crate follows the same design principles as the parent CFD suite:

- SOLID, CUPID, GRASP principles
- Zero-copy/zero-cost abstractions where possible
- Clean architecture with clear separation of concerns
- Extensive use of Rust's type system for safety
