# cfd-1d

1D CFD simulations for microfluidic and millifluidic networks.

## Features

- Channel-based flow networks
- Electrical circuit analogy solvers
- Pressure-driven flow simulations
- Component-based microfluidic devices (pumps, valves, sensors)
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
- `solver`: Various solving algorithms (electrical analogy, Hagen-Poiseuille)
- `resistance`: Resistance models for different channel geometries
- `scheme_integration`: Optional integration with the scheme library

## Design Principles

This crate follows the same design principles as the parent CFD suite:
- SOLID, CUPID, GRASP principles
- Zero-copy/zero-cost abstractions where possible
- Clean architecture with clear separation of concerns
- Extensive use of Rust's type system for safety