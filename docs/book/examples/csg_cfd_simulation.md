# Example: csg_cfd_simulation

**Crate**: `cfd-suite` (workspace root)  
**Run**: `cargo run --example csg_cfd_simulation`  
**Source**: [`examples/csg_cfd_simulation.rs`](../../../examples/csg_cfd_simulation.rs)

## What This Example Demonstrates

Blueprint-native geometry construction → mesh generation → 1D venturi screening
in a single end-to-end pipeline using the `cfd-schematics` → `gaia` → `cfd-1d` path.

| Stage | API |
|---|---|
| Parametric venturi blueprint | `cfd_schematics::venturi_rect(...)` |
| Blueprint validation | `blueprint.validate()` |
| Mesh generation | `BlueprintMeshPipeline::run(&blueprint, &config)` |
| 1D flow screening | `evaluate_venturi_screening(VenturiScreeningInput {...})` |
| Cavitation risk assessment | `assess_venturi_screening(&screening)` |

## Key Code Snippet

```rust
use cfd_schematics::venturi_rect;
use cfd_mesh::application::pipeline::{BlueprintMeshPipeline, PipelineConfig};
use cfd_1d::{evaluate_venturi_screening, VenturiScreeningInput};

let blueprint = venturi_rect("demo", 2.0e-3, 0.7e-3, 1.0e-3, 2.4e-3);
blueprint.validate()?;

let mut mesh = BlueprintMeshPipeline::run(&blueprint, &PipelineConfig {
    circular_segments: 24, axial_rings: 12, ..Default::default()
})?;

let screening = evaluate_venturi_screening(VenturiScreeningInput {
    upstream_pressure_pa:  160_000.0,
    throat_velocity_m_s:   2.85,
    density_kg_m3:         1_025.0,
    vapor_pressure_pa:     3_170.0,
    ..
})?;
```

## Architecture

```
cfd-schematics (parametric blueprint)
  └─ gaia (mesh generation)
       └─ cfd-1d (screening / reduced-order model)
```

## Book Chapter

[← Geometry Construction and CFD Coupling](../geometry_and_meshing.md)

