# Example: bifurcation_schematic

**Source**: `crates/cfd-2d/examples/bifurcation_schematic.rs`  
**Crate**: `cfd-2d`

## Overview

Two-phase schematic-driven bifurcation simulation: design with `cfd-schematics` then simulate with `cfd-2d`.

## Phase 1 — Design (`cfd-schematics`)

- Symmetric bifurcation topology as `NetworkBlueprint`
- Murray's law sizing: `r_d = r_p / 2^(1/3)`
- Renders schematic PNG → `outputs/bifurcation_schematic.png`

## Phase 2 — Simulation (`cfd-2d`)

- `BifurcationSolver2D` (Casson blood model, 60×40 grid)
- Maps flow rate results to `AnalysisOverlay` (Viridis colormap)
- Renders overlay PNG → `outputs/bifurcation_flow_overlay.png`

## Run

```bash
cargo run -p cfd-2d --example bifurcation_schematic
```

## Part Reference

Part VIII — 2-D and Schematic Examples
