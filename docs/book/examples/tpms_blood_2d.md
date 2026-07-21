# Example: tpms_blood_2d

**Source**: `crates/cfd-2d/examples/tpms_blood_2d.rs`  
**Crate**: `cfd-2d`

## Overview

2D blood flow through TPMS-inspired periodic channel constrictions. A sinusoidal wall profile derived from the Gyroid isosurface creates periodically constricting walls.

## Geometry

```text
y_wall(x) = amplitude · (1 + cos(2π·x/λ)) / 2
```

Cells inside the wall region are masked as solid; the SIMPLE solver enforces zero velocity there.

## Physics

- **Rheology**: Carreau-Yasuda blood model (non-Newtonian)
- **Solver**: Staggered-grid FVM SIMPLE with mask-based solid geometry
- **Grid**: Staggered structured 2D

## Run

```bash
cargo run -p cfd-2d --example tpms_blood_2d
```

## Part Reference

Part VIII — 2-D and Schematic Examples
