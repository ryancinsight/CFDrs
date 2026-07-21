# Example: blood_venturi

**Source**: `crates/cfd-2d/examples/blood_venturi.rs`  
**Crate**: `cfd-2d`

## Overview

2D blood flow through a Venturi constriction using the Immersed Boundary Method and Carreau-Yasuda non-Newtonian rheology.

## Key Parameters

| Parameter | Value |
|---|---|
| Domain | 4 cm × 1 cm |
| Throat stenosis | 50% (0.5 cm) |
| Inlet velocity | Parabolic, max 0.5 m/s |
| Rheology | Carreau-Yasuda blood model |
| Solver | SIMPLEC |

## Features

- IBM for geometry enforcement without body-fitted mesh
- Non-Newtonian viscosity updated per step from local strain rate

## Run

```bash
cargo run -p cfd-2d --example blood_venturi
```

## Part Reference

Part VIII — 2-D and Schematic Examples
