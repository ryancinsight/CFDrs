# Example: hemolysis_serpentine_analysis

**Source**: `crates/cfd-1d/examples/hemolysis_serpentine_analysis.rs`  
**Crate**: `cfd-1d`

## Overview

Full-pipeline hemolysis analysis in a serpentine millifluidic channel. Quantifies blood damage using the Giersiepen–Wurzinger power-law model.

## Hemolysis Model

```text
D = C · τ^α · t^β
```
where τ = wall shear stress, t = channel residence time (L/v̄).

## Pipeline

1. Generate serpentine geometry via `cfd-schematics`
2. Convert to 1D network; solve Carreau-Yasuda blood flow
3. Compute wall shear stress and hemolysis index per channel
4. Render colored schematics for hemolysis index and WSS
5. Export JSON results

## Run

```bash
cargo run -p cfd-1d --example hemolysis_serpentine_analysis
```

## Part Reference

Part VIII — 1-D Biomedical Flows
