# Example: cavitation_venturi_analysis

**Source**: `crates/cfd-1d/examples/cavitation_venturi_analysis.rs`  
**Crate**: `cfd-1d`

## Overview

Full pipeline for Sonodynamic Therapy (SDT) device design — generates a bifurcation geometry with frustum channels, solves 1D blood flow, and computes the cavitation number per channel.

## Pipeline

1. Generate bifurcation geometry via `cfd-schematics` (frustum/tapered channels)
2. Convert to 1D network and solve Hagen-Poiseuille blood flow
3. Compute cavitation number `σ = (p − p_v) / (0.5·ρ·v²)` per channel
4. Classify regime: no cavitation / inception / developed / supercavitation
5. Render colored schematics; export JSON report

## Run

```bash
cargo run -p cfd-1d --example cavitation_venturi_analysis
```

## Part Reference

Part VIII — 1-D Biomedical Flows
