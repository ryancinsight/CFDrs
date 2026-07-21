# Example: tpms_blood_1d

**Source**: `crates/cfd-1d/examples/tpms_blood_1d.rs`  
**Crate**: `cfd-1d`

## Overview

Models blood flow through a millifluidic network whose topology derives from a TPMS (Triply Periodic Minimal Surface) Gyroid lattice. The Gyroid's bicontinuous pore network is approximated as a 1D lumped-parameter system.

## Physics

- Hagen-Poiseuille resistance model
- Carreau-Yasuda non-Newtonian blood rheology
- Multi-path bifurcation network: inlet → 2 parallel paths → outlet

## Run

```bash
cargo run -p cfd-1d --example tpms_blood_1d
```

## Part Reference

Part VIII — 1-D Biomedical Flows
