# Example: medical_millifluidic_screening

**Source**: `crates/cfd-1d/examples/medical_millifluidic_screening.rs`  
**Crate**: `cfd-1d`

## Overview

Comprehensive medical-grade millifluidic CFD screening in a 96-well-plate footprint (ANSI/SLAS 1-2004). Combines all medical-relevant analyses in a single device geometry.

## Pipeline

1. Generate serpentine-bifurcation geometry
2. Solve Carreau-Yasuda non-Newtonian blood flow
3. Compute and visualise: hemolysis index (Giersiepen–Wurzinger), wall shear stress, FDA shear limit violations, cavitation risk, flow rate distribution, pressure field
4. Render 6 colored overlays + 1 plain schematic
5. Export full JSON report with per-channel data

## Run

```bash
cargo run -p cfd-1d --example medical_millifluidic_screening
```

## Part Reference

Part VIII — 1-D Biomedical Flows
