# Example: serpentine_venturi_1d_vs_2d

**Source**: `crates/cfd-2d/examples/serpentine_venturi_1d_vs_2d.rs`  
**Crate**: `cfd-2d`

## Overview

Cross-fidelity comparison of a serpentine-Venturi device solved at two fidelity levels: 1D lumped-parameter vs 2D full N-S. Quantifies the fidelity gap.

## Device

A serpentine channel with Venturi constrictions at each U-turn (Dean-flow apex) using the `serpentine_venturi_rect` preset from `cfd-schematics`.

## Comparison

| Fidelity | Method | Code |
|---|---|---|
| 1D | Hagen-Poiseuille network | `network_from_blueprint` |
| 2D | SIMPLE N-S per channel | `Network2DSolver::solve_all` |

Prints a channel-by-channel table: flow rates, pressure drops, wall shear stress, and outlet-flow agreement percentages.

## Dean Number Relevance

```text
De = Re · √(D_h / 2R_c)
```

Secondary Dean vortices pre-focus cells toward the centreline — ideal for cavitation-enhanced CTC lysis in SDT.

## Run

```bash
cargo run -p cfd-2d --example serpentine_venturi_1d_vs_2d
```

## Part Reference

Part VIII — 2-D and Schematic Examples
