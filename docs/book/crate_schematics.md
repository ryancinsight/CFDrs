# Chapter 21 — 2-D and Schematic Examples

This chapter covers schematic-driven 2D CFD — the two-phase pattern of
designing a device in `cfd-schematics` and simulating it in `cfd-2d` or
`cfd-1d`.

## Two-phase pattern

```text
Phase 1 — Design (cfd-schematics)
   NetworkBlueprint → SVG schematic PNG

Phase 2 — Simulation (cfd-2d or cfd-1d)
   Blueprint → network solve → AnalysisOverlay PNG
```

| Example | Crate | Description |
|---------|-------|-------------|
| `bifurcation_schematic` | `cfd-2d` | Bifurcation with Murray's law sizing + Casson 2D solve |
| `venturi_schematic` | `cfd-2d` | Venturi topology + ISO 5167 2D pressure solve |
| `serpentine_mixing_schematic` | `cfd-2d` | Serpentine mixer schematic + overlay |
| `blood_venturi` | `cfd-2d` | IBM + Carreau-Yasuda blood, 50% stenosis |
| `tpms_blood_2d` | `cfd-2d` | TPMS sinusoidal wall constrictions + masked SIMPLE |
| `serpentine_venturi_1d_vs_2d` | `cfd-2d` | 1D vs 2D fidelity comparison |
| `schematic_demo_integration` | `cfd-1d` | `cfd-schematics` → `cfd-1d` conversion bridge |
| `geometry_integration_demo` | `cfd-1d` | Full pipeline: generate → solve → overlay → JSON |

## Run all

```bash
cargo run -p cfd-2d --example bifurcation_schematic
cargo run -p cfd-2d --example venturi_schematic
cargo run -p cfd-2d --example blood_venturi
cargo run -p cfd-1d --example schematic_demo_integration
cargo run -p cfd-1d --example geometry_integration_demo
```

