# cfd-fluidics Architecture

This crate is the dedicated home for **network generation and design**.

It separates topology design from solver implementation to enforce:

- **SoC (Separation of Concerns)**: design modeling vs generation orchestration vs solver adaptation
- **SRP (Single Responsibility Principle)**: each module has one reason to change
- **SSOT (Single Source of Truth)**: `NetworkBlueprint` is the canonical design representation

## Vertical Hierarchical Tree

```text
crates/cfd-fluidics/
├── src/
│   ├── domain/
│   │   ├── model/
│   │   │   ├── blueprint.rs       # SSOT aggregate for topology definitions
│   │   │   ├── ids.rs             # Node/edge identifiers
│   │   │   └── specs.rs           # Node/channel specifications
│   │   └── rules/
│   │       └── validator.rs       # Domain invariants
│   ├── application/
│   │   ├── ports/
│   │   │   └── graph_sink.rs      # Output port abstraction
│   │   └── use_cases/
│   │       └── generate_network.rs# Orchestrates validate -> build
│   ├── infrastructure/
│   │   └── adapters/
│   │       └── cfd1d_graph_sink.rs# `cfd-1d` adapter implementation
│   └── interface/
│       ├── facade/
│       │   └── fluidic_designer.rs# User-facing entrypoint
│       └── presets/
│           ├── bifurcation.rs
│           ├── trifurcation.rs
│           ├── venturi.rs
│           └── serpentine.rs
└── Cargo.toml
```

## Data Flow

1. Presets/facade create a `NetworkBlueprint`
2. `BlueprintValidator` enforces domain constraints
3. `NetworkGenerationService` invokes a `GraphSink`
4. `Cfd1dGraphSink` converts blueprint to `cfd_1d::NetworkGraph<f64>`

## Why this split

- `cfd-1d` remains focused on **physics and solving**.
- `cfd-fluidics` now owns **design-time topology generation**.
- Future adapters (e.g., JSON export, mesh preprocessor, external network tools)
  can be added under `infrastructure/adapters` without touching domain logic.
