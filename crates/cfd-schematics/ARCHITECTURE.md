# cfd-schematics Architecture

This crate is the dedicated home for **network generation and design**.

It separates topology design from solver implementation to enforce:

- **SoC (Separation of Concerns)**: design modeling vs generation orchestration
- **SRP (Single Responsibility Principle)**: each module has one reason to change
- **SSOT (Single Source of Truth)**: `NetworkBlueprint` is the canonical design representation

## Vertical Hierarchical Tree

```text
crates/cfd-schematics/
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
│   │       └── petgraph_graph_sink.rs # Concrete GraphSink for reusable DiGraph output
│   └── interface/
│       └── presets/
│           ├── bifurcation.rs
│           ├── trifurcation.rs
│           ├── venturi.rs
│           └── serpentine.rs
└── Cargo.toml
```

## Data Flow

1. Presets create a `NetworkBlueprint`
2. `BlueprintValidator` enforces domain constraints
3. `NetworkGenerationService` validates and delegates graph materialization to a sink
4. `PetgraphGraphSink` materializes a reusable `DesignGraph` for downstream crates

## Why this split

- `cfd-1d` remains focused on **physics and solving**.
- `cfd-schematics` owns **design-time topology generation**.
- `cfd-schematics` now ships reusable network/design adapters so downstream crates
  can consume shared graph-generation logic directly.
