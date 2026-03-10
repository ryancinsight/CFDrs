# Gap Audit: cfd-schematics

## Phase 1: Foundation (Current)
- Auditing `cfd-schematics`...
- Goal: Identify redundancy, shared component duplication, and structural integrity violations.
- Checking for placeholders, temporary workarounds, and approximations.

## Findings
### 1. Dual-Path Topology Generation (SSOT & DRY Violation)
- `interface::presets` (e.g., `bifurcation.rs`, `serpentine.rs`) directly instantiate `NetworkBlueprint`, manually inserting nodes and channels.
- `topology::presets` combined with `topology::factory::BlueprintTopologyFactory` does exactly the same thing but via declarative `BlueprintTopologySpec`.
- **Gap**: `interface::presets` must be completely removed or refactored into thin wrappers that call `topology::presets` and `BlueprintTopologyFactory`.

### 2. Leaking Physics Constants (SOC Violation)
- `BLOOD_MU` (3.5e-3) and `shah_london_resistance` / `hp_resistance` formulas are hardcoded into schematic generation (`cfd-schematics/src/interface/presets/*` and `factory.rs`).
- **Gap**: Schematic generation (topology/geometry mappings) should not calculate fluidic resistances. That is the domain of `cfd-1d`. The schematic layer should only provide geometry (`length_m`, `width_m`, `height_m`) and let downstream solvers compute resistances based on their own runtime fluid models.
