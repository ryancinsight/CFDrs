# cfd-optim — Agent Reference

> **Role**: SDT millifluidic design optimiser — parametric sweep → top-5 ranked candidates.  
> **Depends on**: `cfd-core`, `cfd-1d`, `cfd-schematics`

---

## Purpose

`cfd-optim` generates and ranks millifluidic device designs for **Sonodynamic Therapy (SDT)**
applications constrained to a 96-well plate footprint.

Two optimisation objectives:
1. **`SdtCavitation`** — maximise hydrodynamic cavitation via venturi geometry  
2. **`UniformExposure`** — maximise uniform acoustic / fluid exposure across a 6×6 well array

The public entry point is:

```rust
let designs: Vec<RankedDesign> = SdtOptimizer::new(OptimMode::SdtCavitation, weights)
    .top_5()?;
```

---

## Module Structure

```
src/
  lib.rs             pub mods + SdtOptimizer re-export
  error.rs           OptimError: SolverFailure, ConstraintViolation, InvalidTopology

  constraints.rs     Hard constraints: footprint, haemolysis limit, manufacturing limits
  design.rs          DesignTopology, DesignCandidate, parameter sweeps
  metrics.rs         SdtMetrics computation: cavitation number, HI, pressure drop, uniformity
  scoring.rs         score_candidate: multi-objective weighted sum
  optimizer.rs       SdtOptimizer, top_5(), OptimStats, RankedDesign
```

---

## Hard Constraints (`constraints.rs`)

| Constraint | Value | Source |
|------------|-------|--------|
| Device footprint | 127.76 × 85.47 mm | ANSI/SLAS 1-2004 (96-well plate) |
| Maximum wall shear stress | ≤ 150 Pa | FDA haemolysis guidance |
| Min channel width | 100 µm | Typical millfluidic fabrication limits |
| Max Reynolds number | ≤ 500 | Laminar flow assumption validity |
| Min cavitation number | σ < 1.0 (SdtCavitation mode) | Below inception threshold |

Any candidate violating a hard constraint is immediately **disqualified** — it never
reaches the scoring phase.

---

## Key Types

### `DesignTopology`

```rust
pub enum DesignTopology {
    SerpentineVenturi { n_bends: u8, throat_ratio: f64 },  // SDT cavitation
    BifurcatingArray  { depth: u8, symmetry: Symmetry },   // uniform exposure
    TrifurcatingArray { depth: u8 },
    StraightParallel  { n_channels: u8 },
    MixedHierarchical { stages: Vec<TopologyStage> },
}
```

### `DesignCandidate`

```rust
pub struct DesignCandidate {
    pub topology: DesignTopology,
    pub dims:     ChannelDimensions,        // width, height, throat diameter
    pub fluid:    FluidCondition,           // flow rate, blood model
}
```

### `SdtMetrics`

```rust
pub struct SdtMetrics {
    pub cavitation_number: f64,   // σ = (p − p_v) / (½ρV²)
    pub haemolysis_index:  f64,   // HI: Giersiepen–Wurzinger
    pub max_wall_shear:    f64,   // Pa — must be < 150 Pa (FDA)
    pub pressure_drop:     f64,   // Pa
    pub uniformity_index:  f64,   // 0–1 (1 = perfectly uniform)
    pub exposure_time:     f64,   // s — acoustic dwell time
}
```

### `SdtWeights`

```rust
pub struct SdtWeights {
    pub cavitation:  f64,    // weight on σ minimisation
    pub uniformity:  f64,    // weight on uniformity_index
    pub haemolysis:  f64,    // weight on HI penalty (negative)
    pub pressure:    f64,    // weight on pressure drop (negative)
}
```

### `RankedDesign`

```rust
pub struct RankedDesign {
    pub rank:      usize,           // 1 = best
    pub candidate: DesignCandidate,
    pub metrics:   SdtMetrics,
    pub score:     f64,             // weighted objective value
}
```

---

## Metrics Computation (`metrics.rs`)

### Haemolysis Index

Computed via `cfd_core::physics::hemolysis::giersiepen_hi(tau_wall, exposure_time)`.
The Giersiepen–Wurzinger correlation is defined and owned by `cfd-core`; see
`crates/cfd-core/agents.md § Physics Subsystems → Haemolysis`.

### Cavitation Number

```
σ = (p∞ − p_vapour) / (½V_throat²)
```
Definition owned by `cfd-core::physics::cavitation`. Computed here from 1D venturi
pressure at the throat via `VenturiSolver1D`.

### Uniformity Index

```
U = 1 − (max(Qᵢ) − min(Qᵢ)) / mean(Qᵢ)    where Qᵢ = flow rate in branch i
```
Computed from `cfd-1d` network flow solution.

---

## Optimisation Loop (`optimizer.rs`)

```
for topology in DesignTopology::all_variants() {
    for dims in parameter_sweep(topology) {
        let candidate = DesignCandidate { topology, dims, fluid };
        if constraints::check(&candidate).is_err() { continue; }
        let network = SchematicsConverter::convert(&candidate.to_blueprint())?;
        let flow_state = NetworkSolver::solve(&network_problem)?;
        let metrics = compute_metrics(&flow_state, &candidate)?;
        let score = score_candidate(&metrics, &weights, mode);
        candidates.push(RankedDesign { rank: 0, candidate, metrics, score });
    }
}
candidates.sort_by(|a, b| b.score.partial_cmp(&a.score));
// assign ranks 1..5, return top 5
```

`OptimStats` records: total evaluated, disqualified (constraint), CPU time.

---

## Optimisation Modes

| Mode | Objective | Dominant Physics |
|------|-----------|-----------------|
| `SdtCavitation` | Minimise σ at throat + penalise HI | Venturi + H-P resistance |
| `UniformExposure` | Maximise U (uniformity) + dwell time | Bifurcating array + H-P |

---

## Relationship to Other Crates

| Crate | Relationship |
|-------|-------------|
| `cfd-1d` | All flow physics computed via Hagen-Poiseuille + Casson blood solver |
| `cfd-schematics` | `DesignCandidate::to_blueprint()` uses schematics presets |
| `cfd-core` | `ConstantPropertyFluid`, `Pressure`, `Velocity`; haemolysis module |
| `cfd-2d` / `cfd-3d` | Not used directly; top-5 designs are forwarded for high-fidelity 2D/3D verification |

---

## Key Mathematical Theorems

| Module | Theorem |
|--------|---------|
| `metrics.rs` | Uniformity: U = 1 − (max − min) / mean |
| `scoring.rs` | Pareto dominance: design A dominates B iff all metrics of A ≥ B's |

> Physics-layer theorems (Giersiepen–Wurzinger, cavitation number, Murray's law) are
> defined in their respective owner crates (`cfd-core`, `cfd-1d`) and used here by delegation.

---

## Prohibited Patterns

- Do not call `cfd-2d` or `cfd-3d` from this crate (solver calls only in 1D)
- Never relax the 150 Pa FDA haemolysis constraint — it is a regulatory hard limit
- `top_5()` must always return exactly 5 results (pad with `None` if fewer qualify)
- Scores must be normalised to [0,1] before weighted summation (no scale bias)
