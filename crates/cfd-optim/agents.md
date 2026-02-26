# cfd-optim — Agent Reference

> **Role**: SDT millifluidic design optimiser — parametric sweep + genetic algorithm → top-k ranked candidates.
> **Depends on**: `cfd-core`, `cfd-1d`, `cfd-schematics`, `cfd-mesh`

---

## Purpose

`cfd-optim` generates and ranks millifluidic device designs for **Sonodynamic Therapy (SDT)**
and **Pediatric Leukapheresis** applications constrained to a 96-well plate footprint
(ANSI/SLAS 1-2004, 127.76 × 85.47 mm).

Two search strategies:
1. **Parametric sweep** (`SdtOptimizer`) — exhaustive grid over 14 fixed topology families + all
   parameter combinations; returns top-k ranked by multi-objective score.
2. **Genetic algorithm** (`GeneticOptimizer`) — real-coded GA (SBX + polynomial mutation) searching
   all 19 topology families including `AdaptiveTree`; returns top-k from any generation.

Public entry points:

```rust
// Parametric
let designs: Vec<RankedDesign> = SdtOptimizer::new(OptimMode::SdtTherapy, weights).top_5()?;
let designs: Vec<RankedDesign> = SdtOptimizer::new(mode, weights).top_k(10)?;

// Genetic
let designs: Vec<RankedDesign> = GeneticOptimizer::new(mode, weights)
    .with_population(80)
    .with_max_generations(150)
    .with_top_k(5)
    .run()?;
```

---

## Module Structure

```
src/
  lib.rs             pub re-exports + legacy compatibility API
  error.rs           OptimError: SolverFailure, ConstraintViolation, InvalidTopology
  constraints.rs     Hard constraints: footprint, haemolysis limit, manufacturing limits
  design.rs          DesignTopology (19 variants), DesignCandidate, parameter sweeps
  metrics.rs         SdtMetrics computation (39 fields): cavitation, HI, pressure, cell separation
  scoring.rs         OptimMode (7 variants), score_candidate(), score_description()
  optimizer.rs       SdtOptimizer, top_5(), top_k(), OptimStats, RankedDesign
  export.rs          save_top5_json(), save_comparison_svg(), save_schematic_svg()
  evo.rs             GeneticOptimizer, MillifluidicGenome (13 genes), decode_genome()
```

---

## Hard Constraints (`constraints.rs`)

| Constraint | Value | Source |
|------------|-------|--------|
| Device footprint | 127.76 × 85.47 mm | ANSI/SLAS 1-2004 (96-well plate) |
| Maximum wall shear stress (main channel) | ≤ 150 Pa | FDA haemolysis guidance |
| Min channel width | 100 µm | Typical millifluidic fabrication limits |
| Haemolysis index per pass | ≤ `HI_PASS_LIMIT` | Giersiepen 1990 |
| Total ΔP | ≤ inlet gauge pressure | Physical feasibility |

Any candidate violating a hard constraint is immediately **disqualified** — it never
reaches the scoring phase.

---

## Key Types

### `DesignTopology` (19 variants, `design.rs`)

| Index | Variant | Short | Venturi count | Notes |
|-------|---------|-------|---------------|-------|
| 0  | `SingleVenturi`                   | `SV` | 1 | Point cavitation |
| 1  | `BifurcationVenturi`              | `BV` | 2 | 2-arm tree |
| 2  | `TrifurcationVenturi`             | `TV` | 3 | 3-arm tree |
| 3  | `VenturiSerpentine`               | `VS` | 1 | Venturi → full-grid serpentine |
| 4  | `SerpentineGrid`                  | `SG` | 0 | Pure serpentine, all 36 wells |
| 5  | `CellSeparationVenturi`           | `CS` | 1 | Cancer focus + venturi SDT |
| 6  | `WbcCancerSeparationVenturi`      | `WC` | 1 | WBC+Cancer→center, RBC→wall |
| 7  | `DoubleBifurcationVenturi`        | `D2` | 4 | [Bi,Bi] → 2×2 grid |
| 8  | `TripleBifurcationVenturi`        | `D3` | 8 | [Bi,Bi,Bi] → octant grid |
| 9  | `DoubleTrifurcationVenturi`       | `T9` | 9 | [Tri,Tri] → 3×3 grid |
| 10 | `BifurcationTrifurcationVenturi`  | `B6` | 6 | [Bi,Tri] → 2×3 grid |
| 11 | `SerialDoubleVenturi`             | `S2` | 2 | 2× venturi in series |
| 12 | `BifurcationSerpentine`           | `BS` | 0 | Bifurcation → 2 serpentine arms |
| 13 | `TrifurcationSerpentine`          | `TS` | 0 | Trifurcation → 3 serpentine arms |
| 14 | `AsymmetricBifurcationSerpentine` | `AB` | 0 | Zweifach-Fung wide+narrow arms |
| 15 | `ConstrictionExpansionArray { n_cycles }` | `CE` | 0 | Wu Z. 2019 WBC margination |
| 16 | `SpiralSerpentine { n_turns }` | `SP` | 0 | Nivedita 2017 inertial spiral |
| 17 | `ParallelMicrochannelArray { n_channels }` | `PM` | 0 | Clinical-scale parallel array |
| 18 | `AdaptiveTree { levels, split_types }` | `AT` | N | GA-only variable-depth tree |

**Key methods**:
- `name()` → &'static str
- `short()` → &'static str
- `has_venturi()`, `has_serpentine()`, `has_distribution()` → bool
- `venturi_count()`, `parallel_venturi_count()`, `serpentine_arm_count()` → usize
- `nominal_well_coverage()` → f64 (fraction of 36 wells covered)

### `DesignCandidate` (`design.rs`)

```rust
pub struct DesignCandidate {
    pub id:                  String,        // e.g. "0001-SV"
    pub topology:            DesignTopology,
    pub flow_rate_m3_s:      f64,           // total inlet flow [m³/s]
    pub inlet_gauge_pa:      f64,           // gauge pressure above ATM [Pa]
    pub throat_diameter_m:   f64,           // venturi throat width [m]; 0 if no venturi
    pub inlet_diameter_m:    f64,           // venturi upstream diameter [m]
    pub throat_length_m:     f64,           // venturi throat length [m]
    pub channel_width_m:     f64,           // rectangular main channel width [m]
    pub channel_height_m:    f64,           // rectangular main channel height [m]
    pub serpentine_segments: usize,         // straight segment count in serpentine
    pub segment_length_m:    f64,           // length of each serpentine segment [m]
    pub bend_radius_m:       f64,           // radius of serpentine 180° bends [m]
    pub feed_hematocrit:     f64,           // blood HCT at inlet (0.0–0.45); default 0.45
}
```

Derived methods: `inlet_pressure_pa()`, `channel_velocity_m_s()`, `per_venturi_flow()`,
`to_blueprint()` → `NetworkBlueprint`, `to_channel_system()` → `ChannelSystem`.

### `SdtMetrics` (`metrics.rs`, 39 fields)

| Group | Fields |
|-------|--------|
| Cavitation | `cavitation_number`, `cavitation_potential`, `throat_shear_rate_inv_s`, `throat_shear_pa`, `throat_exceeds_fda` |
| Main-channel safety | `max_main_channel_shear_pa`, `fda_main_compliant` |
| Haemolysis | `hemolysis_index_per_pass` (Giersiepen–Wurzinger 1990) |
| Distribution | `flow_uniformity`, `well_coverage_fraction`, `mean_residence_time_s` |
| System | `total_pressure_drop_pa`, `total_path_length_mm`, `pressure_feasible` |
| 2-pop separation | `cell_separation_efficiency`, `cancer_center_fraction`, `rbc_peripheral_fraction` |
| 3-pop separation | `three_pop_sep_efficiency`, `wbc_center_fraction`, `rbc_peripheral_fraction_three_pop`, `wbc_equilibrium_pos`, `cancer_equilibrium_pos`, `rbc_equilibrium_pos` |
| Leukapheresis | `wbc_recovery`, `rbc_pass_fraction`, `wbc_purity`, `total_ecv_ml` |

### `OptimMode` (`scoring.rs`, 7 variants)

| Variant | Objective | Dominant Physics |
|---------|-----------|-----------------|
| `SdtCavitation` | Maximise cavitation, min HI, max coverage | Venturi + H-P |
| `UniformExposure` | Maximise flow uniformity + dwell time | Bifurcating array + H-P |
| `Combined { cavitation_weight, exposure_weight }` | Weighted blend | Both |
| `CellSeparation` | Cell separation + cavitation + HI | Inertial lift + Dean drag |
| `ThreePopSeparation` | WBC+Cancer→center, RBC→wall | 3-pop inertial |
| `SdtTherapy` | 35% sep + 30% HI + 20% cav + 15% uniformity | All |
| `PediatricLeukapheresis { patient_weight_kg }` | 40% WBC recovery + 30% RBC removal + 20% WBC purity + 10% throughput | Cell-cell interaction |

Hard constraint for `PediatricLeukapheresis`: ECV ≤ 10% patient TBV; wall shear ≤ 150 Pa.

### `SdtWeights` (`scoring.rs`)

Struct with pub f64 weight fields for each sub-objective in each mode.
`SdtWeights::default()` provides standard clinical weights.

### `RankedDesign` (`optimizer.rs`)

```rust
pub struct RankedDesign {
    pub rank:      usize,
    pub candidate: DesignCandidate,
    pub metrics:   SdtMetrics,
    pub score:     f64,    // weighted objective ∈ [0, 1]
}
```

---

## Parametric Optimisation Loop (`optimizer.rs`)

```
for topology in 14 fixed topologies {
    for (Q, P_gauge, d_throat, w_ch, n_segs) in parameter_sweep(topology) {
        let candidate = DesignCandidate { topology, ... };
        if constraints::check(&candidate).is_err() { continue; }
        let blueprint = candidate.to_blueprint();
        let metrics   = compute_metrics(&blueprint, &candidate)?;
        let score     = score_candidate(&metrics, mode, &weights);
        candidates.push(RankedDesign { rank: 0, candidate, metrics, score });
    }
}
candidates.sort_by(|a, b| b.score.partial_cmp(&a.score));
// assign ranks 1..k, return top k
```

`OptimStats` records: total evaluated, disqualified (constraint), CPU time.

---

## Genetic Algorithm (`evo.rs`)

### Genome

`MillifluidicGenome` has **13 normalised genes** (all ∈ [0, 1]):

| Gene | Parameter | Decode |
|------|-----------|--------|
| 0 | Topology index | `ALL_EVO_TOPOLOGIES[(gene × 18.0) as usize]` (0-18) |
| 1 | Flow rate | log-linear between `Q_MIN` and `Q_MAX` |
| 2 | Inlet gauge pressure | linear between `GAUGE_MIN` and `GAUGE_MAX` |
| 3 | Throat diameter | linear between `THROAT_MIN` and `THROAT_MAX` (venturi only) |
| 4 | Channel width | linear 2–6 mm |
| 5 | Serpentine segment count | `(gene × 10.0 + 2.0) as usize` (2–12) |
| 6 | Segment length fraction | 0.5–1.5 × `TREATMENT_WIDTH_MM` |
| 7 | Bend radius fraction | 0.05–0.25 × channel_width |
| 8 | n_cycles / n_turns / n_channels | topology-specific decode (see below) |
| 9–12 | AdaptiveTree per-level split types | bit `i` of `split_types` = 0 → Bi, 1 → Tri |

Gene 8 topology-specific decode:
- `ConstrictionExpansionArray`: `(gene × 18.0 + 2.0) as usize` → n_cycles ∈ [2, 20]
- `SpiralSerpentine`: `(gene × 18.0 + 2.0) as usize` → n_turns ∈ [2, 20]
- `ParallelMicrochannelArray`: `(gene × 490.0 + 10.0) as usize` → n_channels ∈ [10, 500]

### `ALL_EVO_TOPOLOGIES` (19 entries)

Indices 0–14: same as fixed parametric topologies (excluding AdaptiveTree).
Index 15: `ConstrictionExpansionArray { n_cycles: 10 }`
Index 16: `SpiralSerpentine { n_turns: 8 }`
Index 17: `ParallelMicrochannelArray { n_channels: 100 }`
Index 18: `AdaptiveTree { levels: 0, split_types: 0 }` (decoded from genes 9–12)

### GA Operators

- **Selection**: Tournament (k = 3)
- **Crossover**: Simulated Binary Crossover (SBX, η = 2)
- **Mutation**: Polynomial (η_m = 20, p_m = 1/n_genes)
- **Elitism**: Top 10% preserved each generation
- **Default**: pop = 60, generations = 120, top_k = 5

---

## Metrics Computation (`metrics.rs`)

### Haemolysis Index

Giersiepen–Wurzinger 1990: `HI = 3.62×10⁻⁷ · τ^2.416 · t^0.785`
Implemented via `cfd_core::physics::hemolysis::giersiepen_hi(tau_wall, exposure_time)`.

### Cavitation Number

`σ = (p∞ − p_vapour) / (½ρV²_throat)`
Computed from 1D venturi pressure via `VenturiSolver1D` (in `cfd-1d`).

### Uniformity Index

`U = 1 − (max(Qᵢ) − min(Qᵢ)) / mean(Qᵢ)` where Qᵢ = flow in branch i.
Computed from `cfd-1d` network solution.

### Cell Equilibrium Position (3-pop)

`x̃_eq` for each cell type via `cfd_1d::cell_separation::margination::lateral_equilibrium()`.
For new topologies with `feed_hematocrit > 0.01`, uses `enhanced_lateral_equilibrium()`
from `cfd_1d::cell_separation::cell_interaction` (includes CFL + margination enhancement Γ).

### WBC Recovery / Purity / ECV (Leukapheresis mode)

```
wbc_recovery = fraction of WBCs in center outlet
wbc_purity   = WBCs / (WBCs + RBCs) in center outlet
total_ecv_ml = n_channels × single_channel_volume_ml
```
Used only when mode = `PediatricLeukapheresis`.

---

## Export (`export.rs`)

| Function | Output |
|----------|--------|
| `save_top5_json(designs, path)` | JSON array of `RankedDesign` (all metrics + candidate) |
| `save_comparison_svg(designs, path, mode)` | Bar-chart comparison of top-k designs |
| `save_schematic_svg(candidate, path)` | 2D wave-channel schematic via `cfd-schematics` |

---

## Output Directories

All examples anchor outputs to the crate directory using `env!("CARGO_MANIFEST_DIR")`:

| Example | Output directory |
|---------|-----------------|
| `sdt_therapy` | `cfd-optim/outputs/sdt_therapy/` |
| `sdt_export` | `cfd-optim/outputs/` |
| `sdt_evo` | `cfd-optim/outputs/` |

---

## Relationship to Other Crates

| Crate | Relationship |
|-------|-------------|
| `cfd-1d` | All network flow physics via Hagen-Poiseuille + Casson blood solver |
| `cfd-schematics` | `to_blueprint()` calls composite presets; `to_channel_system()` for SVG |
| `cfd-core` | `ConstantPropertyFluid`, `Pressure`, haemolysis module, cell_interaction |
| `cfd-mesh` | Top-5 designs forwarded via `save_schematic_svg → BlueprintMeshPipeline` for 3D mesh export |

---

## Prohibited Patterns

- Do not call `cfd-2d` or `cfd-3d` from this crate (1D solver only here)
- Never relax the 150 Pa FDA haemolysis constraint — it is a regulatory hard limit
- `top_5()` must always return exactly 5 results (pad with lowest-scoring qualifiers)
- Scores must be normalised to [0, 1] before weighted summation
- `feed_hematocrit` must default to 0.45 for SDT modes; only set lower for leukapheresis
- Gene 8 decode must use topology-specific ranges (not a fixed formula for all topologies)

---

## Output Convention

All examples use `env!("CARGO_MANIFEST_DIR").join("outputs")` — never CWD-relative paths.
`outputs/` is listed in the root `.gitignore`.
