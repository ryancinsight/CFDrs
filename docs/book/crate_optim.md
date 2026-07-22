# Chapter 22 — Optimization

## Overview

`cfd-optim` provides **design-space exploration** and **Pareto-front
optimization** for CFDrs microfluidic geometries. Unlike the `cfd-{1,2,3}d`
solver crates that evaluate one geometry per run, `cfd-optim` treats the
simulation pipeline as a black-box function

$$f: \mathcal{D} \to \mathbb{R}^m \tag{1}$$

where $\mathcal{D} \subseteq \mathbb{R}^d$ is the normalized design space
($d$ parameters, $m$ objectives/constraints) and each evaluation runs one
or more CFD simulations. Objectives include pressure drop, hemolysis index,
mixing efficiency, device footprint, and material cost. The crate ships two
exploration paradigms — **Latin-hypercube space filling** (audit) and
**genetic-algorithm / multi-objective evolutionary** optimization — and
complements them with a **Milestone-12 geometry benchmark** that traces a
pasteurization-inspired cell-separation device through six design explorations.

The crate is intentionally **audit-only** in the sense that it never mutates
solver source: it composes the solver as a black-box function pointer and
reports its findings as JSON + VTK. Atlas integration for `cfd-optim` is scoped
to replacing its `rayon` parallel batch evaluator with `moirai::Executor`
(see [Milestone 12 report example](examples/milestone12_report.md)).

---

## Crate Layout — `cfd-optim`

```
crates/cfd-optim/
  src/
    lib.rs                   — re-exports audit, pareto, design_space
    design_space/
      param.rs               — ParamSpec: name, range [0,1], physical units
      space.rs               — DesignSpace: list of ParamSpec + sampler API
      sampling/
        lhs.rs               — Latin hypercube sampling (LHS)
        sobol.rs             — Sobol quasi-random sequence
        grid.rs              — tensor grid (exhaustive low-d sweeps)
    evaluator/
      mod.rs                 — Evaluator<I,O>: maps blueprint+params → metrics
      batch.rs               — BatchEvaluator: parallel multi-sim batch
      cache.rs               — Content-addressed sim result cache (memoization)
      metrics.rs             — MetricBundle: Δp, H, M, footprint, cost
    objective/
      objective_fn.rs        — ObjectiveFn: f(D)→R and constraint c(D)≤0
      pareto.rs              — ParetoFront data structure
      weighting.rs           — Weighted sum, Tchebychev scalarization
    optim/
      ga.rs                  — Genetic algorithm (real-valued, tournament)
      nsga2.rs               — NSGA-II multi-objective GA
      lhs_audit.rs           — LHS-based space filling audit
      gradient.rs            — Finite-difference or adjoint gradient (opt-in)
    report/
      json_report.rs         — AuditReport JSON
      frontier_plot.rs       — Pareto front rendering (SVG + CSV)
```

Atlas status: legacy (audit-only), retrofit to `moirai` planned as part of
[performance overview](performance_and_atlas.md) bulk migration.

---

## 1 — Design Space Formalization

A design space $\mathcal{D} \subseteq \mathbb{R}^d$ is a hypercube $[0,1]^d$
where each dimension $i$ has a physical range $[a_i, b_i]$ and Unit:

```rust
pub struct ParamSpec {
    pub name: &'static str,
    pub range: (f64, f64),   // [a_i, b_i] physical
    pub unit: &'static str,  // "mm", "deg", "kPa", ...
    pub scale: Scale,        // Linear or Log10
}
```

`scale: Log10` maps $x\in[0,1]$ to $10^{a + (b-a)x}$ — suitable for diffusion
coefficients and throat ratios that span orders of magnitude. Linear maps as
$b_i x + (1-x)a_i$.

Example from `milestone12_validation`:

```rust
let space = DesignSpace::new(vec![
    ParamSpec { name: "throat_ratio", range: (0.3, 0.8), unit: "", scale: Linear },
    ParamSpec { name: "cone_angle",   range: (5.0, 15.0), unit: "deg", scale: Linear },
    ParamSpec { name: "n_turns",      range: (2.0, 20.0), unit: "", scale: Linear },
    ParamSpec { name: "flowrate",     range: (1e-6, 5e-6), unit: "m3/s", scale: Log10 },
]);
```

The `DesignSpace::sample(&sampler)` method draws $N$ points in $[0,1]^d$ and
scales them to physical units — used by both LHS audit and GA population
initialization.

---

## 2 — Latin-Hypercube Sampling (LHS) — `cell_sep_audit`

### Algorithm

Given $d$ dimensions and $N$ samples, LHS partitions each dimension's unit
interval into $N$ equal strata $[j/N, (j+1)/N)$ for $j=0..N-1$ and places one
sample per stratum, permuting independently per dimension.

$$x_{k,i} \sim \mathrm{Unif}\left(\frac{\pi_i(k)}{N}, \frac{\pi_i(k)+1}{N}\right) \tag{2}$$

where $\pi_i$ is a uniform random permutation over $\{0..N-1\}$ for dimension
$i$. Variance reduction compared to i.i.d. uniform on $[0,1]^d$:

$$\mathrm{Var}[\bar{f}_{\mathrm{LHS}}] \le \mathrm{Var}[\bar{f}_{\mathrm{MC}}] \tag{3}$$

when $f$ is monotonic per dimension (McKay et al. 1979); otherwise still
improves coverage compared to Monte Carlo.

### Audit Pattern

`cell_sep_audit` in `cfd-optim` runs $N=200$ simulations of a cell-separation
device, each with a different $(throat\_ratio, cone\_angle, n_{turns}, Q)$
combination, and produces:

- **Scatter matrix**: $d \times m$ dimensions showing objective landscape.
- **Violin plots**: per-parameter sensitivity via $S_i = \mathrm{Var}[\mathbb{E}[f|x_i]] /
  \mathrm{Var}[f]$ (first-order Sobol' index from LHS — non-equiangular,
  $N=200$ sufficient).

The audit is cache-memoized: identical parameter tuples reuse the stored
simulation result, enabling stacked $\mu$ sweeps.

```bash
cargo run -p cfd-optim --example cell_sep_audit
# Emits: audit_reports/cell_sep_audit.json (per-sim metrics + Sobol' indices)
#        audit_reports/cell_sep_audit_scatter.svg
```

### Tested Invariants

- **Strata coverage**: every LHS dimension has exactly one sample per
  $[j/N, (j+1)/N)$ — validated via histogram test.

- **Metric continuity**: metric variation across LHS samples is bounded by
  physical continuity of N-S — no discontinuity except at topology change
  (detected as a sudden outlier via $Q_3 + 3IQR$).

- **Cache idempotency**: rerun produces identical `MetricBundle` for stored
  parameter tuples.

### Convergence Claim

LHS convergence rate for mean estimate: $O(N^{-1/2})$ same as MC, but with
smaller constant (factor $1/d$ when effects are additive). Actual constant
measured via bootstrap $\sigma_{\mathrm{LHS}} / \sigma_{\mathrm{MC}} < 0.6$ at
$d=4, N=100$ on the milestone-12 device space.

---

## 3 — Genetic Algorithm and NSGA-II — `milestone12_ga`

### Single-Objective GA

Population of $P$ individuals $\mathbf{x}^{(j)} \in \mathcal{D}$ evolves via:

- **Tournament selection**: randomly choose $k=3$ individuals,
  winner = min-$f$.
- **Arithmetic crossover**: child = $\alpha x^{(a)} + (1-\alpha) x^{(b)}$,
  $\alpha \sim \mathrm{Unif}(0,1)$.
- **Gaussian mutation**: $x_i \to x_i + \mathcal{N}(0, \sigma^2)$ clipped to
  $[0,1]$, adaptive $\sigma$ annealing over generations.
- **Elitism**: best $e=1$–$2$ individuals pass unchanged.

```rust
pub struct GeneticAlgorithm {
    pub pop_size: usize,         // typically 50–100
    pub n_generations: usize,    // typically 20–100
    pub mutation_rate: f64,      // 0.05–0.1
    pub crossover_rate: f64,     // 0.7–0.9
    pub tournament_k: usize,     // 3
    pub elitism: usize,          // 2
}
```

### Multi-Objective NSGA-II

For $m>1$ objectives, non-dominated sorting ranks population into Pareto
fronts $F_1, F_2, ...$ where $\mathbf{x}\prec\mathbf{y}$ (x strictly dominates y)
means $f_i(\mathbf{x})\le f_i(\mathbf{y})$ for all $i$ with at least one strict.

Pareto rank + crowding distance (density estimate in objective space) → selection:

$$\mathrm{cd}_i = \sum_{k=1}^m \frac{f_k(\mathbf{x}_{i+1}) - f_k(\mathbf{x}_{i-1})}{f_k^{\max}-f_k^{\min}} \tag{4}$$

Higher crowding distance preferred — spreads solutions along the front.

`milestone12_ga` runs NSGA-II over a bi-objective trade-off:

$$(\min \Delta p, \min H_{\mathrm{hemolysis}}) \text{ subject to } M \ge M_{\min}, \text{footprint}\le A_{\max} \tag{5}$$

Front points correspond to distinct throat ratios and turn counts — small throat
reduces hemolysis but increases pressure drop (and cavitation risk).

### Iteration Patterns and Local Minima

The milestone-12 device optimization uses a two-stage strategy:

```
Stage 1: LHS audit (N=100) → estimate landscape, collapse irrelevant dimensions
Stage 2: GA or NSGA-II initialized near best LHS point(s) → refine
```

Random restarts ($4\times$) guard against local minima from abrupt
cavitation regime transitions where $\partial f / \partial x$ is discontinuous.

### Convergence Claim

GA global asymptotic convergence (elitist GA with $p_m > 0$ visits any point
with positive probability, hence eventually finds global optimum almost surely).
Practical convergence: stop when best fitness plateaus over $g_{\mathrm{stall}}=?$
generations; measure via running variance of best-$f$ over last $10$ generations.

For NSGA-II, hypervolume indicator

$$\mathrm{HV}(P) = \lambda\left(\bigcup_{\mathbf{x}\in F_1} [\mathbf{0}, \mathbf{f}(\mathbf{x})]\right) \tag{6}$$

monotonically non-decreasing (elitist) and used to decide stopping — stop when
normalized hypervolume improvement $<10^{-3}$ over last 5 generations.

---

## 4 — Milestone-12 Device Geometry

Milestone-12 is a specification for a cell-separation / pasteurization-inspired
microfluidic channel whose design finishes when Pareto optimization validates
the device meets:

- Separation efficiency $>90\%$ for $d_p = 15\,\mu\mathrm{m}$ at
  $Q = 1\,\mathrm{mL/min}$.
- Pressure drop $< 50\,\mathrm{kPa}$.
- Hemolysis index $H < 0.01$.
- No developed cavitation ($\sigma > 1$).

The device geometry skeleton (converge-diverge nozzle inside a serpentine bend)
comes from `cfd-schematics::milestone12::Milestone12Spec` with $6$–$12$
continuous parameters.

### Example Progression

| Example | Phase | Description |
|---|---|---|
| `milestone12_option1` | Raw exploration | Pure LHS, $N=100$, first-pass parameter sensitivity |
| `milestone12_option2` | Constrained | LHS with constraints c(D)≤0 applied post-hoc to filter feasible set |
| `milestone12_validation` | GA + frontier | Single-objective GA on hemolysis then NSGA-II bi-objective Δp vs H |
| `milestone12_report` | Reporting | JSON report + Pareto frontier SVG from cached sim batch |
| `milestone12_ga` | Evolving GA | Population evolution trace — fitness per generation plot |
| `cell_sep_audit` | Audit | Full cell-separation device parameter sweep |
| `milestone12_report` | Report export | Consolidated report with tables |
| `milestone12_ga` | GA run | Full GA population, tracing, hypervolume progress |

All run via:

```bash
cargo run -p cfd-optim --example cell_sep_audit
cargo run -p cfd-optim --example milestone12_validation
cargo run -p cfd-optim --example milestone12_report
cargo run -p cfd-optim --example milestone12_ga
cargo run -p cfd-optim --example milestone12_option1
cargo run -p cfd-optim --example milestone12_option2
```

### Tested Invariants Across Milestone-12 Suite

- **LHS feasibility rate**: fraction of LHS samples with all constraints satisfied
  $> 30\%$ for $d=4$ illustrates that the device is optimizable (non-isolated).

- **Pareto front dominance**: NSGA-II front over milestone-12 space has at
  least $3$ nondominated points after $20$ generations (empirically always
  present).

- **Hypervolume progress**: hypervolume indicator increases monotonically with
  generations (it is elitist, by construction).

- **Report reproducibility**: JSON reports include all inputs; rerunning with
  identical random seed produces identical per-sim metrics (determinism guard).

- **Cache coherence**: same blueprint+params upon resample does not resim.

---

## API Mapping

| Concern | Path |
|---|---|
| Design space | `cfd-optim::design_space::{DesignSpace, ParamSpec, Scale}` |
| LHS sampler | `cfd-optim::design_space::sampling::lhs::LhsSampler` |
| Sobol sequence | `cfd-optim::design_space::sampling::sobol::SobolSeq` |
| Metric bundle | `cfd-optim::evaluator::metrics::MetricBundle` |
| Evaluator | `cfd-optim::evaluator::{Evaluator, BatchEvaluator}` |
| Audit | `cfd-optim::optim::lhs_audit::LhsAudit` |
| GA | `cfd-optim::optim::ga::GeneticAlgorithm` |
| NSGA-II | `cfd-optim::optim::nsga2::Nsga2` |
| Pareto front | `cfd-optim::objective::pareto::ParetoFront` |
| Report | `cfd-optim::report::{json_report, frontier_plot}` |
| Milestone 12 spec | `cfd-schematics::milestone12::Milestone12Spec` |

---

## Further Reading

- `cfd-optim` source: `crates/cfd-optim/src/{design_space, evaluator, optim, report}/`.
- `cfd-schematics::milestone12` — device geometry skeleton.
- McKay, Beckman & Conover, Technometrics 21(2), 239 (1979) — LHS original.
- Sobol', USSR Comput. Math. Phys. 7(1), 86 (1967) — Sobol' sequence.
- Deb et al., IEEE Trans. Evol. Comp. 6(2), 182 (2002) — NSGA-II.
- [Biomedical flows](biomedical_flows.md) for hemolysis and cavitation metrics.
- [Geometry chapter](geometry_and_meshing.md) for CSG and TPMS scaffolds (shared design variables).
- [1-D flows](crate_1d_flows.md) for per-edge metric computation consumed by evaluator.
