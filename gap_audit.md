## Finding 2026-08-20: CFDrs scope-vs-delivery audit

Static-evidence audit (no build, no test run) of the declared scope in
`README.md`, `docs/book/SUMMARY.md`, and `docs/adr.md` against the tree at
`68e69690` (branch `codex/cfdrs-tvd-test-integration`, 35 pre-existing dirty
files from a peer session, left untouched). Every number below is a command
result, not an estimate. No gate was executed, so nothing here asserts that any
suite passes.

### Measured baseline

| Quantity | Value | How measured |
|---|---|---|
| Workspace packages | 12 (11 libs + xtask) | `cargo metadata --offline --no-deps` |
| Library source | 1154 `.rs`, 268 175 lines under `crates/*/src` | `find`/`cat`/`wc -l` |
| All Rust incl. tests/benches/examples | 1393 `.rs`, 331 108 lines | same |
| Test functions | 3129 `#[test]` | `grep -rn '#\[test\]'` |
| `proptest!` blocks | 60 | `grep -rn 'proptest!'` |
| Cargo targets | 10 lib, 1 cdylib, 2 bin, 112 test, 49 example, 18 bench | `cargo metadata` |
| `todo!` / `unimplemented!` / TODO / FIXME | 0 / 0 / 0 | `grep` over `crates xtask tests benches examples` |
| Files > 500 lines | 135 | `find … -exec wc -l` (Atlas baseline records 134; +1 drift) |
| Crate-level `#![allow(...)]` in `lib.rs` | 292 across 10 of 11 crates | `grep -h '^#!\[allow' crates/*/src/lib.rs \| wc -l` |
| `#[allow(` sites | 97 | `grep -rn '#\[allow('` (Atlas baseline 89; **+8 regression**) |
| Crates with `#![deny(missing_docs)]` | **0 of 11** (all use `warn`) | `grep -l 'deny(missing_docs)' crates/*/src/lib.rs` |
| Production `.unwrap()` | 0 | `grep -rn '\.unwrap()' crates/*/src \| grep -v test` |
| `.expect(` sites in `crates/*/src` | 1701, of which **739 read `expect("expected value")`** | `grep -rho 'expect("[^"]*"' \| sort \| uniq -c` |
| Ignored tests | 32, of which 3 admit unimplemented behaviour | `grep -rn '#\[ignore'` |
| Book chapters | 19 numbered + 2 appendices, 1610 total lines incl. a 249-line glossary | `wc -l docs/book/*.md` |
| ADRs | 1 file (`docs/adr.md`, 1601 lines, 53 table rows), 0 numbered ADR files, no index | `ls docs/adr` → absent |

### Substrate consumption — clean at the manifest layer

No workspace or member manifest declares `ndarray`, `nalgebra`, `rayon`,
`tokio`, `rustfft`, or `approx`. The two occurrences in `Cargo.lock` are
legitimate boundary transitives: `ndarray` ← `numpy` (the PyO3 NumPy bridge in
`cfd-python`) and `rayon` ← `criterion` (dev-only). The ndarray→leto migration
is complete at the dependency-graph level; the only `ndarray::`/`nalgebra::`
string in the tree is inside `xtask/src/migration_audit.rs`, which is the
scanner that enforces their absence.

Atlas providers genuinely bound (import-site counts across `crates`, `xtask`,
`benches`, `examples`, `tests`): eunomia 864, aequitas 518, leto 411,
gaia-mesh 409 (as `cfd_mesh`), leto-ops 92, athena-leto 31, athena-core 17,
moirai 17, iris 15, hephaestus-wgpu 14, apollo-fft 8, melinoe 4, hyperion 4,
themis 3, apollo-nufft 3, hermes-simd 2, coeus-optim 2, consus 3, ritk-vtk 1,
harmonia 1, proteus 1, tyche-core 1. `moirai` is real consumption — six crates
declare it and use `Adaptive`, `fold_reduce_with`, `map_collect_mut_with`,
`ParallelSlice` — not a manifest-only edge. Only 4 `std::thread` sites remain in
library source. Zero raw `wgpu::` references survive; GPU work routes through
`hephaestus-wgpu`, though 14 in-repo `.wgsl` kernels remain
(`crates/cfd-core/src/compute/gpu/kernels/*`, `crates/cfd-math/src/shaders/*`),
which `Cargo.toml:110` already labels an in-flight replacement.

### F-1 [arch][major] Ungated `#[global_allocator]` inside a library crate

`crates/cfd-validation/src/benchmarking/memory.rs:93` installs
`TrackingAllocator` as the process global allocator with no `#[cfg]` guard:

    #[global_allocator]
    static GLOBAL_ALLOCATOR: TrackingAllocator = TrackingAllocator {
        allocator: System,
        stats: MemoryStats::new(),
    };

Every binary that links `cfd-validation` — its 40 test targets, both benches,
its examples, and any downstream consumer — silently routes all allocation
through atomic counters (`unsafe impl GlobalAlloc`, lines 315-345). Two
consequences, both observable without running anything: allocation-heavy
benchmark numbers produced anywhere in that link graph are contaminated by
instrumentation the measured code did not ask for, and a downstream consumer
that declares its own `#[global_allocator]` cannot compile against
`cfd-validation` at all (`E0152`, one allocator per program). This is public
behaviour of a published crate, so removing or gating it is `[major]`.

### F-2 [verification][patch] 54 files / 10 543 lines sit outside every cargo target

The root manifest is a virtual workspace: `grep '^\[' Cargo.toml` returns
`[workspace]`, `[workspace.package]`, `[workspace.dependencies]`,
`[workspace.lints.*]` — and no `[package]`. `Cargo.toml:26` states the cause
outright: "There is no `cfd-suite` package." Target auto-discovery is
per-package, so the root `examples/` (37 files), `benches/` (11 files), and
`tests/` (6 files) are claimed by nothing. Confirmed authoritatively:
`cargo metadata --offline --no-deps` reports **0 targets** whose `src_path`
lies under those three directories, while the 49 example / 18 bench / 112 test
targets it does report all live under `crates/*/`.

`README.md:130-141` presents those same directories as the live set — "37
example programs", "`benches/` holds 11 criterion suites", and "Examples are not
run by the test gates, so `cargo build --examples` is the check that keeps them
from rotting". That check cannot reach them. Neither can `cargo bench
--workspace`, `cargo check --workspace --all-targets` (the CI step at
`.github/workflows/ci.yml:85`), or nextest. 10 543 lines of validation drivers —
including `examples/cavity_validation.rs`, `examples/pipe_flow_validation.rs`,
and `examples/turbulent_channel_flow.rs`, the three named as figure sources in
the book manifest — are unverified by construction.

### F-3 [arch][major] 20 563 lines of committed Python domain logic in `validation/`

`git ls-files validation` returns 126 tracked files: 59 `.py` (20 563 lines),
6 `.pyc` bytecode files under `__pycache__/`, 52 `.xml` and 5 `.png` run
outputs, 1 `.ps1`, 2 `.md`. The Python is not binding glue — it is a parallel
implementation of the crate's own declared scope:
`validation/analytical/poiseuille_2d.py` and
`validation/analytical/bernoulli_venturi.py` restate analytical solutions that
`crates/cfd-validation/src/analytical/{poiseuille_2d,…}.rs` already own, and
`validation/convergence_study.py` restates `crates/cfd-validation/src/convergence/`.
`xtask/src/main.rs` carries the harness for it (`setup_venv`, `install_deps`,
`install_fenics`, `validate(…, plot, quick, category)`), so the second stack is
tooled, not vestigial.

None of it runs in CI — `.github/workflows/ci.yml` has exactly six run steps
(metadata, fmt, check, clippy, two nextest invocations, doctests) plus the
figure job. The run outputs are committed under timestamped directories
(`validation/fluidsim_output/NS2D_64x64_S2pix2pi_2026-02-09_18-10-11/…`),
which is run-output segregation inverted: the outputs are in git and the
verification is not in the gate. `CHECKLIST.md:15-17` shows the `.pyc` caches
were excluded from the *sdist*; they remain tracked in the repository.

### F-4 [arch][minor] Three in-repo SIMD implementations beside a 2-call-site hermes binding

2026 lines of SIMD live in the workspace across three unrelated homes:

- `crates/cfd-core/src/compute/simd/x86.rs` (286 lines) and `aarch64.rs` (197):
  hand-written `std::arch` intrinsics behind six `pub unsafe fn`
  (`advection_avx2`, `advection_sse41`, `diffusion_avx2`, `advection_neon`,
  `diffusion_neon`, `dot_product_neon`), each `#[target_feature(enable = …)]`.
- `crates/cfd-math/src/simd/` (972 lines: `cfd.rs`, `vector.rs`,
  `vectorization.rs`, `ops/`, `tests.rs`).
- `crates/cfd-2d/src/solvers/simd_kernels.rs` (562 lines).

`hermes-simd`, the Atlas SIMD provider, is declared by exactly one crate
(`crates/cfd-math/Cargo.toml:22`) and referenced at exactly two lines —
`crates/cfd-math/src/simd/ops/mod.rs:8` and `:196`, the second being an error
conversion. For an integrator layer this is the substrate-duplication pattern:
the provider is bound as a thin alias while the workspace keeps three
first-party implementations of the capability the provider exists to own.
Related: 17 `unsafe` sites in library source carry only 2 `// SAFETY:` comments;
the six `pub unsafe fn` do carry `# Safety` rustdoc, and the two genuine
`unsafe {}` blocks (`crates/cfd-core/src/physics/fluid_dynamics/operations.rs:41,73`)
are the two that are commented, so the gap is in the intrinsic bodies rather
than at the call sites.

### F-5 [correctness][patch] 739 content-free `expect` sites, including production solver paths

`grep -rho 'expect("[^"]*"' crates/*/src | sort | uniq -c` puts
`expect("expected value")` at 739 occurrences — 43% of all 1701 `.expect(` sites
— against 154 that carry an `invariant:`-prefixed statement. The message states
nothing a reader or a panic report can use. These are not confined to tests:

    crates/cfd-2d/src/physics/momentum/solve.rs:83   self.matrix_u.as_ref().expect("expected value"),
    crates/cfd-2d/src/physics/momentum/solve.rs:84   self.rhs_u.as_ref().expect("expected value"),
    crates/cfd-2d/src/physics/momentum/solve.rs:343  MomentumComponent::U => self.matrix_u.take().expect("expected value"),

That file contains no `mod tests`. Each site is an `Option` field whose
`Some`-ness depends on an earlier assembly call — an undocumented call-order
requirement panicking with an unlabelled message. The repo's 0-production-unwrap
conformance number is therefore satisfied in form: `clippy::unwrap_used` is
denied at `Cargo.toml:150` and the sites were converted to `expect` without the
proof the panic policy asks the message to carry.

### F-6 [verification][patch] The flagship grid-convergence study asserts nothing

`crates/cfd-2d/tests/ghia_cavity_simplec_validation.rs:1088-1123` computes an
observed order per grid pair and prints it. Line 1109 reads
`let _target_order = 2.0;` — bound, underscored, never compared. The only
branch is `if *final_error < expected_error { println!("✓ …") } else {
println!("⚠ …") }`. Replacing the solver body with a constant field would still
let this test pass; it cannot fail on a convergence regression.

The MMS layer is a source-term library rather than a code-verification gate. The
Reynolds-stress "grid convergence" test says so in its own comment —
`crates/cfd-2d/tests/reynolds_stress_comprehensive_tests.rs:583-584`: "This
tests the MMS implementation itself rather than numerical convergence." The one
genuine order-of-accuracy assertion found is
`crates/cfd-2d/tests/poisson_fdm_validation.rs:443-448` (`reduction_factor > 3.0`
for an expected 4× at second order); `crates/cfd-math/src/time_stepping/imex.rs:484`
covers temporal order. So the V&V ladder has real published-benchmark anchors
(Ghia cavity, Armaly backward step in `benchmarks/step.rs`, Taylor-Green,
Poiseuille, Couette, Blasius, Womersley, Stokes, cylinder) and real analytical
oracles, but code verification via observed order is asserted at one spatial and
one temporal site, not across the solver families that declare second order.

Three ignored tests are admissions of unimplemented declared behaviour rather
than slowness:

    crates/cfd-1d/tests/component_validation.rs:185  "Implementation does not clamp parameters to physical bounds"
    crates/cfd-1d/tests/component_validation.rs:234  "Valve resistance behavior does not follow expected monotonic relationship"
    crates/cfd-1d/tests/component_validation.rs:443  "FlowSensor component API differs from expected"

and two more park a migration: `crates/cfd-math/tests/amg_coarsening_tests.rs:5,9`
"pending migration of domain-specific multigrid code to leto-ops API".

### F-7 [verification][patch] The figure SSOT gate checks names, and its provenance is dangling

`xtask/src/check_figures.rs:94-129` — the function CI runs at
`.github/workflows/ci.yml:131` — parses `figures/*.svg` references out of
`SUMMARY.md` and `README.md`, collects the `FIGURE_SPECS` names from
`prebook.rs`, and reports the two set differences. It never opens an SVG, never
compares a hash, and never regenerates anything. A figure whose content diverged
from the data it depicts passes.

`docs/book/figures/MANIFEST.json` records each figure's
`source_example` as `{"crate_name":"cfd-suite","example_name":"…"}` —
`cfd-suite` is the package `Cargo.toml:26` states does not exist, and the named
examples (`pipe_flow_validation`, `cavity_validation`, `turbulent_channel_flow`)
are root-directory files that F-2 shows cargo cannot build. The figure → data
chain is therefore unverifiable end to end: the gate does not check content, and
the recorded producer cannot be run.

### F-8 [pm-hygiene][patch] Two parallel PM systems, a non-conformant ADR store, stale status

Both trees are live and neither defers to the other:

    root:  backlog.md (4487 L)  CHECKLIST.md (4550 L)  gap_audit.md (8705 L)   last touched 2026-08-19
    docs/: docs/backlog.md (4587 L)  docs/checklist.md (5212 L)  docs/gap_audit.md (4505 L)  last touched 2026-08-20 (68e69690)
           plus docs/gap_audit_clean.md (a fourth gap file)

`README.md:146` names the root three. The peer session at HEAD writes the
`docs/` three. `git ls-files` tracks `CHECKLIST.md`; README references
`checklist.md` — identical on this case-insensitive host, a broken reference on
a case-sensitive one.

`CFDRS-VAL-RED-1` (`backlog.md:702`) is the item this audit was asked to locate.
Its status marker reads "(in progress 2026-07-31; owner=Claude atlas session
0161539d)" while its own body at `backlog.md:729` records "CLOSURE:
cfd-validation Nextest 434/434 (1 slow)", with root causes named for both
failures and all four timeouts. The work is reported complete; only the status
field is stale. It is left as found — the item is another owner's.

ADRs do not follow the governance form: there is no `docs/adr/` directory, no
numbered records, and no index. `docs/adr.md` is one 1601-line file whose header
reads `## Status: ACTIVE - Version 1.42.0-SIMD-EXCELLENCE` (the workspace is at
`0.3.0`) above a 53-row decision table. Decisions cannot be cited by number from
board items or commits.

`docs/` additionally holds 52 sprint/audit report files
(`SPRINT_1.45.0_SUMMARY.md` … `SPRINT_1.94.0_…`, `FINAL_VALIDATION_REPORT.md`,
`AUDIT_SUMMARY_2025_11_18.md`), the report-file genre in bulk. `CHANGELOG.md`
opens with 30 lines of internal vocabulary policy before its `# Changelog`
heading; the same block is duplicated verbatim at `gap_audit.md:41`.
Repository-root non-source files include `errors.json` (174 KB of raw
`cargo --message-format=json` build output), `ACCOMPLISHMENTS.md`, and
`ARCHITECTURE.md`.

Book chapter numbering has drifted from `SUMMARY.md`: `core_flows.md` and
`governing_equations.md` both open "Chapter 2"; `pressure_velocity.md` and
`turbulence_multiphase.md` both open "Chapter 3"; `crate_schematics.md` opens
"Chapter 21" where SUMMARY lists it as 13 and its own figure is labelled 13.1;
`crate_optim.md` opens "Chapter 22" against SUMMARY's 19 / figure 19.1;
`cavitation.md`, `matrix_free_operators.md`, `schematic_integration_2d.md`, and
`vascular_bifurcations.md` carry no chapter number at all. Several chapters are
stubs — `crate_schematics.md` is 8 lines (a figure plus a cross-reference),
`core_flows.md` 11 lines.

### Conformance floor: the pedantic baseline is suppressed at crate level

`Cargo.toml:133-145` sets `clippy::all` and `clippy::pedantic` to warn and denies
`unwrap_used`, `print_stdout`, `print_stderr`, `dbg_macro`, with the comment
"the floor itself never `allow`s them — that would erase the ratchet signal."
The floor is nonetheless erased downstream: `crates/cfd-validation/src/lib.rs:1`
is `#![allow(clippy::print_stdout)]` and `crates/cfd-3d/src/lib.rs:6` allows
`print_stderr`. 409 print/`dbg` sites sit in `crates/*/src` (2673 across all
targets). Only 5 `#[expect(…)]` ratchet sites exist against 97 `#[allow(`.
`README.md:113-119` describes this honestly, which is why the crate-level count
in that paragraph is worth keeping accurate — it read 288 and measures 292.

All 12 packages are `edition = "2021"` with `resolver = "2"`; the toolchain pin
is `1.97.0`. `crates/cfd-python` ships `pyproject.toml` and a README but no
`py.typed` and no `.pyi` stubs, so the Python surface is untyped to mypy and
IDEs; it also carries a committed `casson_rheology_validation.png`.

`crates/cfd-validation/src/convergence/richardson.rs` (224 lines) and
`crates/cfd-validation/src/manufactured/richardson/core.rs` (551 lines) are two
implementations of Richardson extrapolation in one crate — `estimate_order`,
`extrapolate`, and `is_asymptotic` exist in both, the second returning
`Result<_, String>` (9 stringly-typed error returns in `cfd-validation`, 56
workspace-wide).

### Completeness

66% of declared scope delivered and verified. Denominator: the capability set
named in `README.md` (11 crates and their stated solver families), the 19
chapters plus 2 appendices of `docs/book/SUMMARY.md`, and the 53 accepted
decisions in `docs/adr.md`. Weighting per the Atlas rubric — capabilities
implemented without stubs 40% (scored ≈0.92: no `todo!`/`unimplemented!`
anywhere, every declared solver family present, minus five ignored tests that
admit unfinished behaviour), verification depth 25% (≈0.55: strong analytical
and published-benchmark anchors, but observed-order code verification asserted
at two sites, the flagship convergence study vacuous, 10 543 lines of drivers
and 20 563 lines of Python validation outside every gate, figure provenance
dangling), documentation 20% (≈0.55: complete book skeleton with stub chapters
and drifted numbering, zero crates at `deny(missing_docs)`, non-conformant ADR
store), conformance floor 15% (≈0.30: 135 oversized files, 292 crate-level
allows, +8 `#[allow]` regression against the Atlas baseline, edition 2021).

## ATLAS-CFDRS-BACKWARD-STEP-108 — provider-owned geometry and shear (hosted closure pending 2026-08-19)

## Finding 2026-08-19: exact default-head Rust and figure gates pass

CFDrs default `931ee3a0130a5238461a1ee9547e12aef11e90bf` passes hosted run
`32222487306`. The Rust workspace job completes format, check, Clippy, tests,
the numerical-fidelity suite, doctests, and native fontconfig setup. The figure
SSOT job also passes. This is behavioral and figure evidence for the exact
provider head, not proof of Pages deployment or PyPI publication; those remain
separate gates. The local standalone locked package path remains blocked by the
Atlas development overlay, so no local package result is claimed.

The first implementation placed a new streamfunction–vorticity solver in
`cfd-validation`, which violated provider-first ownership because `cfd-2d`
already owns the SIMPLE Navier–Stokes field and mask path. The follow-up moves
the step geometry mask, explicit inlet/outlet/no-slip contract, signed
downstream lower-wall shear samples, and interpolated reattachment crossing to
`crates/cfd-2d/src/solvers/ns_fvm/backward_step.rs`. The validation benchmark is
now a thin adapter and no longer owns a second field solver or matrix type.
The provider applies a normalized parabolic profile only on fluid inlet cells;
solid inlet cells remain zero during every SIMPLE iteration.

The provider now uses the Armaly inlet hydraulic diameter
`D = 2 * (channel_height - step_height)` in its viscosity calculation, and the
validation result carries the configured Reynolds number and geometry. The
adapter exposes only the published two-dimensional anchors `Re=100,
x_r/h=2.84` and `Re=389, x_r/h=7.83`; it does not interpolate the nonlinear
literature curve. The exact release filter passes 14/14 at run
`314957f4-817b-44bb-bea6-0f1138f7b6c7`; the default solve takes 5.79 s and
returns `x_r/h = 2.0016`, which remains at the edge of the existing 30%
acceptance band. The debug termination and the remaining SIMPLE/grid-fidelity
error are still open. A 96-cell probe exceeded the committed slow budget and
did not close the fidelity gap, so its workload was not retained. Full
Atlas-locked compilation remains separately blocked by the shared overlay's
peer-dirty Asclepius checkout requiring `aequitas ^0.1.0` while the current
local graph provides `0.2.0`.

## ATLAS-CFDRS-BACKWARD-STEP-108 — field-derived reattachment (in progress 2026-08-17)

The provider now owns the backward-facing-step SIMPLE field solve, geometry
mask, and explicit boundary contract. Its normalized parabolic inlet is
applied only to fluid cells, and its lower-wall velocity gradient supplies
signed downstream shear samples. Reattachment is the first downstream
negative-to-non-negative shear crossing, with linear interpolation between
grid columns. Fields without such a crossing return a typed invalid-input
error; the consumer adapter no longer owns a streamfunction field solver or a
correlation result.

The focused source and integration tests were updated to assert the crossing,
boundary values, finite physical position, and exact error for a field without
reattachment. The remaining runtime defect was in the provider execution path:
each momentum solve rebuilt GMRES workspace and copied CSR structure, while the
backward-step pressure correction needed a measured inner sweep cap. The
provider now reuses `KrylovWorkspace`, consumes Leto's mutable CSR parts
directly, exposes a validated consumer-owned SOR override, and the benchmark
declares `pressure_sor_relaxation = 1.7` with `pressure_sweep_cap = 33`. The
default Newtonian split remains unchanged for other geometries; the affected
validation binaries pass 16/16 locally, including the formerly timing-out
integration case in 11.407 s. The cfd-math library passes 202/202, cfd-2d's
focused SIMPLE suite passes 5/5, and warning-denied Clippy passes. Hosted
exact-head Rust and Pages verification remains open.

## CFDrs book and example gate — Apollo public-bound dependency

The book's executable-test failure was caused by untyped diagrams, equations,
commands, and non-self-contained workspace API excerpts. Those blocks now have
explicit fence languages; `mdbook test docs/book` and `mdbook build docs/book`
pass locally. This does not substitute for compiling the linked examples.

The linked Cargo example gate is currently blocked at `cfd-3d`: its spectral
consumer names the public `apollo_fft::PlanScratch` bound, while the integrated
Apollo head `c87a1abe` does not expose it. Apollo remote head
`81583aab8b3eb48c96d138e3980e2c554d9d83fa` contains the provider-owned
re-export and public module change. CFDrs remains unadvanced until that
upstream head is merged; no compatibility shim is permitted.

## Bounded Newton-Krylov recovery integration — CFDrs (2026-08-19)

- `crates/cfd-1d/src/solver/core/newton_fallback.rs` is now declared in the
  solver module graph and is called once when the existing bounded-amplitude
  stagnation detector classifies the Picard trajectory as stalled.
- The recovery budget derives its warm-up and Newton/Krylov limits from
  `SolverConfig.max_iterations`, so the fallback does not create a second
  unbounded solver attempt. `cfd-math::JfnkSolver::solve_checked` propagates
  residual callback failures as typed errors instead of panicking.
- Static source reachability, package formatting, and diff checks are clean.
  Locked local Cargo verification is blocked before compilation by the shared
  Atlas overlay/lock mismatch; provider hosted exact-head Rust and book-figure
  gates remain open until the source branch is checked.

## Lint-floor gate — partial closure (2026-08-06)

- The root workspace now owns the Atlas lint floor; every CFDrs workspace
  package and `xtask` inherit it.
- Passing evidence: `cargo clippy -p xtask --all-targets`,
  `cargo clippy -p cfd-core --lib`, `cargo nextest run -p cfd-core --lib`
  (246/246), `cargo test --doc -p cfd-core` (3/3), and
  `cargo run --manifest-path xtask/Cargo.toml -- legacy-migration-audit`
  (zero legacy dependencies, zero legacy source tokens, clean allowlist).
- Residual: the workspace all-target Clippy gate is not green because existing
  cfd-math unwrap/output sites, cfd-schematics missing docs, cfd-core test/bench
  lint debt, and unrelated format debt still need ratchet increments. This
  increment does not claim full workspace closure.

## cfd-math coarsening ordering slice — committed follow-up (2026-08-06)

- Multigrid coarsening no longer unwraps floating-point partial comparisons;
  finite values sort before unordered values, preserving deterministic
  diagnostics for NaN inputs. The regression test covers finite/NaN ordering.
- `cargo nextest run -p cfd-math --lib` passes 198/198. Focused cfd-math
  library Clippy remains red at 48 existing diagnostics after this slice.

## cfd-math hierarchy and storage invariants — follow-up (2026-08-06)

- Multigrid hierarchy/interpolation state now uses invariant-checked
  expectations instead of bare unwraps. JFNK and spectral operations use named
  C-contiguous storage helpers or explicit connectivity invariants.
- `cargo nextest run -p cfd-math --lib` passes 198/198. Focused cfd-math
  library Clippy remains red at 22 existing diagnostics after this slice.

## cfd-math diagnostics and DG output — closure slice (2026-08-06)

- Performance-monitor mutex accesses now carry invariant diagnostics and its
  calibration messages use structured tracing. DG progress, warnings, and
  completion metrics also use structured tracing instead of stdout writes.
- `cargo clippy -p cfd-math --lib` passes with the Atlas floor; cfd-math
  Nextest remains 198/198. Workspace all-target closure is still open on
  cfd-schematics docs, cfd-core test/bench lint debt, and format debt.

## cfd-schematics topology model documentation — partial closure (2026-08-07)

- The exported topology specification contract now documents its model types,
  fields, enum variants, aliases, and route lookup methods in
  `src/topology/model.rs`.
- `cargo clippy -p cfd-schematics --lib` still fails only on the existing
  documentation floor at this stage; diagnostics decrease from 712 to 611.
- Residual: 611 missing-documentation diagnostics remain across the package's
  other modules. No crate-wide lint suppression was added.

## cfd-schematics constants and config manifests — partial closure (2026-08-07)

- The `ConstantsRegistry` getters and fields, adaptive primitive constants, and
  public config module manifests now carry API documentation.
- `cargo clippy -p cfd-schematics --lib` still fails only on the existing
  documentation floor at this stage; diagnostics decrease from 611 to 534.
- Residual: 534 missing-documentation diagnostics remain across geometry,
  domain, interface, infrastructure, and topology modules outside this scope.

## cfd-schematics geometry builders — partial closure (2026-08-07)

- The public node and channel builder setters now document their domain
  effects, including names, geometry, visual roles, therapy zones, and Venturi
  metadata.
- The package documentation-floor residual decreases from 534 to 524
  diagnostics. No runtime behavior or builder defaults changed.

## cfd-schematics geometry-generator entrypoints — partial closure (2026-08-07)

- The metadata configuration, generation entry points, and fluent builder now
  document their public contracts, including metadata, topology, lineage, and
  rendering inputs.
- The package documentation-floor residual decreases from 524 to 508
  diagnostics. No runtime behavior or generation defaults changed.

## cfd-schematics linear geometry generators — partial closure (2026-08-07)

- The series and parallel geometry builders now document their specification-
  driven blueprint construction contracts.
- The package documentation-floor residual decreases from 508 to 506
  diagnostics. No runtime behavior or generated geometry changed. The current
  working tree reports 492 because peer-owned documentation changes are also
  present in `domain/model/blueprint/analysis_impl.rs` and remain uncommitted.

## cfd-schematics selective-tree generator — partial closure (2026-08-07)

- The public selective-tree path specification, topology variants, request
  fields, and generator entrypoint now document their physical and topology
  contracts.
- The package documentation-floor residual decreases from 506 to 468
  diagnostics. The current working tree reports 454 because peer-owned
  documentation changes remain uncommitted in
  `domain/model/blueprint/analysis_impl.rs`.

## Delivery closure — external RecurseML status (2026-08-06)

CFDrs PR #325 is merged at `fa29c5174c29ac84f5c14e385b4c73866164f712`.
The typed Aequitas metric implementation and repository-owned verification are
green: the book-figure SSOT check passes, and the focused cfd-schematics and
cfd-1d gates are recorded below. The historical `recurseml/analysis` error
for the pre-merge range `77e8a77f..a7159b63` was report-only and exposed no Rust
diagnostic; it is retained as historical evidence, not an open delivery
blocker. No CFDrs Aequitas metric gap remains in this slice.

## Venturi geometry metadata metrics (CFDRS-AEQ-MET-47, 2026-08-05)

A fresh public-contract scan found a residual unit gap below the schematic-mesh
geometry closure: `VenturiGeometryMetadata` exposed channel widths, lengths,
angles, and throat position as scalar metadata, while `ChannelVenturiSpec`
returned scalar cavitation dose and pressure-drop results. The direct 1D
coefficient builder, 2D projection, schematic interchange, mesh converter,
examples, and physics fixtures consumed those fields as untyped values.

The gap is closed. `VenturiGeometryMetadata` now carries Aequitas `Length`,
`Angle`, and `Dimensionless`; `ChannelVenturiSpec` carries typed `Length` and
returns `Dimensionless` cavitation dose and `Pressure` pressure drop from typed
`VolumetricFlowRate`, `MassDensity`, and `Dimensionless` inputs. FDA cavitation
compliance now carries typed Mach-index results and accepts typed pressure and
density. Angle fields use canonical names and store Eunomia base radians.
Scalar extraction is limited to 1D/2D numerical formulas, mesh/interchange
conversion, and other explicit representation boundaries. No compatibility
adapter remains.

The Venturi geometry is real-valued under Eunomia. No complex or imaginary SI
quantity applies; complex values remain reserved for genuine phasor/Fourier
fields. Topology authoring structs remain a separate follow-up domain: their
raw degree/metre inputs are converted at the typed metadata boundary.

Static verification passes with `cargo metadata --offline --locked
--no-deps`, targeted rustfmt, `git diff --check`, and residue scans. The
warning-denied `cfd-schematics` Clippy gate passes, and its focused Nextest
passes 158/158 (`47490ed6-44c0-448c-af2c-dad01cdf5c18`). The `cfd-1d` library
Nextest passes 498/498 (`0ab2db66-61ea-4687-b1cd-0ae0b4d2c7a1`) and its
`blueprint_metadata_physics` integration target passes 5/5
(`fc331715-f569-434e-ac39-878a71527903`). Focused source checks pass for
`cfd-1d`, `cfd-2d`, and `cfd-schematic-mesh` (library/example). The full
`cfd-schematics` doctest gate ran 15/16; the unrelated
`SmoothTransitionConfig` doctest executable was quarantined by Windows
Defender with OS error 225 and `Trojan:Win32/Wacatac.C!ml` (ThreatID
`2147749372`). A targeted retry was blocked by the shared target queue and
timed out before producing a second Rust result. No Defender exclusion, test
weakening, or source workaround was applied. The doctest result is therefore
a host verification residual, not a remaining metric gap. The unused local
`[patch]` and `profile.test.env` warnings are existing stack configuration
warnings outside this slice; provider/linker failures remain separately
tracked.

## Cross-consumer refresh and cfd-3d facade closure (2026-07-31)

The CFDrs, Helios, and Kwavers re-audit found no new missing Aequitas
dimension in the already-closed CFDrs metric families: schematic volumes,
analytical validation, shared fluid properties, cavitation, vascular, and
transient contracts remain typed at their public physical boundaries. The
real-valued CFD/FEM contracts continue to use Eunomia real scalar traits;
Fourier and phasor fields remain the only complex numerical boundaries, so no
imaginary-unit Aequitas quantity is introduced.

The audit did find and close a source defect in `cfd-3d`: its exported
`turbulence` module contained no-op k-epsilon, k-omega-SST, and Smagorinsky
implementations beside the real Eunomia-backed models under
`physics::turbulence`. The placeholder file is deleted, the crate-level path
now exposes the canonical module, and value-semantic turbulence tests exercise
that public path. The historical `cfd-3d` compilation checklist is reconciled
with the current source. The Atlas-overlay locked form remains blocked before
rustc by its mutable provider lock refresh; the standalone lock and compile
gates are recorded above.

The adjacent `cfd-3d::multiphase::exchange` boundary is also now typed with
Aequitas `Dimensionless`, `MassDensity`, and `DynamicViscosity`. The
interpolation formula extracts base scalars once and returns typed phase
properties; f64 and f32 value regressions cover liquid/gas limits. The mixture
remains real-valued under Eunomia and has no complex or imaginary physical
quantity.

The canonical turbulence metric boundary is now closed by
`CFDRS-AEQ-MET-44`: `cfd-core::TurbulenceModel` returns Aequitas
`KinematicViscosity<T>` (`m²/s`) and `SpecificEnergy<T>` (`J/kg`) values, and
all canonical cfd-3d closures wrap their real scalar results at that boundary.
Solver state remains scalar where transport kernels require it; formulas and
dense-field assertions extract base scalars explicitly. The Aequitas provider
uses the coherent `J/kg` semantic alias, while Eunomia complex values remain
reserved for genuine phasor/Fourier fields. No imaginary physical unit is
introduced.

The configured Atlas-overlay locked compile and focused standalone Nextest
remain environment-limited as described above. This is not a remaining source
metric gap; the hosted standalone graph and Pages/figure gates pass.

Verification residuals are separate from the metric audit: the focused offline
package check is green, the cfd-core/cfd-3d doctests pass, and the migrated
metric suite passes 70/70. The cfd-3d library clippy gate is warning-clean;
full all-targets clippy still reports 47 pre-existing test/validation lint
findings in untouched modules. The cfd-math source lint blockers encountered
while checking the dependency graph were fixed at source, and its test target
compiles with `cargo check --tests`; its Nextest link step remains blocked by
the local MinGW/Clang linker invocation without a diagnostic. Neither residual
is a missing Aequitas metric implementation.

## Schematic volume metric audit (2026-07-31)

`CFDRS-AEQ-MET-43` addresses the remaining untyped volume-summary boundary in
`cfd-schematics` and its `cfd-schematic-mesh` consumer. Public summaries and
mesh traces now carry Aequitas `Length`, `Area`, `Volume`, and `Dimensionless`
values instead of unit-suffixed `f64` fields. The previous parallel
millimetre/microlitre fields are removed; microlitres remain an explicit report
conversion from `Volume` through Eunomia-backed Aequitas units.

Mesh signed-volume values enter through `Volume::from_unit::<CubicMillimeter>`;
relative-error percentages extract base values only inside the percentage
formula. This keeps geometry-provider scalars at the mesh boundary and typed
values through the public diagnostic contract. The slice is real-valued under
Eunomia `RealField`; no imaginary-unit metric exists. Eunomia complex values
remain valid for genuinely complex numerical fields elsewhere, but do not
apply to geometric volume.

The focused package check passes for `cfd-schematics` and
`cfd-schematic-mesh`. Nextest run `0383cb4d-2e80-43e5-a0e0-683025defcbd`
passes 207/207. Targeted Rustfmt checks for the changed summary and example,
`git diff --check`, and the typed-field residue scan pass. Existing
package-wide format drift outside this slice remains separately visible in the
repository baseline.

## Legacy analytical benchmark consolidation (2026-07-31)

`CFDRS-AEQ-MET-42` closes the legacy `analytical_benchmarks` duplication
gap. Its raw-scalar Couette, Poiseuille, and Taylor-Green implementations and
their local tests were removed; `tests/physics_validation.rs` now exercises
the canonical Aequitas-backed `analytical` models. The Ghia lid-driven-cavity
reference is owned by `benchmarks::cavity::LidDrivenCavity` as normalized
dimensionless benchmark data rather than physical-unit model state.

The migrated consumer constructs `Velocity`, `Length`, `PressureGradient`,
`DynamicViscosity`, `KinematicViscosity`, and `Time` explicitly and matches
the typed `PoiseuilleFlowRate` and `TaylorGreenKineticEnergy` result enums.
Scalar extraction remains at analytical mesh coordinates and report assertions.
The canonical models remain real-valued under Eunomia `RealField`; no complex
or imaginary-unit physical quantity applies. Targeted Rustfmt and
`git diff --check` pass. The focused cfd-validation gate remains blocked
before the crate by peer-dirty `cfd-math::linear_solver::block_preconditioner`
unresolved `leto-ops::{OwnedNumericLu, SymbolicLu, factor_symbolic}` imports
on lines 26-28 and missing `factor_sparse_with_symbolic` on line 769; this is
a provider/peer integration residual, not a MET42 diagnostic.

## Non-Newtonian analytical metric audit (2026-07-31)

`CFDRS-AEQ-MET-41` closes the non-Newtonian analytical-validation slice.
`PowerLawPoiseuille` now stores channel geometry and pressure gradient as
Aequitas `Length` and `PressureGradient`; its centerline/profile velocity,
per-width flow rate, wall stress, wall shear rate, and generalized Reynolds
number return `Velocity`, `AreaPerTime`, `Pressure`, `ReciprocalTime`, and
`Dimensionless`. `CassonPoiseuille` carries typed geometry and plug radius and
returns typed velocity, wall stress, and per-width flow rate. The rheology
trait exchanges typed shear rate, stress, and dynamic viscosity.

The power-law consistency coefficient remains a formula-bound
`PowerLawConsistency` newtype. Its SI unit is `Pa·sⁿ`, which varies with the
runtime exponent and cannot be represented by one fixed Aequitas dimension;
assigning `DynamicViscosity` to every exponent would be dimensionally false.
Scalar extraction is confined to the constitutive formula, numerical Simpson
integration, and `AnalyticalSolution` mesh-coordinate boundaries. The direct
Newtonian-limit, typed-derived-metric, shear-thinning, Casson plug, wall, and
flow-rate regressions cover the changed contracts.

The models remain real-valued under Eunomia `RealField`; no complex or
imaginary-unit physical quantity applies. Complex values remain reserved for
phasor/Bessel/spectral boundaries. Targeted Rustfmt and `git diff --check`
pass. The pinned cfd-validation check remains blocked before the crate by
peer-dirty `cfd-math::linear_solver::block_preconditioner` unresolved
`leto-ops::{OwnedNumericLu, SymbolicLu, factor_symbolic}` imports on lines
26-28 and missing `factor_sparse_with_symbolic` on line 769; this is a
provider/peer integration residual, not a MET41 diagnostic.

## Blasius analytical metric audit (2026-07-31)

`CFDRS-AEQ-MET-40` closes the Blasius analytical-validation slice.
`BlasiusBoundaryLayer` now stores free-stream velocity, kinematic viscosity,
fluid density, and streamwise position as Aequitas `Velocity`,
`KinematicViscosity`, `MassDensity`, and `Length`. Local Reynolds number,
boundary-layer/displacement/momentum thickness, shape factor, wall shear
stress, and skin friction return `Dimensionless`, `Length`, `Length`,
`Length`, `Dimensionless`, `Pressure`, and `Dimensionless` respectively.

The wall-shear formula now derives dynamic viscosity as `rho * nu` rather than
silently treating kinematic viscosity as dynamic viscosity at unit density.
Similarity variables are dimensionless and velocity-at-coordinate returns
typed `Velocity`; scalar extraction remains at interpolation and
`AnalyticalSolution` mesh-coordinate boundaries. The model remains real-valued
under Eunomia `RealField`; no complex or imaginary-unit physical quantity
applies. The typed-field residue scan, direct thickness/scaling/pressure
regressions, targeted Rustfmt, and `git diff --check` pass. The pinned
cfd-validation check remains blocked before the crate by peer-dirty
`cfd-math::linear_solver::block_preconditioner` unresolved
`leto-ops::{OwnedNumericLu, SymbolicLu, factor_symbolic}` imports on lines
26-28 and missing `factor_sparse_with_symbolic` on line 769; this is a
provider/peer integration residual, not a MET40 diagnostic.

## Taylor-Green analytical metric audit (2026-07-31)

`CFDRS-AEQ-MET-39` closes the Taylor-Green analytical-validation slice.
`TaylorGreenVortex` now stores length, velocity, kinematic viscosity, and
density as Aequitas quantities. Reynolds number, decay rate, kinetic energy,
and enstrophy retain `Dimensionless`, `ReciprocalTime`, a typed dimensional
energy enum, and `ReciprocalTimeSquared` respectively.

`TaylorGreenDimension` makes the 2D/3D branch explicit. The 2D kinetic-energy
result is `TaylorGreenKineticEnergy::PerDepth(Force)`, while the 3D result is
`TaylorGreenKineticEnergy::Volumetric(Energy)`. Aequitas
`ReciprocalTimeSquared` owns the enstrophy dimension. Scalar extraction remains
at formula, mesh-coordinate, and benchmark/report boundaries.

The model remains real-valued under Eunomia `RealField`; no complex or
imaginary-unit physical quantity applies. Aequitas provider commit `f67462a`
adds the reciprocal-time-squared aliases and dimension-law regression. The
typed-field residue scan, direct 2D/3D metric regressions, benchmark call-site
migration, targeted Rustfmt, and `git diff --check` pass. The pinned
cfd-validation check remains blocked before the crate by peer-dirty
`cfd-math::linear_solver::block_preconditioner` unresolved
`leto-ops::{OwnedNumericLu, SymbolicLu, factor_symbolic}` imports on lines
26-28 and missing `factor_sparse_with_symbolic` on line 769; this is a
provider/peer integration residual, not a MET39 diagnostic.

## Stokes analytical metric audit (2026-07-31)

`CFDRS-AEQ-MET-38` closes the Stokes analytical-validation slice.
`StokesFlow` now stores sphere radius, free-stream velocity, dynamic viscosity,
and fluid density as Aequitas quantities. Drag force, drag coefficient,
Reynolds number, and the spherical stream function return `Force`,
`Dimensionless`, `Dimensionless`, and `VolumetricFlowRate` respectively.
Scalar extraction remains only at the analytical formula and mesh-coordinate
boundaries.

The sphere stream function has volumetric-flow dimensions (`m³/s`) because
its definition is velocity times area. The contract remains real-valued under
Eunomia `RealField`; no complex or imaginary-unit physical quantity applies.

The typed-field residue scan, direct Stokes-law regression, targeted Rustfmt,
and `git diff --check` pass. The pinned cfd-validation test-target check is
still blocked before `cfd-validation` by peer-dirty
`cfd-math::linear_solver::block_preconditioner` unresolved
`leto-ops::{OwnedNumericLu, SymbolicLu, factor_symbolic}` imports on lines
26-28 and missing `factor_sparse_with_symbolic` on line 769; this is a
provider/peer integration residual, not a MET38 diagnostic.

## Couette and Poiseuille analytical metric audit (2026-07-31)

`CFDRS-AEQ-MET-37` closes the Couette and Poiseuille analytical-validation
slice. `CouetteFlow` now stores velocity, gap height, pressure gradient, and
dynamic viscosity as Aequitas quantities and returns typed reciprocal time,
pressure, and dimensionless Reynolds metrics. `PoiseuilleFlow` carries typed
velocity, characteristic length, pressure gradient, and viscosity, returns
typed velocity and Reynolds metrics, and exposes a geometry-specific typed
flow-rate enum.

Aequitas now provides `AreaPerTime` for planar flow per unit width (`m²/s`).
`PoiseuilleFlowRate::PerWidth` uses that quantity and
`PoiseuilleFlowRate::Volumetric` uses `VolumetricFlowRate` (`m³/s`); the
runtime geometry enum no longer forces either result through a raw scalar or
an incorrect fixed dimension. Scalar extraction remains at analytical
formula and mesh-coordinate boundaries.

Couette and Poiseuille remain real-valued under Eunomia `RealField`. No
complex or imaginary-unit metric applies; Eunomia complex values remain
reserved for phasor/Bessel formula boundaries.

The typed-field residue scan, direct literature regressions, targeted
Rustfmt, and `git diff --check` pass. The pinned cfd-validation test-target
check and focused Nextest remain blocked before `cfd-validation` by peer-dirty
`cfd-math::linear_solver::block_preconditioner` unresolved
`leto-ops::{OwnedNumericLu, SymbolicLu, factor_symbolic}` imports on lines
26-28 and missing `factor_sparse_with_symbolic` on line 769; this is a
provider/peer integration residual, not a MET37 diagnostic.

## Temperature-model metric refresh (2026-07-30, CFDRS-AEQ-MET-31)

The public-surface scan found that `PolynomialViscosity`,
`ArrheniusViscosity`, `AndradeViscosity`, and `SutherlandViscosity` still
stored density, reference and offset temperatures, viscosity parameters,
specific heat, thermal conductivity, and sound speed as raw generic scalar
fields. Their public storage and direct calculation methods now use Aequitas
`MassDensity`, `ThermodynamicTemperature`, `TemperatureDifference`,
`ReciprocalTemperature`, `DynamicViscosity`, `SpecificHeatCapacity`,
`ThermalConductivity`, and `Velocity`. Temperature is converted to a scalar
only at the existing Eunomia formula boundary; the `Fluid` trait adapter
converts its legacy scalar condition at entry.

Polynomial coefficient vectors remain scalar formula data because coefficient
unit dimensions vary with polynomial order; treating the vector as one typed
quantity would be dimensionally false. The model contract remains real-valued
under Eunomia `RealField`; complex values and imaginary-unit quantities do not
apply to these real constitutive models.

Focused temperature-model Nextest passes 7/7, cfd-core test-target check and
warning-denied Clippy/all-targets pass, cfd-core doctests pass 3/3, targeted
rustfmt and diff checks pass, and the public typed-field residue scan is clean.
The remaining Aequitas model-property gaps are the ideal-gas and the larger-
vessel blood rheology families, which remain separate dependency-ordered
items.

## Microvascular blood metric refresh (2026-07-30, CFDRS-AEQ-MET-33)

The blood audit found that `FahraeuasLindqvist` still stored vessel diameter,
feed hematocrit, and plasma viscosity as raw generic scalars. The calculator
now stores `Length`, `Dimensionless`, and `DynamicViscosity`, and returns typed
apparent viscosity and tube hematocrit. The empirical Pries and Secomb
correlations extract base scalars only at their formula and micrometre
conversion boundaries. The cfd-1d and PyO3 adapters construct the typed core
contract explicitly and extract only for their existing scalar interfaces.

The model remains real-valued under Eunomia `RealField`: ordered diameter and
hematocrit validation plus real-valued powers do not admit Eunomia
`Complex<T>`. Complex values and an imaginary-unit SI material quantity remain
at phasor, Bessel, Womersley, and spectral boundaries rather than being
introduced into a real blood-property contract.

Focused blood-model Nextest passes 3/3, cfd-core/cfd-1d/cfd-python test-target
checks pass, cfd-1d cell-separation Nextest passes 155/155, and the typed
Fåhræus-Lindqvist public-field residue scan is clean. The remaining Aequitas
model-property gaps are the ideal-gas and larger-vessel blood rheology
families; those are not silently claimed by this microvascular slice.

## Ideal-gas metric refresh (2026-07-30, CFDRS-AEQ-MET-34)

The ideal-gas audit found that `cfd-core::physics::fluid::newtonian::IdealGas`
still stored its fixed-dimension gas constant, heat capacity, reference
viscosity, reference temperature, Sutherland offset, and conductivity
coefficient as raw generic scalars. The model now carries those values with
Aequitas `SpecificHeatCapacity`, `DynamicViscosity`,
`ThermodynamicTemperature`, and `TemperatureDifference`; `properties_at`
returns the existing typed `MassDensity`, `DynamicViscosity`,
`SpecificHeatCapacity`, `ThermalConductivity`, and `Velocity` state.
Pressure is typed at the equation-of-state boundary. The gas constant and
conductivity coefficient use `SpecificHeatCapacity` because both are
dimensionally `J/(kg·K)`; no misleading consumer-local semantic alias was
added.

Scalar extraction is confined to the ideal-gas, Sutherland, conductivity, and
sound-speed formula boundaries. The model remains real-valued under Eunomia
`RealField`: ordering and positivity checks plus real-valued powers do not
admit `Complex<T>`. Complex values and an imaginary-unit SI material quantity
remain at phasor, Bessel, Womersley, and spectral boundaries.

Focused cfd-core Nextest passes 16/16 selected tests, including the ideal-gas
value and invalid-input regressions; the cfd-core test-target check, warning-
denied Clippy, cfd-core doctests 3/3, no-default-features rustdoc, targeted
rustfmt, diff checks, and the typed `IdealGas` public-field residue scan all
pass. The default-feature rustdoc closure remains subject to the unrelated
peer `hephaestus-wgpu` `Send + Sync` error at
`application/elementwise_seam.rs:413`. The remaining CFDrs Aequitas model-
property gap is the larger-vessel blood rheology family tracked by MET35.

## Non-Newtonian metric refresh (2026-07-30, CFDRS-AEQ-MET-32)

The public-surface scan found that the Bingham, Casson, Carreau–Yasuda,
Power-law, and Herschel–Bulkley model structs still stored fixed-dimension
physical parameters as raw generic scalars. Their density, yield stress,
viscosity, shear-rate, relaxation-time, exponent, temperature, activation-
energy, thermal, and acoustic fields now use Aequitas quantities. Scalar
extraction is restricted to the existing formula boundary and legacy solver
trait adapters.

Consistency-index fields remain scalar formula data because the `Pa·s^n`
dimension varies with the runtime flow exponent; assigning one fixed unit would
be dimensionally false. The models remain real-valued under Eunomia
`RealField`; complex values and imaginary-unit SI properties do not apply to
these constitutive state contracts.

Focused non-Newtonian Nextest passes 4/4, cfd-core test-target check and
warning-denied all-targets Clippy pass, cfd-core doctests pass 3/3, targeted
rustfmt and diff checks pass, and the public typed-field residue scan is clean.
No-default-features cfd-core rustdoc also passes. Default-feature rustdoc
remains blocked by the unrelated peer `hephaestus-wgpu` bound error at
`application/elementwise_seam.rs:413`; that peer defect is outside this clean
slice and is not masked by the no-default-features result. At that stage in
the ordered migration, the remaining Aequitas model-property gaps were the
ideal-gas and blood constitutive families; later MET34 and MET35 entries close
those families.

## Solid material metric refresh (2026-07-29, CFDRS-AEQ-MET-28)

The source audit found a public cfd-core material family that was not present
in the earlier Aequitas rows: `SolidProperties` and `ElasticSolid` exposed
density, Young's modulus, thermal conductivity, specific heat capacity, and
thermal expansion as raw scalars. The closure carries `MassDensity`,
`Pressure`, `ThermalConductivity`, `SpecificHeatCapacity`, and
`ReciprocalTemperature` through the public contract. Poisson's ratio and the
derived shear-modulus ratio remain `Dimensionless`; the shear-modulus result
is a typed `Pressure`.

This is a real-valued isotropic-material contract (`T: FloatElement + Copy`),
so Eunomia `Complex<T>` and an imaginary-unit Aequitas quantity are not
applicable. Complex constitutive data, if introduced later, requires a
separate provider-backed contract rather than widening this real trait.

The focused cfd-core check and material tests are the acceptance gates. The
broader fluid `FluidState`/`FluidProperties` and remaining constitutive model
fields still exposed raw SI scalars across many implementors. Subsequent MET29
through MET35 entries close those dependency-ordered model-property families;
the solid-material result itself remains unchanged.

## Eunomia complex compatibility refresh (2026-07-28)

CFDrs uses Eunomia `Complex<T>` inside Womersley/Bessel and spectral stability
formulas, but the audited public physical contracts carry real Aequitas
quantities and return real observables. No public complex pressure, impedance,
or other phasor quantity remains untyped, so no Aequitas complex-unit provider
extension is required in CFDrs. The imaginary component is formula data at the
existing numerical boundary, not a missing SI dimension.

## Blueprint cross-fidelity trace refresh (2026-07-28)

The residual public boundary audit found `cfd-3d::blueprint_integration`
exposed reference density, dynamic viscosity, total flow, channel volume,
pressure drop, velocity, nodal pressure, and nodal flow residuals as raw
`f64` values. Those fields are public computational metrics, not private
solver storage. They now use Aequitas `MassDensity`, `DynamicViscosity`,
`VolumetricFlowRate`, `Volume`, `Pressure`, and `Velocity`; scalar extraction
is limited to the cfd-1d/cfd-2d solver adapters and the existing millimetre
mesh trace boundary. Percent errors, normalized pressure-drop coefficients,
grid sizes, and solver tolerances remain dimensionless or structural values.

The public field names no longer encode units, and all in-tree blueprint
integration tests construct and assert the typed contracts. A follow-up scan
of the separate cfd-2d cell-tracking surface found one additional public
physical family; it is recorded and closed by CFDRS-AEQ-MET-27 below.

The focused value gate is `blueprint_integration` Nextest 6/6 and the package
test-check is green. Library-only warning-denied Clippy is green after fixing
the unit-result, return-binding, and `map_or` diagnostics encountered in the
dependency path. The all-targets Clippy command remains blocked by 47
pre-existing diagnostics in peer-edited cfd-3d validation/test modules; this
is verification debt, not a metric residual, and no peer WIP was modified.

## 2D cell-tracking metric refresh (2026-07-29, CFDRS-AEQ-MET-27)

The follow-up public-surface scan found that `cfd-2d::solvers::cell_tracking`
still exposed cell positions, trajectory time, velocity-field coordinates and
velocities, cell diameter/density, fluid viscosity/density, hydraulic geometry,
and bifurcation geometry as raw scalars. These are physical API contracts, not
dense solver storage or dimensionless routing results.

The closure carries `Length`, `Velocity`, `Time`, `MassDensity`, and
`DynamicViscosity` through the public tracker, and represents routing fractions
and lift coefficients with `Dimensionless`. `TrackedPosition` replaces the
heterogeneous `[x, y, t]` array. Scalar extraction is confined to staggered
grid interpolation, particle-force/integration formulas, and the Pries
formula boundary. Eunomia complex values do not enter this real-valued model,
so no imaginary-unit or complex Aequitas extension is required.

The focused `cargo check -p cfd-2d --tests --offline` gate passed. Cell-tracking
Nextest run `bd8e2be6-da1b-49db-85d5-5715fbdb5638` passed 5/5 tests, with 593
tests skipped by the filter. `cargo test --doc -p cfd-2d --offline` passed 1
doctest with 2 ignored. Warning-denied Clippy is not clean because the
peer-owned `crates/cfd-2d/src/physics/momentum/solve.rs:6` retains an unused
`cfd_math::iterative::IterativeLinearSolver` import. The remaining cfd-math
peer edits and the separate cfd-3d runtime residual are not part of this metric
closure.

## Solver runtime refresh (2026-07-28)

The focused reproduction of `venturi::validation::tests::test_venturi_blood_flow`
still exceeds the committed 30-second Nextest budget on the unchanged FEM
solver path: the established restart-200 setting terminates at 30.107 seconds.
A bounded restart-128 experiment terminates at 30.053 seconds and is rejected;
it does not close the budget and changes no production setting. The broad
cfd-3d validation result remains 291/292 with the same Venturi timeout at
30.663 seconds. The runtime row remains `CFDRS-RUNTIME-001`; it is a solver
performance residual, not an untyped Aequitas metric. The concurrent cfd-math
module-tree edits and 47 all-target Clippy diagnostics remain untouched.

## Aequitas public metric gap audit (2026-07-24)

### Verification refresh (2026-07-27, CFDRS-AEQ-MET-23)

The live cell-separation audit found the remaining public physical boundary
after MET-09: equilibrium lateral position, residual force, Dean drag, direct
margination/cell-interaction inputs, Fahraeus/CFL diameters and widths,
viscosities, shear rate, plasma-skimming diameters, and cross-junction
geometry/flow inputs still crossed the public API as raw SI scalars. The
provider now owns `Force`/`Newton`; the migrated contracts use `Length`,
`MassDensity`, `DynamicViscosity`, `Velocity`, `ReciprocalTime`, and
`VolumetricFlowRate`. Scalar extraction is confined to validation and the
numerical force, transport, viscosity, and conductance formulas. All
unit-suffixed identifiers were removed from the migrated public family without
a compatibility facade. The focused `cfd-1d` package check and its full
Nextest suite pass (736/736, three skipped); focused cfd-validation cell
separation, multi-layer junction, and physics-model validation passes 57/57;
and cfd-1d doctests pass 8/8 with three ignored. Warning-denied cfd-1d
Clippy is green. Shared Atlas patch and `profile.test.env` warnings are
environmental and do not affect these gates.

### Verification refresh (2026-07-28, CFDRS-AEQ-MET-24)

The live cfd-3d VOF audit found a real residual after MET-08: public
cavitation configuration and bubble-dynamics contracts exposed surface
tension, bubble radius, nuclei number density, relaxation time, vapor and
liquid density, vapor pressure, and sound speed as raw SI scalars. The
contracts now carry Aequitas `SurfaceTension`, `Length`, `NumberDensity`,
`Time`, `MassDensity`, `Pressure`, and `Velocity`. The public cavitation step
accepts `Time`; scalar extraction is confined to Rayleigh-Plesset and damage
formula kernels, mesh spacing, and dense pressure/density field boundaries.
Dimensionless inception, void-fraction, and damage statistics remain scalar
because their contracts are not dimensional quantities. No scalar facade was
added. The focused cfd-3d check passes through a command-line local-provider
overlay, and the affected cfd-3d Nextest targets pass 83/83. The standalone
CFDrs lock remains dirty in peer-owned work and the ordinary child build still
resolves duplicate local/git provider identities; those are integration
residuals, not metric gaps.

### Verification refresh (2026-07-28, CFDRS-AEQ-MET-25)

The follow-up cfd-core scan found a shared cavitation boundary that MET-24 did
not cover. Rayleigh-Plesset, Venturi, cavitation-number, nuclei-transport,
phase-transfer, damage, biological-damage, and regime-analysis contracts now
carry Aequitas `Length`, `MassDensity`, `DynamicViscosity`, `SurfaceTension`,
`Pressure`, `Velocity`, `Frequency`, `Time`, `Angle`, `Volume`,
`ThermalDiffusivity`, `MassDensityRate`, `ThermodynamicTemperature`, `Energy`,
and `Dimensionless` quantities where the metric has a declared physical
dimension. Scalar extraction remains at analytical equations, dense fields,
mesh/GPU storage, and model-specific dimensionless coefficients, fractions,
probabilities, and indices. The public cfd-3d closure seam no longer contains
raw pressure/density scalars or zero-valued collapse-rate placeholders: it now
delegates to the cfd-core Rayleigh-Plesset and Schnerr-Sauer models and returns
typed errors for invalid closure inputs.

The cfd-core test-target check passes through the command-line local-provider
overlay, and cfd-core Nextest passes 202/202 with no skips. The cfd-3d
test-target check also passes through the same overlay. The broad cfd-3d
validation run passes 291/292; the sole timeout is the pre-existing
`venturi::validation::tests::test_venturi_blood_flow` runtime residual at
30.663 seconds against the committed 30-second budget. All cavitation, VOF,
closure, robustness, and validation tests in that run pass; the timeout is
tracked under `CFDRS-RUNTIME-001`, not this metric contract. Touched files pass
direct rustfmt and `git diff --check`. The peer-owned CFDrs lock and the
ordinary standalone duplicate local/git provider identity remain integration
residuals until the shared provider graph is reconciled; they are not metric
gaps.

### Verification refresh (2026-07-26)

The former Coeus and missing-cutile path blockers are stale: CFDrs has no
Coeus dependency in its current manifest, and locked `cargo check -p cfd-2d
-p cfd-3d -p cfd-validation` passes through the current local Atlas graph. The
focused producer suite (`cfd-core`, `cfd-1d`, and `cfd-optim`) passes
1,127/1,127 with three skips. The cfd-2d library gate passes 517/517 with one
skip. The broader solver-validation gate retains the eight named solver-heavy
tests at or above the committed 30-second budget; those runtime defects remain
open and are not classified as missing metric types.

Production-library warning-denied Clippy is green. The full all-target warning
gate still reports 47 pre-existing lint findings in untouched cfd-3d solver and
test surfaces; four isolated test findings were fixed in `2474604c`, and the
remaining baseline is not an Aequitas metric gap.

This increment warm-starts cfd-2d momentum solves, preserves face-pressure CSR
topology across iterations, rejects unrecoverable large pressure-solver
failure instead of accepting an incomplete correction, fixes SIMPLEC/PIMPLE
state and symmetry/slip boundary handling, routes cfd-3d bifurcation iterations
through the FEM Picard warm-start path, caches FEM and shear geometry, reuses
AMG hierarchy/workspace state, and removes temporary sparse/GMRES hot-path
allocations. cfd-math Nextest passes 357/357; GMRES-focused tests pass 21/21;
SpMV-focused tests pass 7/7; and the 1 mm cfd-3d Venturi regression passes 1/1
in 18.350 s. The exact 5 mm cfd-3d Venturi regression still times out at
30.042 s, so the broad runtime item remains open. No additional Aequitas metric
definition is required for CFDrs in this audit pass.

### Verification refresh (2026-07-27)

The live source audit closed the solver-state boundary identified after
`CFDRS-AEQ-MET-16`. `Network::pressures` and `Network::flow_rates`, their
accessors and setters, `NetworkState` pressure/flow/time state, and the public
network analysis results now use Aequitas quantities. In-tree cfd-2d,
cfd-optim, cfd-validation, test, benchmark, and example callers extract base
scalars only at numerical formula, mesh/GPU, or explicit reporting boundaries.

Residual vectors and `last_residual_norm` remain scalar by design: their units
depend on the assembled equation and scaling policy, so assigning one SI
dimension would misrepresent the solver contract. The classification and
conversion rules are recorded in
[`network-state-metrics.md`](docs/atlas-migration/network-state-metrics.md).
No Aequitas provider addition is required; the existing `Pressure`,
`VolumetricFlowRate`, `Time`, and derived hydraulic dimensions cover the
consumer contract.

The focused cfd-1d gate passes 731/731 with three skips in 23.458 s, and cfd-1d
doctests pass 8/8 with three ignored. The package test-target and example checks
also pass. The warning-denied library gate is currently blocked by a concurrent
peer deletion of the cfd-math `linear_solver` and `interpolation` module trees;
that WIP produces unresolved imports before cfd-1d lint can complete. The
broader 825-test validation gate remains open on its unchanged solver-runtime
budget. No CFDrs compatibility shim is introduced.

### Verification refresh (2026-07-27, CFDRS-AEQ-MET-18)

The live metadata audit found a residual public boundary that the previous
solver-state closure did not cover. `NodeProperties` stored node pressure and
temperature as raw `T`, while `NetworkMetadata` stored total volume and both
pressure/temperature ranges as raw `T` values. The implementation now carries
those fields as Aequitas `Pressure`, `ThermodynamicTemperature`, and `Volume`
quantities. Builder setters accept the same typed contracts, and defaults and
range/volume preservation are covered by value-semantic regressions.

`NodeProperties::metadata` remains `HashMap<String, T>` because its keys have
no declared physical dimension; assigning one Aequitas dimension would misstate
the extensibility contract. No compatibility facade or scalar mirror was
introduced. `cargo check -p cfd-1d --offline` passes, and the configured
`cargo nextest run -p cfd-1d --offline` lane passes 735/735 with 3 skipped.
`cargo test --doc -p cfd-1d --offline` passes 8/8 with 3 ignored. The
warning-denied Clippy gate reaches the peer-owned `cfd-math` graph and stops on
the existing missing documentation for
`linear_solver::dense_bridge::solve_leto_csr_with_leto_dense_array`; this is
outside the metadata slice and remains an active peer residual.

### Verification refresh (2026-07-27, CFDRS-AEQ-MET-19)

The transient droplet quantity migration is now behaviorally verified after
the peer cfd-math error-contract repair. `cargo check -p cfd-1d --tests
--offline` passes; `cargo nextest run -p cfd-1d --test
transient_droplet_parity --offline` passes 9/9; and
`cargo nextest run -p cfd-1d --test transient_literature_validation --offline`
passes 5/5. The selected tests retain volume conservation, flow-weighted
branching, pressure-driven advection, occupancy transitions, and the
literature transport relations. The accidental peer-staged cfd-math
additions were preserved in the shared history and their Leto error contract
was repaired forward in `32be436b`; no droplet compatibility facade exists.

The composition temporal/control boundary is now typed. Event timestamps and
simulation intervals use `Time`, hematocrit and CFL use `Dimensionless`, flow
overrides and snapshot flow use `VolumetricFlowRate`, and pressure overrides
use `Pressure`. Solver residuals remain a separate classification because
their units depend on the assembled equation and scaling policy.

### Verification refresh (2026-07-27, CFDRS-AEQ-MET-20)

The transient composition control migration is behaviorally verified. The
locked cfd-1d test-target check passes; `cargo nextest run -p cfd-1d --test
transient_composition_parity --offline` passes 21/21, including the typed
control-value regression; transient droplet parity passes 9/9; and literature
validation passes 5/5. Scalar extraction remains at event sorting, solver
application, and CFL formula boundaries. No compatibility facade was added.
The existing `MixtureComposition::fractions` map was the remaining public
dimensionless-storage residual for a follow-on representation slice; solver
residuals remain equation-dependent and are not assigned an SI dimension.

### Verification refresh (2026-07-27, CFDRS-AEQ-MET-21)

The remaining public transient time boundary is now typed. Requested timepoint
arguments for composition and droplet simulation use `Vec<Time<T>>`, and
`SimulationTimeConfig` returns typed result and calculation timepoint vectors.
Private transport kernels retain scalar vectors only after an explicit
conversion at sorting, tolerance, and solver boundaries. No public `Vec<T>`
timepoint contract remains in the audited composition or droplet simulators.

The cfd-1d test-target check passes; package Nextest passes 736/736 with 3
skips; warning-denied all-target Clippy passes; and the transient composition,
droplet, and literature tests remain included in that package gate. Mixture
Solver residuals remain equation-dependent and are not assigned an SI unit.

### Verification refresh (2026-07-27, CFDRS-AEQ-MET-22)

The final audited transient composition metric boundary is now typed.
`MixtureComposition::fractions`, blood-hematocrit construction and accessors,
weighted public blends, approximate-equality tolerances, and edge/node
concentration queries now carry Aequitas `Dimensionless<T>`. Scalar extraction
is confined to normalization, mixture arithmetic, solver transport, and
value-semantic assertions. No public raw dimensionless fraction map or
hematocrit/concentration accessor remains in the audited cfd-1d transient
composition surface. Solver residual norms remain equation-dependent and are
not assigned an SI dimension.

The cfd-1d test-target check passes; package Nextest passes 736/736 with 3
skips; warning-denied all-target Clippy passes; and doctests pass 8/8 with 3
ignored. No compatibility facade was added.

| ID | Evidence | Closure |
|---|---|---|
| `CFDRS-AEQ-MET-27` | `cfd-2d::solvers::cell_tracking` exposed physical positions, velocities, time, cell/material properties, and bifurcation geometry as raw scalars. | **IMPLEMENTED in this increment.** Public contracts use Aequitas `Length`, `Velocity`, `Time`, `MassDensity`, and `DynamicViscosity`; dimensionless routing uses `Dimensionless`; scalar extraction remains at interpolation, numerical, and Pries formula boundaries. Focused cfd-2d checks and cell-tracking Nextest provide the value-semantic gate. See [`cell-tracking-physical-metrics.md`](docs/atlas-migration/cell-tracking-physical-metrics.md). |
| `CFDRS-AEQ-MET-21` | Public transient composition and droplet simulation APIs accepted requested, calculated, and returned timepoints as `Vec<T>` after event/configuration time had been typed. | **IMPLEMENTED and package-verified.** Public timepoint vectors and timing accessors now use `Time<T>`; scalar conversion is private to sorting, tolerance, and solver boundaries. Test-target check passes; cfd-1d Nextest passes 736/736 with 3 skips; warning-denied all-target Clippy passes; no compatibility facade exists. Mixture fraction storage remains separately tracked as dimensionless representation. See [`transient-composition-metrics.md`](docs/atlas-migration/transient-composition-metrics.md). |
| `CFDRS-AEQ-MET-22` | Public transient mixture fraction storage, blood-hematocrit construction/accessors, weighted blends, tolerances, and node/edge concentration queries still exposed raw dimensionless scalars. | **IMPLEMENTED and package-verified.** Public maps, constructors, accessors, blends, and tolerances now use `Dimensionless<T>`; scalar extraction is limited to normalization, arithmetic, transport, and assertions. Test-target check passes; cfd-1d Nextest passes 736/736 with 3 skips; warning-denied all-target Clippy passes; doctests pass 8/8 with 3 ignored. No compatibility facade exists. Solver residuals remain equation-dependent. See [`transient-composition-metrics.md`](docs/atlas-migration/transient-composition-metrics.md). |
| `CFDRS-AEQ-MET-20` | Transient composition events, timing/control configuration, and snapshots exposed activation time, hematocrit, flow, pressure, CFL, and snapshot flow/time as raw scalars. | **IMPLEMENTED and focused-verified.** Public contracts now carry Aequitas `Time`, `Dimensionless`, `VolumetricFlowRate`, and `Pressure`; simulator conversion is confined to numerical boundaries. The cfd-1d test-target check passes; composition parity passes 21/21; droplet parity passes 9/9; literature validation passes 5/5; and the typed control regression proves value preservation. Solver residuals remain a separate equation-dependent classification; no compatibility facade exists. See [`transient-composition-metrics.md`](docs/atlas-migration/transient-composition-metrics.md). |
| `CFDRS-AEQ-MET-18` | `cfd-1d` `NodeProperties` and `NetworkMetadata` exposed pressure, temperature, total volume, and pressure/temperature ranges as raw scalar values after MET-17. | **IMPLEMENTED in this increment.** Public metadata contracts now carry Aequitas `Pressure`, `ThermodynamicTemperature`, and `Volume`; typed builder setters, default values, and range preservation are covered by value-semantic tests. `HashMap<String, T>` remains dimension-unknown metadata by contract. Locked `cfd-1d` check passes; Nextest passes 735/735 with 3 skips; doctests pass 8/8 with 3 ignored. Warning-denied Clippy reaches the peer-owned cfd-math missing-docs residual and is not a clean gate for this slice. |
| `CFDRS-AEQ-MET-19` | The transient droplet public boundary exposed droplet volume, injection/snapshot time, normalized positions, occupancy spans, and split thresholds as raw scalar values. | **IMPLEMENTED and focused-verified.** Public and internal droplet state now carries Aequitas `Volume`, `Time`, and `Dimensionless` quantities, with scalar extraction confined to transport formulas. The same revision passes `cargo check -p cfd-1d --tests --offline`; focused Nextest passes 9/9 droplet-parity tests and 5/5 literature-validation tests. Composition event/time contracts remain a separately tracked residual. See [`transient-droplet-metrics.md`](docs/atlas-migration/transient-droplet-metrics.md). |
| `CFDRS-AEQ-MET-17` | `cfd-1d` `Network` and `NetworkState` exposed pressure, volumetric-flow, and simulation-time state as raw scalar vectors; network analysis returned additional physical hydraulic metrics as raw scalars. Residual norms have model-dependent units. | **IMPLEMENTED in the current typed network-state increment; provider dimensions in `f19ba15`.** Public state and analysis contracts now use Aequitas quantities, and all in-tree callers are migrated without adapters. Residuals remain scalar under the documented equation-dependent classification. Locked library, test-target, and example checks pass; cfd-1d Nextest passes 731/731 with 3 skips in 23.458 s; doctests pass 8/8 with 3 ignored. Warning-denied library Clippy is pending the concurrent cfd-math module deletion being reconciled. See [`network-state-metrics.md`](docs/atlas-migration/network-state-metrics.md). |
| `CFDRS-AEQ-MET-16` | `cfd-1d` network `Edge` and `EdgeProperties` exposed flow rate, linear hydraulic resistance, quadratic hydraulic-loss coefficient, and parallel conductance as raw `T` values after geometry and vascular-result cutovers. `Network` pressure/flow/residual vectors were a separate solver-state boundary with additional unit semantics. | **VERIFIED in `a50a9e91`; provider dimensions in `f19ba15`.** Aequitas `VolumetricFlowRate`, `HydraulicResistance`, `QuadraticHydraulicResistance`, and `HydraulicConductance` now remain at public edge contracts; base scalars are extracted only inside resistance, validation, junction-loss, analysis, and matrix-assembly kernels. Locked `cargo check -p cfd-1d` passes; Nextest passes 731/731 with 3 skipped; cfd-1d doctests pass 8/8 with 3 ignored. The remaining solver-state boundary was audited and closed by MET-17. |
| `CFDRS-AEQ-MET-15` | The remaining vascular metric gap covered Murray optimal bifurcation geometry and Olufsen structured-tree terminal radius/impedance after the Womersley and network boundary migration. | **IMPLEMENTED in `e3b664e5`.** `OptimalBifurcation` now uses Aequitas `Length`, `Angle`, `VolumetricFlowRate`, `Dimensionless`, `DynamicViscosity`, and `Pressure`; `OlufsenParameters` uses `Length` and returns `HydraulicResistance`. Validation and PyO3 consumers convert at explicit boundaries. Provider aliases are from `446eb9f`. Locked `cfd-1d` check passes; Nextest passes 731/731 with 3 skipped; focused Womersley adversarial tests pass 2/2; locked `cfd-validation` check passes. See [`vascular-metrics.md`](docs/atlas-migration/vascular-metrics.md). |
| `CFDRS-AEQ-MET-14` | `cfd-1d` vascular `WomersleyNumber`, `WomersleyFlow`, `WomersleyProfile`, `VesselSegment`, `Bifurcation`, and `BifurcationNetwork` expose length, radius, pressure, density, viscosity, frequency, flow, and derived vascular results as raw scalar values. | **IMPLEMENTED.** Aequitas owns the named pressure-gradient, hydraulic-resistance, hydraulic-inertance, and compliance dimensions in `446eb9f`. The Womersley and bifurcation contracts now carry typed physical inputs/results and extract base scalars only at analytical kernels. Murray/Olufsen completion is recorded in `CFDRS-AEQ-MET-15`. |
| `CFDRS-AEQ-MET-13` | `cfd-1d` `ChannelType::Curved` and `Micromixer` exposed curvature radius, hydraulic diameter, and path length as raw SI scalars at public construction and storage boundaries. | **IMPLEMENTED.** Curved-channel radius and micromixer hydraulic diameter/path length now use Aequitas `Length`; scalar extraction remains at resistance and dynamic-parameter boundaries. In-tree constructors and value-semantic tests are migrated. `cargo check -p cfd-1d --tests` passes; Nextest passes 731/731 with 3 skipped; doctests pass 8/8 with 3 ignored; Rustdoc exits 0 with 11 pre-existing link warnings. Warning-denied Clippy remains blocked only by the pre-existing `cfd-math::matrix_zeros` dead-code warning, and no MET-13 test exceeds the committed runtime budget. See [`curved-micromixer-metrics.md`](docs/atlas-migration/curved-micromixer-metrics.md). Vascular and Womersley metrics remain outside this slice. |
| `CFDRS-AEQ-MET-12` | `cfd-1d` `CrossSection`, `ChannelGeometry`, `Edge`, and `EdgeProperties` exposed channel dimensions, area, hydraulic diameter, and length as raw SI scalars across public network construction and analysis. | **IMPLEMENTED.** Cross-section dimensions/custom area, channel length, edge area, and edge-property geometry now use Aequitas `Length` and `Area`. Scalar extraction remains at resistance, junction-loss, transient-transport, and analysis kernels; blueprint conversion, examples, and test fixtures use explicit adapters. The locked check passes and Nextest passes 729/729 with 3 skipped. Doctest, Rustdoc, Clippy, and runtime-budget limits remain recorded below. See [`channel-geometry-metrics.md`](docs/atlas-migration/channel-geometry-metrics.md). |
| CFDRS-AEQ-MET-11 | cfd-1d RectangularChannel, CircularChannel, PorousMembrane, OrganCompartment, and ChannelProperties exposed linear geometry, roughness, area, or volume as raw SI scalars; Component::volume returned Option<T>. | **IMPLEMENTED.** Component geometry now uses Aequitas Length, area methods return Area, Component::volume returns Volume, and ChannelProperties stores Length. Factory and setter scalar inputs convert immediately at their dynamic-parameter boundary; resistance models extract base scalars only for numerical formulas. cfd-1d check passes; Nextest passes 729/729 with 3 skipped; doctests pass 8/8 with 3 ignored; Rustdoc completes with 11 existing link warnings. Clippy remains blocked by the pre-existing cfd-math matrix_zeros dead-code warning. See docs/atlas-migration/component-geometry-metrics.md. |
| CFDRS-AEQ-MET-10 | cfd-1d SurfaceProperties stored roughness, contact angle, and surface energy as raw SI scalars. cfd-core WettingProperties, FluidSolidInterface, and InterfaceProperties likewise exposed surface tension and static/advancing/receding angles without Aequitas dimensions. | **IMPLEMENTED.** Public channel surface contracts now use Aequitas Length, Angle, and EnergyPerArea; material interfaces use SurfaceTension and Angle, and adhesion returns EnergyPerArea. Scalar extraction is confined to the Darcy resistance and cosine-law boundaries. cfd-core Nextest passes 259/259; cfd-1d Nextest passes 729/729 with 3 skipped; doctests pass 11 with 3 ignored; focused shape-factor validation passes 1/1. The non-Newtonian validation binary exceeds the committed 30-second budget, and Clippy remains blocked by the pre-existing cfd-math `matrix_zeros` warning. See docs/atlas-migration/surface-wetting-metrics.md. |
| `CFDRS-AEQ-MET-09` | `cfd-1d` cell-separation and kappa-aware cascade APIs exposed cell diameter, density, treatment/recovery hydraulic diameter, parent inflow velocity, Zweifach–Fung channel diameter, and optimization stage widths as raw SI scalars. The values cross public constructors, routing validation, and `cfd-optim` blueprint summaries. | **IMPLEMENTED in this increment.** `CellProperties`, `PeripheralRecovery`, `CascadeStage`, the public Zweifach–Fung functions, and `StageBlueprintSeparationSummary` now carry Aequitas `Length`, `MassDensity`, and `Velocity` values. Scalar extraction is confined to validation and numerical formula boundaries; all in-tree tests and blueprint construction are migrated. See [`cell-separation-physical-metrics.md`](docs/atlas-migration/cell-separation-physical-metrics.md). Touched-file rustfmt, scoped diff checks, residue scans, `cfd-1d` Nextest (728/728, 3 skipped), `cfd-optim` Nextest (137/137), and focused cell-separation validation Nextest (16/16) pass. The full `cfd-validation` package gate remains blocked by the committed 30-second budget: `test_venturi_flow_3d` and `microventuri_35um_case_produces_converged_informative_2d_result` timed out at approximately 31 seconds; those tests are outside this slice. |
| `CFDRS-AEQ-MET-23` | The remaining public cell-separation family exposed equilibrium position, residual force, Dean drag, direct fluid/geometry inputs, Fahraeus/CFL diameters and widths, viscosity/shear metrics, plasma-skimming diameters, and cross-junction geometry/flow inputs as raw SI scalars after MET-09. | **IMPLEMENTED; focused value-semantic gates pass.** Aequitas adds `Force`/`Newton`; the public cell-separation family now uses typed `Length`, `MassDensity`, `DynamicViscosity`, `Velocity`, `ReciprocalTime`, and `VolumetricFlowRate` contracts. Scalar extraction remains at validation and numerical formula boundaries; callers, tests, and cross-fidelity validation are migrated with no compatibility facade. Rustfmt, metadata, diff checks, residue scans, `cfd-1d` check, full cfd-1d Nextest (736/736, three skipped), focused cfd-validation Nextest (57/57), cfd-1d doctests (8/8, three ignored), and warning-denied cfd-1d Clippy pass. See [`cell-separation-force-metrics.md`](docs/atlas-migration/cell-separation-force-metrics.md). |
| `CFDRS-AEQ-MET-08` | `cfd-core` selective cavitation, `cfd-1d` Venturi screening, and `cfd-optim` Venturi placement/blueprint metrics discarded Aequitas types at public boundaries for pressure, density, velocity, length, viscosity, radius, and surface tension. | **FOCUSED VERIFIED.** The public physical contracts carry Aequitas quantities and all in-tree constructors are migrated. Scalar conversion remains only in numerical formula kernels and the documented serialized report DTO boundary. Provider support is in Aequitas commits `07e2252` and `6dc68c4`. The producer suite passes 1,127/1,127 with 3 skips; warning-denied all-targets Clippy and doctests pass for cfd-core/cfd-1d/cfd-optim. The broad cfd-3d/cfd-validation gate remains open on eight solver-heavy tests exceeding 30 seconds. See [`venturi-physical-metrics.md`](docs/atlas-migration/venturi-physical-metrics.md). |

The remaining raw SI fields in `cfd-optim::metrics::SdtMetrics` are an explicit
serialized display-unit DTO: its module contract states that typed values are
assembled upstream and converted once for reporting. They are not an
unclassified provider-boundary gap. A fresh source scan after MET-10 identified
public component geometry and network channel-property contracts that still
stored lengths, areas, volumes, and roughness as raw scalars; those are now
closed by MET-11 and MET-12. A follow-up scan now finds remaining public
geometry outside the channel/network boundary: `ChannelType::Curved` radius,
`Micromixer` hydraulic diameter/length, and vascular vessel/Womersley metrics.
MET-13 closes the first two. MET-14 and MET-15 close the vascular boundary,
including Murray and Olufsen metrics. The post-migration scan finds no raw `T`
fields for the audited vascular physical inputs or derived physical results;
scalar extraction remains at numerical formula boundaries only.

| ID | Evidence | Closure |
|---|---|---|
| `CFDRS-AEQ-MET-07` | `cfd-1d/src/physics/hemolysis/mod.rs` exposed wall shear stress and exposure duration as raw `f64` arguments and fields, while the returned Giersiepen/Taskin indices were dimensionless. `cfd-1d` flow analysis, `cfd-optim` reporting, and `cfd-validation` passed those scalars directly. | **IMPLEMENTED and focused-verified.** Giersiepen and Taskin accept Aequitas `Pressure` and `Time`; `HemolysisExposure` stores the same typed inputs, all in-tree callers are migrated, and the formula owner remains cfd-core/local model code. The producer suite passes 1,127/1,127 with 3 skips. The broader cfd-validation gate retains the eight documented 30-second runtime timeouts, outside the typed hemolysis contract. See [`hemolysis-exposure-metrics.md`](docs/atlas-migration/hemolysis-exposure-metrics.md). |
| `CFDRS-AEQ-MET-06` | `cfd-3d::cascade` exposed channel geometry, flow rate, outlet pressure, wall shear, pressure drop, and maximum velocity as raw SI scalars. The inlet calculation already constructed Aequitas area, flow, and velocity internally, so the public boundary discarded the provider types. | **IMPLEMENTED and source-verified.** `CascadeChannelSpec`, `CascadeConfig3D`, `ChannelResult3D`, and `CascadeResult3D` carry Aequitas `Length`, `VolumetricFlowRate`, `Pressure`, and `Velocity`. Serde keeps the established SI scalar wire keys through explicit representation adapters; FEM and mesh code convert only at the scalar numerical boundary. Locked cfd-3d/cfd-validation check passes. The broader package suite remains runtime-blocked by eight solver-heavy tests at 30 seconds. See [`cascade-physical-metrics.md`](docs/atlas-migration/cascade-physical-metrics.md). |

### Verification refresh for CFDRS-AEQ-MET-07

The typed hemolysis slice is verified by locked checks for `cfd-1d`, `cfd-optim`,
and `cfd-validation`, Nextest (728/728 passed, 3 skipped), doctests (8/8 passed,
3 ignored), and production-library warning-denied Clippy. `cargo doc --no-deps`
completes with pre-existing unrelated rustdoc link warnings. All-targets Clippy
reports pre-existing test/bench lint debt outside this slice; it is not claimed
as resolved here.

The result field names retain their established serialized keys for wire
compatibility; their Rust types now carry the unit contract. The remaining
`cfd-3d` venturi/bifurcation generic geometry kernels are separate scalar
algorithm boundaries and are not silently classified as closed by this slice.

- 2026-07-22 (resolved in CFD-BOOK-CLOSEOUT-1): the stale book commit expanded
  source-backed chapter indexes with public types and behavioral contracts that
  do not exist in CFDrs, and its SUMMARY linked directly to
  `../../../parity_artefacts/INDEX.md`. mdBook treated that path as a source
  page and overwrote tracked archive HTML outside the book on every build.
  The fix-forward retains the twelve non-duplicated pages backed by real examples, restores
  the expanded chapters to their prior source-grounded content, consolidates
  linear-algebra parity on Leto Ops' analytical oracle, and routes archive
  navigation through a local book page. Exact scans find no Rust definitions
  for the rejected contracts (`NonDimParams`, `BoundaryKind`, `GhiaOracle`,
  `ShearReport`, `ScreeningConfig`, `GiersiepenWurzinger`,
  `HemolysisPath`, or `ParetoFront`). The parity HTML blob is byte-identical
  before and after mdBook
  (`85af4889c39f6d03d78b0dfceeb217f5d260efb5`). Evidence tier:
  source-definition audit, cumulative-diff review, successful book rebuild,
  documented example compilation, warning-denied Clippy, configured Nextest
  177/177, and doctests 16/16.

- 2026-07-22 (resolved in CFD-SCHEMATIC-PATH-1): the stale
  `codex/cfd-example-paths` branch contained a valid native-path boundary that
  never reached main. Current main still required ten lossy or fallible UTF-8
  conversions around renderer calls. The recovery ports only that contract
  onto the current tree: renderer traits borrow `Path`, plotting facades accept
  `AsRef<Path>`, sidecar naming stays in `OsStr`, and every live caller passes
  its native path directly. Evidence tier: exact branch/content comparison,
  affected package/example compilation, warning-denied Clippy, and a clean
  renderer conversion scan. Configured Nextest passes all 177
  `cfd-schematics` tests, including native non-UTF-8 path format detection.

- 2026-07-21 (resolved in CFD-IRIS-COLOR-1): `cfd-schematics` duplicated
  Iris's normalized color-law role with a consumer-owned enum and three local
  formulas. Each edge or node lookup also rebuilt a value vector and rescanned
  its full map. The duplicate laws and wrapper enum are deleted; callers use
  Iris `NamedColorMap` directly. `AnalysisOverlay` now lends or owns maps with
  `Cow`, rejects non-finite values at construction, and stores one finite range
  per map. The old render cost was `Theta(E^2 + V^2)` range work with `E + V`
  transient allocations and `Theta(max(E, V))` transient elements; the new
  cost is `Theta(E + V)` construction, zero transient range allocations, and
  expected `O(1)` map/color lookup. This is an asymptotic and allocation proof,
  not a measured speedup claim. Evidence tier: source-level ownership and
  residue audit; focused value-semantic Nextest 176/176; warning-denied
  all-target/all-feature Clippy; affected example compilation; 16 passing
  doctests; warning-denied Rustdoc; and an executed, visually inspected Venturi
  pressure overlay. Major SemVer classification was attempted but its isolated
  temporary graph cannot build cfd-core because existing CFDrs direct pins and
  Proteus/Hephaestus transitive pins select distinct Aequitas and Leto source
  identities. That provider-pin coherence gap is independent of Iris. Kwavers
  volume rendering remains a separately claimed consumer migration.

- 2026-07-20 (resolved in CFD-LAPLACIAN-PROVIDER-1): cfd-math directly
  implemented the two-dimensional CPU Laplacian and cfd-core carried another
  copy as a GPU test oracle, although Hephaestus already owned the WGPU stencil.
  The CPU solver evaluated `-∇²` while the GPU solver evaluated `∇²`.
  Leto now owns the validated spacing, boundary, polarity, and native-precision
  CPU operation; Hephaestus consumes that contract; both CFD solver operators
  select negative polarity. The local formulas are deleted. Evidence tier:
  provider type unification; exact full-grid CPU and real-WGPU regressions;
  configured Nextest 622/622; all-feature and CPU-only checks;
  warning-denied Clippy/Rustdoc; six runnable doctests; and the updated example
  check. Three-dimensional, variable-coefficient, and SIMD diffusion operators
  remain separate contracts outside this slice. `cargo-semver-checks` was
  attempted but blocked in nightly Rustdoc by a long-lived shared-target Leto
  IDE check; the public constructor break remains classified `[major]`.

- 2026-07-17 (resolved stale work; upstream gap open): removing rsparse by
  routing `DirectSparseSolver` through unpreconditioned GMRES is not a valid
  provider migration. The solver chain already attempts GMRES after its exact
  sparse-LU tier, while cfd-2d invokes the direct tier specifically after GMRES
  stagnation or breakdown; the substitution therefore destroys failure-mode
  independence and contradicts the public direct-solver contract. Leto 0.38
  exposes sparse CG and GMRES but no sparse direct factorization. CFDrs retains
  rsparse until upstream item `LETO-SPARSE-DIRECT-1` provides a generic sparse
  direct API and differential conformance. Evidence tier: source-level
  dependency/call-graph inspection, exact tree equivalence to `main`, focused
  value-semantic Nextest (4/4 cfd-math and 1/1 cfd-2d), and warning-denied
  cfd-math Clippy.

- 2026-07-17 (resolved): `GpuContext` acquires and queries through
  Hephaestus's `ComputeDeviceAcquisition` and `ComputeDeviceCapabilities`
  seams after provider release 0.16.1 repaired typed downlevel acquisition.
  The derived seven-storage-binding request preserves the full downlevel
  descriptor; raw adapter and feature methods are deleted. Nextest serializes
  only provider-acquiring tests through `gpu-device`, eliminating process-level
  WGPU device races while CPU tests remain concurrent. Evidence tier:
  compile-time API removal and empty source scan; value-semantic typed-limit
  regression; cfd-core GPU 245/245, cfd-math GPU 362/362, cfd-2d GPU 570/570
  (27 pre-existing skips), root integration 26/26; warning-denied touched
  targets; doctest/rustdoc; and SemVer's expected major-only classification.
  The root all-target example lint baseline is independently tracked by
  CFD-EXAMPLE-CLIPPY-1.

- 2026-07-17 (resolved): root `cfd-suite --all-targets` Clippy no longer
  reports the 29 diagnostics formerly distributed across seven validation
  examples. Four retained examples now execute provider-owned cfd-1d/cfd-2d
  calculations; three unreferenced static reports are deleted rather than
  presenting hardcoded validation outputs. Evidence tier: executable examples
  plus warning-denied all-target Clippy.

- 2026-07-17 (open): `BifurcationSolver3D` builds an unlabeled SDF volume mesh
  but integrates daughter flow only across `outlet_0` and `outlet_1` labels.
  The resulting zero daughter flows are deterministic, so the invalid root FEM
  example is deleted. CFD-3D-BIFURCATION-BOUNDARIES-1 owns the upstream mesh
  terminal-facet contract and cfd-3d flow regressions. Evidence tier: direct
  executable reproduction and source inspection of mesh construction and
  label-based integration.

- 2026-07-17: `cfd-core::compute::gpu::GpuContext::synchronize` now delegates
  completion to Hephaestus `ComputeDevice::synchronize`; `GpuContext` no longer
  exposes raw WGPU device, queue, or limit fields; and cfd-2d creates its
  Poisson solver through `GpuPoissonSolver::from_context`. Evidence tier:
  compile-time provider integration plus GPU-enabled value-semantic regression
  coverage (244/244 cfd-core and 2/2 accelerated cfd-2d nextest),
  warning-denied cfd-core all-target Clippy, and exact source audits with no
  `device.poll`, old Poisson constructor, context device/queue access, or
  public raw-buffer accessor. The remaining adapter/feature introspection risk
  is resolved by the 0.3.0 typed capability boundary above.

- **Verification closure (2026-07-17)**: `cargo doc -p cfd-core --no-deps
  --features gpu --locked` completes warning-clean after the final raw-buffer
  visibility change. The final source state is verified by cfd-core GPU
  nextest 244/244, cfd-2d accelerated nextest 2/2, warning-denied cfd-core/
  cfd-2d Clippy, and the package documentation gate.

- **SemVer classification (2026-07-17)**: Git-baseline semver checks identify
  the intentional removals as three breaking API classes under a minor-change
  assumption. The explicit major-change classification passes. CFDrs remains
  pre-1.0, so the workspace advances from `0.1.0` to `0.2.0` and records the
  migration in the `0.2.0` changelog section without retaining a compatibility
  surface.

- 2026-07-17: Preserved the stale peer's valid Leto source revision
  `6aedde0c7835238867d6f3cd17b030f7e69cb6f2`, which is merged on Leto `main`,
  and advanced its Moirai companion pin to merged `main`
  `5ead788c70c728d971237d7afa0b915ea7cf87e3`. Locked metadata, all-feature
  `cfd-schematics` check, and warning-denied Clippy pass; evidence tier is
  compile-time integration.

- 2026-07-16: Updated the workspace Moirai source pin to merged `main`
  `5ead788c70c728d971237d7afa0b915ea7cf87e3`. Locked metadata resolves Moirai
  0.4 and Themis 0.10; `cfd-schematics` all-feature check, warning-denied
  Clippy, focused nextest, doctests, and docs pass. The source- and
  test-level evidence is compile-time integration plus value-semantic test
  coverage.
- 2026-07-16: Removed the `cfd-schematics` strict-Clippy baseline in the
  touched test/example cone. Direct geometry/phase values use exact
  bit-pattern assertions, and the Venturi example exposes named physical
  fields instead of positional tuple entries. The workspace-wide formatter
  remains blocked by unrelated pre-existing formatting in
  `crates/cfd-schematics/src/error.rs`; touched files pass `rustfmt`.

- 2026-07-10: Removed stale allowlist entries for `cfd-1d/src/scalar.rs` and
  `cfd-3d/src/scalar.rs`. Both seams are provider-native and contain no legacy
  dependency tokens; the active `cfd-core` compute-dispatch diff is untouched.

## Sprint 2026-07-07: cfd-1d/cfd-3d Eunomia identity seam

- **Resolved direct scalar dependency**:
  `Cargo.toml`, `crates/cfd-1d/Cargo.toml`, and
  `crates/cfd-3d/Cargo.toml` no longer declare direct `num-traits`
  dependencies for the 1D/3D solver scalar seams.
- **Resolved identity ownership**:
  `Cfd1dScalar` and `Cfd3dScalar` now expose `zero()` and `one()` through the
  Eunomia `NumericElement` constants already required by the crate-local
  scalar contracts, removing `num_traits::{Zero,One}` as a supertrait
  requirement.
- **Evidence tier**: compile-time integration, empirical nextest coverage,
  touched-file formatting, and static source/manifest audit. In
  `D:/atlas/repos/CFDrs`, touched-file rustfmt passed; `rustup run nightly
  cargo check -p cfd-1d -p cfd-3d` passed; direct residue scan found no
  `num_traits` or direct `num-traits` hits in the touched manifests/scalar
  cones; and `rustup run nightly cargo nextest run -p cfd-1d -p cfd-3d
  --status-level fail` passed 1122/1122 with one existing slow 3D
  mesh-convergence validation.
- **Residual risk**: package-wide fmt and all-targets clippy remain blocked by
  pre-existing unrelated formatting/lint debt outside this slice. Lockfile
  `num-traits` entries, if present, are transitive provider dependencies owned
  by upstream crates rather than direct CFDrs scalar-seam dependencies.

---

## Sprint 2026-07-05: public sparse/linear-solver Leto boundary

- **Resolved public sparse storage boundary**:
  `crates/cfd-math/src/sparse/mod.rs` now exposes
  `leto_ops::CsrMatrix<T>` as the public `SparseMatrix<T>` alias. Sparse
  builders, assembly helpers, sparse operations, and sparse tests construct
  and operate on Leto CSR directly rather than converting through
  nalgebra-sparse.
- **Resolved public solver vector/matrix boundary**:
  `LinearOperator`, `Preconditioner`, `LinearSolver`, direct solver,
  solver-chain, and migrated solver fixtures use Leto CSR and Leto
  `Array1<T>` boundaries in the requested cone. `cfd-validation::numerical`
  now stores computed and analytical validation vectors as `leto::Array1<T>`
  and uses Leto CSR for linear-solver validation test cases.
- **Evidence tier**: compile-time provider integration, package empirical
  regression tests, clippy, rustdoc, and static source audit. In
  `D:/atlas/repos/CFDrs`, cfd-math check passed; cfd-math test-target check
  passed; cfd-math all-target clippy passed; cfd-math doc passed; cfd-math
  nextest passed 361/361; cfd-validation check passed; cfd-validation
  all-target clippy passed; cfd-validation doc passed after fixing a stale
  intra-doc link; the targeted residue scan found no
  `nalgebra_sparse::CsrMatrix`, public nalgebra-sparse re-export, `DVector`,
  `row_offsets()`, `try_from_csr_data`, `CooMatrix`, or `nalgebra::` matches
  under the migrated sparse/linear-solver/validation files.
- **Residual risk**: The requested sparse/linear-solver public boundary has no
  known nalgebra sparse/vector holdouts in the scanned cone. Full
  `cfd-validation` package nextest remains blocked by the existing venturi
  cross-fidelity convergence tests
  `option2_selected_45um_geometry_routes_to_fallback_and_converges` and
  `microventuri_35um_case_produces_converged_informative_2d_result`, which are
  outside this boundary. Broader CFDrs provider migration still has nalgebra
  residue in other crates and contexts outside this slice.

---

## Sprint 2026-07-04: Solver Chain and FEM Consumer Leto Vector Boundary

- **Resolved chain vector API**:
  `crates/cfd-math/src/linear_solver/chain.rs` now exposes
  `LinearSolverChain::solve` and `solve_with_guess` with `leto::Array1<T>`
  RHS/result vectors instead of nalgebra `DVector`.
- **Resolved 2D direct-fallback consumers**:
  `crates/cfd-2d/src/linear_solver_bridge.rs` is the single cfd-2d conversion
  boundary from nalgebra work vectors into the Leto-backed
  `DirectSparseSolver`; momentum and pressure fallback paths, plus the
  momentum regression test, route through it.
- **Resolved 3D FEM assembly consumers**:
  `crates/cfd-3d/src/fem/leto_bridge.rs` is the single FEM conversion boundary
  for nalgebra work vectors crossing `SparseMatrixBuilder::build_with_rhs` and
  `LinearSolverChain`. `FemSolver` and `ProjectionSolver` now use the
  crate-level `Cfd3dScalar` seam, which carries the Leto real-scalar provider
  bound.
- **Evidence tier**: compile-time provider integration, focused empirical
  nextest, clippy on touched library surfaces, static residue scan, and diff
  hygiene. In `D:/atlas/repos/CFDrs`, `rustup run nightly cargo fmt -p
  cfd-math -p cfd-1d -p cfd-2d -p cfd-3d --check`, `cargo check -p cfd-math
  --no-default-features --lib`, `cargo check -p cfd-1d --no-default-features
  --lib`, `cargo check -p cfd-2d --no-default-features --lib`, `cargo check
  -p cfd-3d --no-default-features --lib`, `cargo nextest run -p cfd-math
  --no-default-features chain direct_solver core_solver simple_gmres
  --status-level fail` (4/4), `cargo clippy -p cfd-math
  --no-default-features --all-targets -- -D warnings`, `cargo clippy -p
  cfd-2d --no-default-features --lib -- -D warnings`, and `cargo clippy -p
  cfd-3d --no-default-features --lib -- -D warnings` passed.
- **Residual risk**: `cargo clippy -p cfd-2d --no-default-features
  --all-targets -- -D warnings` is blocked in `cfd-validation`, not cfd-2d:
  validation benchmark modules still pass nalgebra `DVector` to public
  `cfd_math::sparse::spmv`, and generic validation/1D literature paths need
  the Leto scalar bound propagated after the cfd-1d network solver seam moved.
  The broader iterative solver/preconditioner traits still expose nalgebra
  `DVector` until that trait family moves to Leto arrays.

---

## Sprint 2026-07-04: cfd-1d Eunomia/Leto Scalar Boundary

- **Resolved**: `cfd-1d` no longer exposes scattered nalgebra `RealField`
  imports as its domain scalar contract. The crate now routes domain, network,
  component, resistance, vascular, solver, transient, and analysis generic
  bounds through `Cfd1dScalar`, which explicitly combines the remaining
  nalgebra linear-system backend requirement with the Eunomia scalar provider
  contract consumed by migrated `cfd-core` APIs.
- **Geometry contract corrected**: `NetworkDomain::contains_1d` now accepts
  Leto `Point1<T>`, matching the migrated `cfd-core::geometry::Domain`
  contract.
- **Evidence tier**: static source audit, compile-time integration, full
  empirical nextest, and scoped clippy. `rustup run nightly cargo fmt -p
  cfd-1d --check`, `cargo check -p cfd-1d --no-default-features --lib`, `cargo
  nextest run -p cfd-1d --no-default-features --status-level fail` (725/725, 3
  skipped), and `cargo clippy -p cfd-1d --no-default-features --lib -- -D
  warnings` passed. A direct scan found nalgebra `RealField` only in
  `crates/cfd-1d/src/scalar.rs`, where it documents the remaining matrix
  backend boundary.
- **Residual risk**: `cargo check -p cfd-2d --no-default-features --features
  gpu --lib` now reaches `cfd-2d` and fails on `cfd-2d`'s own nalgebra
  `RealField` bounds around migrated `cfd-core` boundary/fluid APIs. `cargo
  clippy -p cfd-1d --no-default-features --all-targets -- -D warnings` remains
  blocked by pre-existing lint debt in tests, examples, and cell-separation
  modules outside this scalar-boundary slice.

---

## Sprint 2026-07-04: cfd-core GPU Poisson Hephaestus Kernels

- **Resolved**: `crates/cfd-core/src/compute/gpu/poisson_solver.rs` no longer
  owns raw WGPU compute pipelines, bind-group layouts, parameter buffers,
  staging buffers, manual `map_async` readback, or `futures`/`mpsc` mapping
  channels. Jacobi, red-black, and residual entry points now dispatch through
  Hephaestus `WgslMultiStorageKernel`.
- **Provider boundary tightened**: Poisson field/source/residual storage now
  uses `WgpuDevice`'s `ComputeDevice` upload, allocation, and download
  contracts. The public constructor still accepts the existing WGPU device and
  queue handles, then immediately wraps them in a Hephaestus provider.
- **Shape contract corrected**: The solver now stores `nx`, `ny`, `dx`, and
  `dy` from construction and validates `phi.len() == source.len() == nx * ny`
  before dispatch. The previous implementation inferred a square grid from
  `phi.len()`, ignoring the constructor geometry.
- **Evidence tier**: static source audit, compile-time integration, full
  empirical nextest, and scoped clippy. `rustup run nightly cargo fmt -p
  cfd-core --check`, `cargo check -p cfd-core --features gpu`, `cargo check -p
  cfd-core --no-default-features`, `cargo clippy -p cfd-core --features gpu
  --all-targets -- -D warnings`, and full `cargo nextest run -p cfd-core
  --features gpu --status-level fail` (231/231) passed. A direct scan of
  `poisson_solver.rs` found no `create_buffer_init`, `create_buffer(`,
  `params_buffer`, `ComputePipeline`, `BindGroupLayout`, `map_async`,
  `futures::channel`, `poll(wgpu::PollType`, or `std::sync::mpsc` residue.
- **Residual risk**: `cfd-2d` accelerated Poisson consumer verification now
  reaches `cfd-2d` and remains blocked by `cfd-2d` Eunomia/nalgebra trait-bound
  errors before the consumer reaches the GPU Poisson path. Broader CFDrs GPU
  cleanup still has raw WGPU orchestration in non-Poisson kernels and tests.

---

## Sprint 2026-07-04: cfd-validation Direct Vector2 Provider Closure
- **Resolved**: `crates/cfd-validation/src/manufactured/navier_stokes.rs`, `crates/cfd-validation/src/conservation/{momentum.rs,angular_momentum.rs,mod.rs}`, `crates/cfd-validation/src/benchmarks/vorticity_stream.rs`, and MMS/conservation validation tests no longer import or use nalgebra `Vector2`. These surfaces now use `leto::geometry::Vector2`; component access uses Leto `[0]`/`[1]` indexing.
- **Boundary**: This closes direct nalgebra `Vector2` ownership in cfd-validation source and tests. It does not remove nalgebra from cfd-validation because dense `DMatrix`, 3D `Vector3`, and geometry point/vector surfaces remain separate Leto/Gaia provider migrations.
- **Evidence tier**: static source audit, compile-time integration, focused empirical nextest, and scoped clippy. `cargo fmt -p cfd-validation --check` passed. `cargo check -p cfd-validation --no-default-features --tests` passed. `cargo clippy -p cfd-validation --no-default-features --lib -- -D warnings` passed. Focused `cargo nextest run -p cfd-validation --no-default-features manufactured mms conservation taylor momentum angular --status-level fail` passed 104/104 tests. A direct scan over `crates/cfd-validation/{src,tests}` found no `nalgebra::Vector2`, `use nalgebra::Vector2`, or `use nalgebra::{..., Vector2}` residue.
- **Residual risk**: Full `cargo clippy -p cfd-validation --no-default-features --lib --tests -- -D warnings` remains blocked by pre-existing unrelated lints in cross-fidelity, schematics, and other validation tests outside this provider slice. Remaining provider work should move cfd-validation dense matrices to Leto storage and geometry types to Gaia/Leto before nalgebra can be removed from its manifest.

---

## Sprint 2026-07-04: cfd-2d Direct Vector2 Provider Closure
- **Resolved**: `crates/cfd-2d/src/physics/vorticity_stream.rs`, `crates/cfd-2d/src/piso_algorithm/corrector.rs`, `crates/cfd-2d/src/solvers/lbm/solver.rs`, `crates/cfd-2d/src/physics/immersed_boundary.rs`, `crates/cfd-2d/examples/blood_venturi.rs`, and `crates/cfd-2d/benches/solver_benchmarks.rs` no longer import or use nalgebra `Vector2`. These surfaces now use `leto::geometry::Vector2`; component access uses Leto `[0]`/`[1]` indexing. `crates/cfd-validation/src/benchmarks/vorticity_stream.rs` was updated as the downstream public-API consumer of `VorticityStreamSolver::velocity_field`.
- **Boundary**: This closes direct nalgebra `Vector2` ownership for cfd-2d source/tests/examples/benches. It does not remove nalgebra from cfd-2d because scalar `RealField`, boundary `Vector3`, dense `DVector`/`DMatrix`, and nalgebra-sparse matrix surfaces remain separate provider migrations. `physics::immersed_boundary` still uses nalgebra `DMatrix` for force/velocity matrices pending a Leto dense-storage migration.
- **Evidence tier**: static source audit, compile-time integration, focused empirical nextest, and scoped clippy. `cargo fmt -p cfd-2d -p cfd-validation --check` passed. `cargo check -p cfd-2d --no-default-features --examples --benches` passed. `cargo check -p cfd-validation --no-default-features` passed. `cargo clippy -p cfd-2d --no-default-features --example blood_venturi --bench solver_benchmarks -- -D warnings` passed. `cargo clippy -p cfd-validation --no-default-features --lib -- -D warnings` passed. Focused `cargo nextest run -p cfd-2d --no-default-features vorticity corrector lbm immersed --status-level fail` passed 44/44 tests. A direct scan over cfd-2d source/tests/examples/benches plus the touched cfd-validation benchmark found no `nalgebra::Vector2`, `use nalgebra::Vector2`, or `use nalgebra::{..., Vector2}` residue.
- **Residual risk**: Full `cargo clippy -p cfd-2d --no-default-features --all-targets -- -D warnings` remains blocked by pre-existing unrelated lints in examples/tests/modules outside this provider slice. The provider migration still needs cfd-math Leto linear-solver storage, cfd-core boundary vector migration, and Eunomia scalar contract replacement before cfd-2d can remove direct nalgebra/nalgebra-sparse manifest ownership.

---

## Sprint 2026-07-04: cfd-2d Problem/Streamtube Atlas Provider Seam
- **Resolved**: `crates/cfd-2d/src/problem.rs` now stores incompressible
  problem and solution velocity fields with `leto::geometry::Vector2` and
  routes initial pressure, velocity-magnitude maxima, and pressure maxima
  through `crates/cfd-2d/src/scalar.rs`/Eunomia instead of local nalgebra
  vector storage or direct `T::zero()` folds.
  `crates/cfd-2d/src/physics/streamtube/partitioning.rs` no longer uses direct
  `num_traits::{Float,FromPrimitive}`, `T::from_f64(...).unwrap()`,
  `T::zero()`, `T::one()`, `Float::abs`, `Float::sqrt`, or scalar `.abs()` in
  the touched APIs and tests; constants, absolute values, and square roots now
  route through `eunomia::{FloatElement,NumericElement}` and the crate-local
  scalar adapter.
- **Boundary**: This slice is limited to the problem setup and streamtube
  partitioning scalar/vector-provider seam. `problem.rs` still carries
  `nalgebra::RealField` because `cfd_core::physics::{boundary,fluid}` types
  are still nalgebra-bound upstream; removing that bound requires an upstream
  cfd-core provider migration. It does not remove the cfd-2d manifest's direct
  `num-traits` dependency because direct residues remain outside this slice.
- **Evidence tier**: static source audit, compile-time integration, focused
  empirical nextest, and lib clippy. `cargo fmt -p cfd-2d --check` passed.
  `cargo check -p cfd-2d --no-default-features` passed. `cargo clippy -p
  cfd-2d --no-default-features --lib -- -D warnings` passed. Focused `cargo
  nextest run -p cfd-2d --no-default-features problem streamtube separating
  --status-level fail` passed 4/4 tests. `git diff --check` passed for the
  touched problem/streamtube and PM artifact files. A direct-provider scan over
  both touched files found no `num_traits`, `FromPrimitive`, `Float::`,
  `T::from_*`, `.to_f64()`, `T::zero()`, `T::one()`, scalar `.abs()`, or
  local nalgebra `Vector2` residue.
- **Residual risk**: Direct `num-traits` remains in cfd-2d immersed-boundary
  tests, momentum setup/interpolation/boundary, turbulence validation, and
  f64-only/test scalar surfaces; full cfd-2d direct `num-traits` removal
  remains a larger crate-level Eunomia migration before the manifest
  dependency can be dropped. Full `problem.rs` nalgebra-bound removal is blocked
  by upstream `cfd-core` boundary/fluid contracts.

---

## Sprint 2026-07-04: cfd-3d Level-Set Eunomia Scalar Seam
- **Resolved**: `crates/cfd-3d/src/level_set/{weno,advection,solver}.rs` no
  longer import or bound direct `num_traits::{FromPrimitive,Float}`. WENO5-Z
  weights, SSPRK3 coefficients, transport input validation, narrow-band limits,
  and reinitialization/Godunov math now route through `cfd-3d::scalar`, backed
  by Eunomia `FloatElement` and `NumericElement`.
- **Boundary**: This is the level-set module scalar-provider cleanup. Direct
  cfd-3d `num-traits` ownership is now closed by the later root
  lib-test/manifest slice. This entry still preserves the current
  `nalgebra::Vector3` boundary pending the larger Leto/Gaia geometry
  migration.
- **Evidence tier**: static source audit, compile-time integration, focused
  empirical nextest, and lib clippy. `cargo fmt -p cfd-3d --check` passed.
  `cargo check -p cfd-3d --no-default-features` passed. Focused `cargo
  nextest run -p cfd-3d --no-default-features level_set --status-level fail`
  passed 13/13 tests. `cargo clippy -p cfd-3d --no-default-features --lib --
  -D warnings` passed. A targeted scan over `crates/cfd-3d/src/level_set` and
  `crates/cfd-3d/src/scalar.rs` found no direct `num_traits`,
  `FromPrimitive`, `ToPrimitive`, `num_traits::Float`, `T::zero()`,
  `T::one()`, or `Float::` residue.
- **Residual risk**: `cargo clippy -p cfd-3d --no-default-features
  --all-targets -- -D warnings` is still blocked by pre-existing lint debt in
  unrelated cfd-3d test/module code (`poiseuille_test`, `fem_tests`,
  `smagorinsky_test`, `blueprint_integration`, `vof_tests`,
  `robustness_tests`, `bifurcation`, `venturi`, `trifurcation`, and VOF
  modules). Full provider completion still requires the remaining cfd-3d
  direct `num-traits` cleanup plus Leto/Gaia replacement of nalgebra geometry
  and storage surfaces.

---

## Sprint 2026-07-04: cfd-1d Direct num-traits Removal and Resistance/Vascular Eunomia Cleanup
- **Resolved**: `cfd-1d` no longer declares or directly references
  `num-traits`. The resistance scalar contract now uses Eunomia
  `FloatElement`/`NumericElement` for scalar construction, diagnostics, and
  transcendental/math operations. The touched hydraulic resistance,
  serpentine, slug-flow, Bessel/Womersley, structured-tree, bifurcation,
  network blueprint/sink, solver-analysis, and package-test seams no longer
  import `num_traits`, use `FromPrimitive`/`ToPrimitive`, call
  `T::from_f64`/`T::from_usize`/`T::from_u32`, or bridge through
  `nalgebra::try_convert`.
- **Boundary**: This closes direct `num-traits` ownership for the `cfd-1d`
  crate. It does not remove transitive `num-traits` through `approx`,
  `nalgebra`, `nalgebra-sparse`, `half`/Eunomia/Leto/Hephaestus,
  `num-complex`, or other provider stacks. It also does not replace the
  remaining nalgebra/nalgebra-sparse storage and solve boundaries; those stay
  as Leto-backed dense/sparse migration work.
- **Evidence tier**: static source audit, compile-time integration, empirical
  nextest, and dependency-tree audit. `cargo fmt -p cfd-1d --check` passed.
  `cargo check -p cfd-1d` passed. `cargo nextest run -p cfd-1d` passed
  725/725 tests with 3 skipped. A full scan over `crates/cfd-1d/Cargo.toml`,
  `crates/cfd-1d/src`, and `crates/cfd-1d/tests` found no direct
  `num_traits`, `num-traits`, `FromPrimitive`, `ToPrimitive`, `T::from_f64`,
  `T::from_usize`, `T::from_u32`, `nalgebra::try_convert`, or
  `.to_f64().unwrap` residue.
- **Residual risk**: `cargo clippy -p cfd-1d --all-targets -- -D warnings`
  remains blocked by existing all-target lint debt outside this provider
  cleanup, including `blueprint_solve_trace.rs`, `adversarial_tests.rs`,
  `resistance_model_validation.rs`, `medical_millifluidic_screening.rs`,
  `geometry_integration_demo.rs`, cell-separation tests/modules,
  droplet-regime tests, entrance-model tests, matrix-assembly tests, and
  venturi coefficient tests. Broader Atlas migration work remains for
  Leto/nalgebra-sparse storage replacement and Hephaestus higher-level GPU
  kernel ownership.

---

## Sprint 2026-07-04: cfd-1d Domain Components Eunomia Boundary Cleanup
- **Resolved**: `cfd-1d` domain components no longer carry direct
  `num_traits` conversion/math bounds. The component trait pressure-drop
  calculation uses Eunomia `NumericElement::abs`; factory defaults and
  component constants use the existing Atlas provider conversion seam; and
  channel, membrane, mixer, pump, valve, and sensor implementations no longer
  import `FromPrimitive` or `Float`.
- **Boundary**: This slice covers
  `crates/cfd-1d/src/domain/components/{mod,channels,factory,membranes,mixers,pumps,sensors,valves}.rs`.
  It preserves the current nalgebra `RealField` and resistance-model
  boundaries for later Leto/Eunomia work.
- **Evidence tier**: compile-time integration, empirical nextest, static
  source audit, and lint regression audit for the touched file. `cargo fmt -p
  cfd-1d --check` passed. `cargo check -p cfd-1d` passed. `cargo nextest run
  -p cfd-1d` passed 725/725 tests with 3 skipped. A focused component scan
  found no direct `num_traits`, `FromPrimitive`, `ToPrimitive`, `T::from_f64`,
  `T::from_usize`, `Float::`, or generic `.abs()` residue. `cargo clippy -p
  cfd-1d --all-targets -- -D warnings` no longer reports the touched
  `channels.rs` item-order lint.
- **Residual risk**: Full `cfd-1d` all-target clippy remains blocked by
  existing lint debt outside this slice in tests/examples, domain-network,
  cell-separation, vascular, solver-core, and benches. Broader `cfd-1d`
  direct-provider residue remains in solver-core, domain-network,
  vascular Bessel/Womersley, tests, and remaining nalgebra/nalgebra-sparse
  storage boundaries.

---

## Sprint 2026-07-04: cfd-1d Channel/Branching/Analysis Eunomia Boundary Cleanup
- **Resolved**: The next coherent `cfd-1d` provider seam no longer depends on
  direct `num_traits` bounds. Channel flow-regime classification, channel
  flow-resistance constants and powers, channel geometry perimeter math,
  Poiseuille shape factors, branching network solver bounds, and network
  pressure/flow/resistance/performance analysis paths now route scalar
  construction, powers, square roots, absolute values, and scalar-to-f64
  display/oracle conversion through `SafeFromF64`, Eunomia `FloatElement`, and
  Eunomia `NumericElement`.
- **Boundary**: This slice covers `crates/cfd-1d/src/domain/channel`, the
  `domain/junctions/branching` solver/physics/validation cone, and
  `solver/analysis` aggregates/analyzers. It does not remove the direct
  `cfd-1d` manifest `num-traits` dependency because other domains still use
  `FromPrimitive`, `ToPrimitive`, and `Float`.
- **Evidence tier**: compile-time integration, empirical nextest, and static
  source audit. `cargo check -p cfd-1d` passed. `cargo nextest run -p cfd-1d`
  passed 725/725 tests with 3 skipped. Focused source scans found no direct
  `num_traits`, `FromPrimitive`, `ToPrimitive`, `T::from_f64`, generic
  `.to_f64()`, or `Float::` residue in the touched channel/branching/analyzer
  cone or in `solver/analysis`.
- **Residual risk**: Full `cfd-1d` all-target clippy remains blocked by
  unrelated existing lint debt in examples, tests, and cell-separation/
  resistance modules. The separate domain-components direct `num-traits`
  residue is closed by the later domain-components slice; broader direct
  provider residue remains in solver-core, domain-network,
  vascular Bessel/Womersley, and other package areas outside this slice.

---

## Sprint 2026-07-04: cfd-1d Murray's-Law Eunomia Boundary Cleanup
- **Resolved**: `cfd-1d` Murray's-law vascular geometry no longer depends on
  `num_traits::FromPrimitive`. `MurraysLaw` and `OptimalBifurcation` route
  scalar constants, power functions, absolute value, and inverse cosine
  through Eunomia `FloatElement`/`NumericElement`. This uses the current
  Eunomia `FloatElement::acos` surface, which already has value-semantic
  tests in the provider checkout.
- **Boundary**: This slice covers only
  `crates/cfd-1d/src/physics/vascular/murrays_law/{law,bifurcation}.rs` plus
  the required bound propagation in `physics/vascular/bifurcation.rs`.
  `cfd-1d` still declares `num-traits` because many other modules still use
  `FromPrimitive`, `ToPrimitive`, and `Float`.
- **Evidence tier**: static source audit plus formatting. `cargo fmt -p
  cfd-1d --check` passed, and a focused scan over the Murray's-law files found
  no `num_traits`, `FromPrimitive`, `ToPrimitive`, or unqualified generic
  `acos`/`powf`/`powi`/`abs` residue. The provider inverse-cosine contract was
  re-verified in the Eunomia checkout with `cargo nextest run -p eunomia acos`
  (2/2).
- **Residual risk**: `cargo check -p cfd-1d` remains blocked by unrelated
  dirty-tree errors outside this slice, including missing `SafeFromF64` bounds,
  generic `abs`/`powf` ambiguity, and stale `to_f64().unwrap_or(...)` call
  sites. Full `cfd-1d` nextest evidence is pending behind that cleanup.

---

## Sprint 2026-07-03: cfd-core Fåhræus-Lindqvist Eunomia Scalars
- **Resolved**: `FahraeuasLindqvist<T>` no longer depends on direct
  `num_traits::FromPrimitive`, direct generic `T::from_f64()`, direct
  `T::zero()`, direct `T::one()`, direct generic `powf()`, direct generic
  `exp()`, or direct generic `abs()` for local scalar construction and
  microvascular viscosity formulas. Pries/Secomb exponent formulas, `mu_45`,
  relative-viscosity clamping, and tube hematocrit now route scalar constants
  and math through Eunomia `FloatElement`/`NumericElement`.
- **Boundary**: This slice is limited to
  `crates/cfd-core/src/physics/fluid/blood/fahraeus_lindqvist.rs`. Together
  with the prior Cross/Casson/Carreau slices, local blood-model scalar
  `num_traits` construction is closed. The broader `Fluid<T>` trait still
  carries an inherited `nalgebra::RealField` boundary for later provider work.
- **Evidence tier**: compile-time provider integration, existing
  value-semantic blood tests, and static source audit. `cargo check -p
  cfd-core` passed. `cargo nextest run -p cfd-core fahraeus_lindqvist` passed
  3/3. `cargo nextest run -p cfd-core blood` passed 24/24. Touched-file
  rustfmt and touched-file `git diff --check` passed. Focused
  `fahraeus_lindqvist.rs` residue scan found no direct `num_traits`,
  `FromPrimitive`, generic `T::from_f64`, `T::zero`, `T::one`, generic `powf`,
  generic `exp`, or generic `abs` residue. Broader blood residue scan now
  matches only concrete `f64` helper/test expressions.
- **Residual risk**: `cargo clippy -p cfd-core --all-targets -- -D warnings`
  remains blocked by pre-existing unrelated lints in
  `crates/cfd-core/src/physics/boundary/applicator.rs` and
  `crates/cfd-core/src/physics/fluid_dynamics/rhie_chow.rs`.

## Sprint 2026-07-03: cfd-core Casson/Carreau Blood Eunomia Scalars
- **Resolved**: `CassonBlood<T>`, `CarreauYasudaBlood<T>`, and
  `BloodModel<T>` no longer require direct `num_traits::FromPrimitive` for
  local scalar construction or model dispatch. Casson constants, hematocrit
  scaling, temperature correction, square-root apparent-viscosity formula, and
  validation constants now route through Eunomia `FloatElement`/
  `NumericElement`. Carreau-Yasuda constants, zero/one identities, and real
  powers remain on Eunomia after the dispatch-bound fix.
- **Boundary**: This slice is limited to
  `crates/cfd-core/src/physics/fluid/blood/{casson,carreau_yasuda,mod}.rs`.
  The concrete `temperature_viscosity_factor(f64)` helper remains intentionally
  concrete. Fåhræus-Lindqvist and the broader `Fluid<T>` trait still carry
  residual blood-fluid provider migration work.
- **Evidence tier**: compile-time provider integration, existing
  value-semantic blood tests, and static source audit. `cargo check -p
  cfd-core` passed. `cargo nextest run -p cfd-core casson` passed 12/12.
  `cargo nextest run -p cfd-core carreau_yasuda` passed 4/4. `cargo nextest
  run -p cfd-core blood` passed 24/24. Touched-file rustfmt and touched-file
  `git diff --check` passed. Focused touched-file residue scan found no direct
  `num_traits`, `FromPrimitive`, generic `T::from_f64`, `T::zero`, `T::one`,
  generic `powf`, generic `sqrt`, or generic `exp` residue; the only remaining
  match is the intentionally concrete `f64` temperature helper.
- **Residual risk**: `cargo clippy -p cfd-core --all-targets -- -D warnings`
  remains blocked by pre-existing unrelated lints in
  `crates/cfd-core/src/physics/boundary/applicator.rs` and
  `crates/cfd-core/src/physics/fluid_dynamics/rhie_chow.rs`.

## Sprint 2026-07-03: cfd-core Cavitation Eunomia Scalar Cone
- **Resolved**: The touched cavitation Rayleigh-Plesset, biological damage,
  regime-analysis, cavitation-number, and material-damage surfaces no longer
  depend on `nalgebra::RealField`, direct `num_traits::FromPrimitive`, direct
  `T::zero()`, direct `T::one()`, or direct `T::from_f64()` construction.
  Scalar constants, powers, square roots, exponentials, finite checks, and
  scalar min/max now route through Eunomia `FloatElement`/`NumericElement`.
- **Coverage added**: Closed-form value tests now cover cavitation-number
  definition, zero-velocity large-index behavior, pressure-recovery scaling,
  Hammitt erosion pressure-ratio power, Rayleigh collapse impact pressure, and
  pit-depth hardness normalization.
- **Boundary**: This slice is limited to
  `crates/cfd-core/src/physics/cavitation/{rayleigh_plesset,bio_damage,number,damage}.rs`
  and `crates/cfd-core/src/physics/cavitation/regimes/`. Remaining cavitation
  scalar holdouts at that point were `models.rs`, `venturi.rs`,
  `nuclei_transport.rs`, and `heterogeneous_nucleation.rs`; the later Venturi
  slice removed `venturi.rs` from the active residual list.
- **Evidence tier**: compile-time provider integration, value-semantic focused
  tests, and static source audit. `cargo check -p cfd-core` passed. `cargo
  nextest run -p cfd-core cavitation` passed 35/35 tests. Focused residue scan
  across the migrated cavitation files found no `nalgebra`, `RealField`,
  `num_traits`, `FromPrimitive`, `ToPrimitive`, direct `T::from_f64`, direct
  `T::from_u64`, direct `T::zero`, direct `T::one`, direct `powf`, or direct
  `powi` residue. Touched-file rustfmt and touched-file `git diff --check`
  passed. Full `cargo clippy -p cfd-core --all-targets -- -D warnings`
  remains blocked by unrelated existing lints in
  `physics/boundary/applicator.rs` and `physics/fluid_dynamics/rhie_chow.rs`.

## Sprint 2026-07-03: cfd-core Material Eunomia Traits
- **Resolved**: `SolidProperties`, `InterfaceProperties`, `ElasticSolid`,
  `WettingProperties`, and `FluidSolidInterface` no longer depend on
  `nalgebra::RealField`; scalar math and constant construction now use Eunomia
  `FloatElement`/`NumericElement`.
- **Boundary**: `MaterialDatabase` still requires `RealField` because it stores
  `Box<dyn Fluid<T>>`; fluid, hemolysis, fluid-dynamics, boundary, geometry,
  mesh, and solver nalgebra surfaces remain open.
- **Evidence tier**: compile-time integration plus empirical focused tests and
  static source audit. Focused material scan shows no `nalgebra::RealField`,
  `RealField`, `num_traits`, `FromPrimitive`, `ToPrimitive`, or `Float` residue
  in the migrated solid/interface files. `cargo check -p cfd-core` passed, and
  `cargo nextest run -p cfd-core material` passed 4/4 tests. No runtime
  performance claim is made.

## Sprint 2026-07-03: cfd-core Velocity Leto Vector
- **Resolved**: `cfd-core::physics::values::Velocity` now stores
  `leto::geometry::Vector3<T>` instead of `nalgebra::Vector3<T>`, and its
  generic contract is Eunomia `FloatElement`/`NumericElement` instead of
  `nalgebra::RealField`. `PhysicalParameters::gravity` now uses the same Leto
  vector type and no longer requires `RealField` for its own methods.
- **Provider extension**: Leto geometry now derives Serde for `Point2`,
  `Point3`, `Vector3`, `UnitVector3`, and `Isometry3`, preserving CFDrs'
  serialized value-object boundary while using provider-owned vector storage.
- **Boundary**: `ProblemAggregate` and `SimulationAggregate` still carry
  `RealField` because their `Domain<T>` and fluid-property contracts have not
  been migrated in this slice. Material, hemolysis, and fluid-dynamics
  nalgebra-bound surfaces remain open.
- **Evidence tier**: compile-time integration plus static source audit.
  Touched-file rustfmt passed. Focused residue scan shows no nalgebra
  `Vector3`/`RealField` in `Velocity` or `PhysicalParameters`. `cargo check -p
  cfd-core` passed after compiling the modified local Leto provider. `cargo
  nextest run -p cfd-core --lib` passed 183/183 tests. Downstream `cargo check
  -p cfd-2d`, `cargo check -p cfd-3d`, and `cargo check -p cfd-validation`
  passed. No runtime performance claim is made.

## Sprint 2026-07-03: cfd-core Physics Value Eunomia Scalars
- **Resolved**: Scalar-only physics value wrappers no longer use
  `nalgebra::RealField` or `nalgebra::ComplexField`. `Temperature`,
  `Pressure`, `ReynoldsNumber`, and `DimensionlessNumber` now depend on
  Eunomia `FloatElement`/`NumericElement` for scalar construction, zero,
  absolute value, and square root. Immediate aggregate owners now declare the
  same bound where they store these value wrappers.
- **Boundary**: `Velocity` remains on `nalgebra::Vector3`, so its `RealField`
  bound is intentionally preserved for the Leto/Gaia vector replacement slice.
  Additional `RealField` use in material, hemolysis, and broader aggregate
  contracts remains open.
- **Evidence tier**: compile-time integration plus static source audit.
  Touched-file rustfmt passed. A focused scan of the four scalar value-wrapper
  files returns no `nalgebra::RealField`, `RealField`, or `ComplexField`
  matches. `cargo check -p cfd-core` passed. `cargo nextest run -p cfd-core
  --lib` passed 183/183 tests. No runtime performance claim is made.
- **Residual risk**: Full nalgebra removal from `cfd-core` still requires
  replacing vector/matrix contracts and material/hemolysis scalar bounds with
  Atlas-owned providers.

## Sprint 2026-07-03: cfd-1d Non-Python ndarray Path Removal
- **Resolved**: `cfd-1d` no longer declares the unused `sprs` dependency that
  pulled `ndarray v0.17.2` into the active 1D/3D dependency graph. The root
  workspace no longer declares unused `ndarray` as a shared dependency.
  cfd-1d tests that construct `ConstantPropertyFluid::water_20c()` now declare
  the required `eunomia::FloatElement` bound explicitly.
- **Dependency audit**: `cargo tree -p cfd-1d -i ndarray` and `cargo tree -p
  cfd-3d -i ndarray` report no matching package. `cargo tree --workspace -i
  ndarray` shows the remaining workspace path as `ndarray v0.16.1 -> numpy
  v0.22.1 -> cfd-python`.
- **Evidence tier**: manifest/lock static audit plus compile-time integration
  and empirical focused tests. Touched-file rustfmt passed. `cargo update -p
  sprs` removed `sprs`, `ndarray v0.17.2`, and stale transitive packages.
  `cargo check -p cfd-1d` passed. `cargo nextest run -p cfd-1d` passed 725/725
  tests with 3 skipped.
- **Residual risk**: cfd-1d still uses `nalgebra`/`nalgebra-sparse`; replacing
  those surfaces with Leto sparse/vector types remains open. Workspace
  `ndarray` remains through the Python `numpy` boundary, which needs a separate
  binding API decision.

## Sprint 2026-07-03: cfd-math Eunomia Multigrid Interpolation
- **Resolved**: `cfd-math::linear_solver::preconditioners::multigrid::interpolation`
  has been changed to remove direct `num_traits::{FromPrimitive, ToPrimitive}`
  imports and direct scalar conversion/fallback paths. Interpolation scalar
  constants, index-distance conversion, quality row-sum extraction, constant
  preservation error, sparsity ratio, and absolute-value dispatch now route
  through Eunomia helpers.
- **Boundary**: This slice preserves the existing nalgebra `DVector` and
  sparse/vector surfaces. GMG transfer, AMG residual bounds, and full Leto
  migration remain separate provider slices.
- **Evidence tier**: compile-time integration plus empirical focused
  interpolation tests and static source audit. Touched-file rustfmt passed.
  `cargo check -p cfd-math` passed. `cargo nextest run -p cfd-math
  interpolation` passed 14/14 tests. Static scan found no direct `num_traits`,
  `FromPrimitive`, `ToPrimitive`, direct `T::from_f64`, direct
  `T::from_usize`, conversion fallback, `from_f64_or`, `SafeFromF64`, stale
  `rayon`, direct `as f64`, or `.to_f64()` fallback hits in
  `multigrid/interpolation.rs`.
- **Residual risk**: `multigrid/gmg` and `amg.rs` still retain direct provider
  residue; nalgebra sparse/vector surfaces remain open for Leto migration.

## Sprint 2026-07-03: cfd-math Eunomia Multigrid Smoothers
- **Resolved**: `cfd-math::linear_solver::preconditioners::multigrid::smoothers`
  no longer uses direct `T::from_f64(...).unwrap_or_else(...)` conversion
  fallbacks or nalgebra absolute-value dispatch for smoother diagonal
  thresholds, Chebyshev eigenvalue defaults, Chebyshev recurrence constants, or
  smoother update thresholds. The touched AMG owner paths now construct
  coarsening thresholds, smoother relaxation values, and complexity filters
  through Eunomia.
- **Boundary**: This slice preserves the existing nalgebra `DVector` and
  `SparseMatrix` surfaces and keeps AMG's `FromPrimitive` bound because deeper
  coarsening/interpolation contracts still require it.
- **Evidence tier**: compile-time integration plus empirical focused
  value-semantic smoother tests and static source audit. Touched-file rustfmt
  passed. `cargo check -p cfd-math` passed. `cargo nextest run -p cfd-math
  test_gauss_seidel_smoother test_jacobi_smoother
  test_symmetric_gauss_seidel test_sor_smoother test_chebyshev_smoother`
  passed 5/5 tests. Static scan found no direct scalar-conversion fallback,
  stale `SafeFromF64`, `from_f64_or`, direct `T::from_usize`, or stale `rayon`
  hits in the touched smoother/AMG files, and no direct `num_traits` provider
  residue in `smoothers.rs`.
- **Residual risk**: `amg.rs` still imports and bounds direct
  `num_traits::FromPrimitive` for deeper coarsening/interpolation routines.
  Other multigrid modules, the raw GPU operator path, and nalgebra sparse/vector
  surfaces remain open provider-migration work.

## Sprint 2026-07-03: cfd-math Eunomia Stability Scalars
- **Resolved**: `cfd-math::time_stepping::stability` no longer imports
  `num_traits::ToPrimitive` or uses direct
  `T::from_f64(...).unwrap_or_else(num_traits::Zero::zero)` fallback
  conversions for analyzer defaults, CFL thresholds, RK order checks, or von
  Neumann amplification outputs.
- **Boundary**: This slice preserves the existing nalgebra `DMatrix`/`DVector`
  Butcher-tableau API. Leto replacement for that matrix/vector surface remains
  a separate migration item.
- **Evidence tier**: compile-time integration plus empirical value-semantic
  stability tests and static source audit. Touched-file `rustfmt --check`
  passed. `cargo check -p cfd-math` passed. `cargo nextest run -p cfd-math
  stability` passed 5/5 tests. `git diff --check` passed. Static scan found no
  `num_traits`, `ToPrimitive`, `FromPrimitive`, `T::from_f64`, or conversion
  fallback hits under `crates/cfd-math/src/time_stepping/stability`.
- **Residual risk**: Package-level `cargo fmt --package cfd-math --check` is
  still blocked by unrelated existing formatting drift outside the touched
  stability files. Broader cfd-math still contains direct `num-traits`,
  nalgebra matrix/vector surfaces, and GPU provider gaps.

## Sprint 2026-07-03: cfd-core Eunomia Physics Value Boundary
- **Resolved**: `cfd-core` physics value objects (`Velocity`, `Temperature`, `Pressure`, `ReynoldsNumber`, `DimensionlessNumber`) and their dependent management aggregates no longer use direct `num_traits::FromPrimitive` bounds or `T::from_f64(...).unwrap_or_else(T::zero/one)` scalar constant conversions.
- **Boundary**: `nalgebra::RealField` remains in this slice because the touched value objects still store `nalgebra::Vector3` values and use nalgebra vector/scalar methods. The cfd-core manifest still retains `num-traits` because `management::conversion` and other modules outside the touched cone still use direct `num-traits`.
- **Evidence tier**: compile-time integration plus empirical value-semantic regression tests and static source audit. `cargo check -p cfd-core --no-default-features` passed. `cargo nextest run -p cfd-core --no-default-features` passed 144/144 tests. Touched-file `rustfmt --check` passed. Static scan over the touched cone found no `num_traits`, `FromPrimitive`, `T::from_f64`, or silent fallback conversion hits.
- **Residual risk**: Full cfd-core Eunomia migration still needs `management::conversion`, remaining numeric helper modules, and eventual Leto replacement for nalgebra-owned vector/matrix state. Package-level `cargo fmt --check --package cfd-core` is still blocked by pre-existing unrelated formatting drift in `error.rs`, `physics/boundary/ghost_cells.rs`, and `physics/fluid_dynamics/operations.rs`.

---

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

---

# Finding 2026-07-10: cfd-1d double-trifurcation Picard non-convergence (test regression)

- **Symptom**: `cfd-suite` integration test `tests/cross_fidelity_blueprint.rs::cross_fidelity_blueprint_complex_branching` fails. `solve_reference_trace` → cfd-1d `Network2DSolver` returns `MaxIterationsExceeded: Maximum iterations (10000) exceeded` on the `double_trifurcation_cif_venturi_rect` network (blood ρ=1060, µ=3.5e-3, Q=1e-7 m³/s). The sibling `cross_fidelity_blueprint_bifurcation` passes.
- **Isolation**: Reproduces with `-p cfd-suite --no-default-features` (gpu OFF), so it is independent of the concurrent cfd-core `gpu/turbulence_compute` directory-module refactor that currently breaks the `gpu` feature build (`E0583 turbulence_compute` — peer WIP, not this finding).
- **Locus**: `crates/cfd-1d/src/solver/core/mod.rs` Picard fixed-point loop (line ~358) with Anderson acceleration + Newton fallback; `has_converged_dual` never satisfied within 10000 iters for this stiff nonlinear network. Solver correctly returns a typed error (no silent bad result).
- **Hypothesis (unverified)**: migration regression from the Leto-CSR assembly push (`d58d1fe3`/`1d768895`) subtly altering the assembled matrix, or a genuine conditioning/stiffness limit the current Picard+Anderson path cannot crack (cf. recent peer work `0d101352` "enhance Anderson QR collapse detection" — actively-evolving convergence path).
- **DoR to resolve**: (1) differential-test the assembled cfd-1d matrix + rhs for this network against the pre-migration commit (parent of `d58d1fe3`) to classify regression vs. genuine stiffness; (2) capture the Picard residual trajectory (diverge / plateau / slow-decay) via a bounded-iteration diagnostic; (3) if genuine stiffness, verify the Newton fallback engages on Picard stagnation. Blocked from immediate fix: convergence path is under active concurrent peer edit (`0d101352`); coordinate before touching `solver/core`.
- **Evidence tier**: empirical (reproduced deterministically, gpu-independent). Not test-gaming: the test asserts real mass-conservation physics; the fix must be in the solver/assembly, never a weakened tolerance or raised iteration cap.
