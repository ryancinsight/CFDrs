# CFDrs

A modular computational fluid dynamics framework in Rust, targeted at microfluidic
and millifluidic devices. It spans three fidelities over one shared schematic model:
1D pipe networks (Kirchhoff / Hagen-Poiseuille), 2D incompressible Navier-Stokes on
structured and unstructured grids, and 3D FEM/spectral/multiphase solvers.

Numeric types are generic over a scalar trait rather than fixed to `f64`, and
geometry, units, and linear algebra come from the shared Atlas stack
(`eunomia`, `leto`, `aequitas`, `apollo`, `moirai`, `hephaestus`).

## Documentation

- [Published CFDrs book](https://ryancinsight.github.io/CFDrs/) — hosted mdBook site.
- [Book source](docs/book/) — Markdown chapters and local build configuration.
- API documentation: `cargo doc --workspace --no-deps --open`.

## Workspace crates

| Crate | Purpose |
|---|---|
| [`cfd-core`](crates/cfd-core/) | Core abstractions, plugin system, fluid properties, boundary conditions, canonical error types |
| [`cfd-math`](crates/cfd-math/) | Linear solvers, multigrid, sparse matrices, time stepping, quadrature, SIMD kernels |
| [`cfd-io`](crates/cfd-io/) | CSV, binary, and checkpoint I/O; optional HDF5 and VTK backends |
| [`cfd-schematics`](crates/cfd-schematics/) | Schematic network design, geometry generation, metadata export, heatmap visualization |
| [`cfd-schematic-mesh`](crates/cfd-schematic-mesh/) | Schematic-to-mesh bridge: blueprint to watertight surface mesh, shell meshing, well plates |
| [`cfd-1d`](crates/cfd-1d/) | 1D network solvers, resistance models, transient composition and droplet tracking |
| [`cfd-2d`](crates/cfd-2d/) | 2D FDM/FVM/LBM, SIMPLE/SIMPLEC/PISO/PIMPLE, turbulence models, schematic projection |
| [`cfd-3d`](crates/cfd-3d/) | 3D FEM and spectral solvers, IBM, level-set/VOF multiphase, cavitation |
| [`cfd-optim`](crates/cfd-optim/) | Design-space search and constraint scoring for millifluidic layouts |
| [`cfd-validation`](crates/cfd-validation/) | Analytical benchmarks, MMS, convergence studies, error metrics, reporting |
| [`cfd-python`](crates/cfd-python/) | PyO3 bindings exposing selected solvers to Python |

`xtask/` holds workspace automation (wheel builds, validation runs, figure checks) and
is not published.

## Requirements

- Rust 1.97.0 — pinned in [`rust-toolchain.toml`](rust-toolchain.toml); `rustup` selects it automatically.
- A GPU backend is optional. The default build enables `cfd-core/gpu`, which pulls in
  `hephaestus-wgpu`; `--no-default-features` builds a CPU-only workspace.

## Quick start

```toml
[dependencies]
cfd-core = { git = "https://github.com/ryancinsight/CFDrs" }
cfd-2d   = { git = "https://github.com/ryancinsight/CFDrs" }
```

```rust
use cfd_2d::grid::{Grid2D, StructuredGrid2D};
use cfd_core::physics::fluid::{ConstantFluid, ConstantPropertyFluid};

fn main() -> Result<(), cfd_core::error::Error> {
    // Water at 20 °C. Properties are range-checked during construction.
    let water = ConstantPropertyFluid::<f64>::water_20c()?;
    let nu = water.dynamic_viscosity().into_base() / water.density().into_base();

    // 64x64 structured grid over the unit square [0,1] x [0,1].
    let grid = StructuredGrid2D::<f64>::unit_square(64, 64)?;
    assert_eq!(grid.num_cells(), 64 * 64);

    assert!((nu - 1.004e-6).abs() < 1e-9);
    Ok(())
}
```

This snippet is the crate-level doctest in [`crates/cfd-2d/src/lib.rs`](crates/cfd-2d/src/lib.rs),
so `cargo test --doc -p cfd-2d` fails if it stops compiling. It covers construction and
property access only — driving a solver to convergence additionally needs a boundary
condition set and an iteration loop; see the `examples/` directory below.

### Install the hooks

```bash
git config core.hooksPath .githooks
```

Git never applies tracked hooks on its own, so this is a one-time step per
clone. The `pre-push` hook runs `scripts/lockfile.py --check`, which is the
same check CI runs. It matters most when working inside the Atlas stack: the
stack's `[patch]` overlay makes cargo resolve first-party dependencies to
local paths and write a `Cargo.lock` with every `source = "git+..."` line
stripped. That lock resolves fine under the overlay and fails every
`--locked` job in CI, so without the hook the corruption is invisible until a
runner reports it. Repair with `python3 scripts/lockfile.py --regenerate`.

## Building

```bash
# Default build (GPU feature on)
cargo build --release

# CPU-only
cargo build --release --no-default-features
```

Optional features:

| Feature | Crate | Effect |
|---|---|---|
| `gpu` (default) | `cfd-core`, `cfd-math`, `cfd-2d` | wgpu device acquisition via `hephaestus-wgpu` |
| `simd` | `cfd-core` | Architecture-conditional SIMD dispatch |
| `experimental` | `cfd-core` | Unstabilized APIs |
| `hdf5` | `cfd-io` | HDF5 output via `consus-hdf5` |
| `vtk` | `cfd-io` | VTK read/write via `ritk-vtk` |

There is no MPI or distributed-memory backend. Parallelism is shared-memory
(`moirai`) plus optional GPU offload.

## Running the gates

```bash
cargo fmt --all --check
cargo clippy --workspace --all-targets -- -D warnings
cargo nextest run --workspace          # per-test budget: 15s slow / 30s terminate
cargo test --doc --workspace           # nextest does not run doctests
cargo doc --workspace --no-deps
cargo run -p xtask -- check-figures    # book figure SSOT check, mirrors CI
```

The per-test budget is committed in [`.config/nextest.toml`](.config/nextest.toml).
CI ([`.github/workflows/ci.yml`](.github/workflows/ci.yml)) runs the figure check;
[`book-pages.yml`](.github/workflows/book-pages.yml) builds and deploys the book, and
[`rust-release.yml`](.github/workflows/rust-release.yml) validates and publishes
individual crates on tagged releases.

Ten of the eleven crates declare `#![warn(clippy::pedantic)]` and then pre-suppress it
with crate-level `#![allow(...)]` attributes — 292 of them across those ten `lib.rs`
files. `cfd-schematic-mesh` declares neither, only `#![forbid(unsafe_code)]`. Those
exemptions are debt being burned down, not a clean pedantic baseline: read the
`lib.rs` header before treating a warning-free clippy run as evidence of pedantic
conformance.

## Examples and benchmarks

```bash
cargo build --examples                 # check every example still compiles
cargo run --release --example cavity_validation
cargo bench --workspace
```

`examples/` declares 37 example programs covering cavity and pipe-flow validation,
blood rheology, cavitation, spectral and FEM solvers, CSG meshing, and GPU detection;
`crates/cfd-2d/examples/` adds the book's worked cases. `benches/` holds 11 criterion
suites. Examples are not run by the test gates, so `cargo build --examples` is the
check that keeps them from rotting.

## Repository layout

- `crates/` — workspace members (table above)
- `examples/`, `benches/`, `tests/` — runnable programs, criterion suites, integration tests
- `docs/book/` — mdBook source; `docs/` also holds ADRs, SRS, and PRD
- `xtask/` — workspace automation
- `backlog.md`, `checklist.md`, `gap_audit.md` — active work tracking and known gaps
- `CHANGELOG.md` — externally observable changes

Known limitations and open validation gaps are tracked in
[`gap_audit.md`](gap_audit.md) and [`backlog.md`](backlog.md) rather than restated here,
so they cannot drift out of date.

## Contributing

- Match the surrounding module structure; new shared logic goes in the deepest crate
  that all its consumers already depend on.
- Value-semantic assertions only — `is_ok()` alone is not a test.
- Public items need Rustdoc; runnable examples become doctests.
- Run the gate sequence above before opening a PR.

## References

- Patankar, S. V. (1980). *Numerical Heat Transfer and Fluid Flow*.
- Versteeg, H. K. & Malalasekera, W. (2007). *An Introduction to Computational Fluid Dynamics*.
- Leonard, B. P. (1979). A stable and accurate convective modelling procedure. *Comput. Methods Appl. Mech. Eng.* 19(1), 59–98.

## License

MIT OR Apache-2.0
