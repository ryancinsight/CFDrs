# CFDrs Computational Fluid Dynamics Simulation Suite

CFDrs is a high-performance, pure-Rust computational fluid dynamics library
built on the [Atlas physics stack](https://github.com/ryancinsight).  It
covers the full simulation pipeline -- governing equations, discretisation
schemes, solver kernels, validation oracles, and GPU acceleration -- with a
consistent zero-copy, zero-cost-abstraction design.

## Architecture at a Glance

```
                +-----------------------------------------------------+
                |  cfd-validation  |  cfd-optim  |  cfd-python       |
                |  benchmarks, planning, Python bindings            |
                +-----------------------------------------------------+
                                 |
                +-----------------------------------------------------+
                |  cfd-3d  |  cfd-2d  |  cfd-1d  |  cfd-schematics    |
                |  dimensional scheme matrices, biomedical networks  |
                +-----------------------------------------------------+
                                 |
                +-----------------------------------------------------+
                |  cfd-core  |  cfd-math  |  cfd-io                  |
                |  traits, scalar/math utilities, persistence        |
                +-----------------------------------------------------+
                                 |
                +-----------------------------------------------------+
                |  Atlas: leto | hephaestus | coeus | eunomia | hermes |
                |  moirai | mnemosyne | themis | ritk                 |
                +-----------------------------------------------------+
                                 |
                +-----------------------------------------------------+
                |  CPU: leto + hermes + moirai + eunomia              |
                |  GPU: hephaestus-wgpu + coeus + ritk                |
                |  managed alloc: mnemosyne + themis                 |
                +-----------------------------------------------------+
```

## Atlas Dependencies

| CFDrs need | Atlas crate |
|---|---|
| Scalar field traits | `eunomia` (RealField, FloatElement) |
| Linear algebra, SpMV | `leto` / `leto-ops` (Point, Vector, Array1/2, CsrMatrix) |
| GPU tensors | `hephaestus-core` / `hephaestus-wgpu` |
| Autodiff tensors | `coeus` |
| SIMD kernels | `hermes-simd` |
| Parallel iteration | `moirai-parallel` |
| Memory allocation | `mnemosyne` / `themis` |
| Optimisation / GA | `apollo-optim` |

## Getting Started

```
cargo run --example cavity_validation                  # lid-driven cavity, Re ~ 100
cargo run --example pipe_flow_validation               # Hagen-Poiseuille 1-D benchmark
cargo run --example mesh_3d_integration                # structured 16 x 4 mesh + CSG
cargo run -p cfd-validation --example richardson_convergence   # validation oracle
```

Note: `cfd-suite` examples run at the workspace root via
`cargo run --example <name>`; `cfd-validation` is a workspace
member, so the run form is `cargo run -p cfd-validation
--example richardson_convergence`.  See the `**Run**:` line on
each example's `docs/book/examples/<name>.md` page for the
canonical invocation.

## Chapters

This book progresses from the architectural overview
([Chapter 1](foundations.md)) through governing equations,
turbulence, biomedical flows, numerical methods, geometry, 3-D
flows, validation, performance, and optimisation.  Each chapter
links to runnable examples that exercise the described functionality.

## Figures

The figures linked from the chapter tables of contents are
deterministically generated SVG assets committed with the chapter
text.  Each entry names the example that the figure mirrors;
reproducibility is enforced by the `xtask prebook` manifest at
`figures/MANIFEST.json`.

Run commands below mirror the `**Crate**:` / `**Run**:` headers on
the linked example `.md` page.  All entries point to either
`cfd-suite` (workspace root, run via `cargo run --example ...`) or
`cfd-validation` (member crate, run via `cargo run -p cfd-validation
--example ...`).

- [Plane Poiseuille parabolic profile (1-D analytical)](figures/poiseuille_parabolic_profile.svg)
  -- mirrors `cfd-suite` example `pipe_flow_validation` (the
  Hagen-Poiseuille analytical solution verified by the FEM pipe
  flow example).
- [Lid-driven cavity streamfunction contours (Re ~ 100)](figures/cavity_streamfunction_contour.svg)
  -- mirrors `cfd-suite` example `cavity_validation`.
- [Semi-log residual convergence (SIMPLE vs PISO)](figures/residual_convergence_semilog.svg)
  -- mirrors `cfd-suite` example `turbulent_channel_flow` (pressure
  coupling solver residual trajectory).
- [Structured 16 x 4 channel mesh layout](figures/channel_mesh_layout.svg)
  -- mirrors `cfd-suite` example `mesh_3d_integration` (the structured
  rectangular mesh that integrates with the CSG geometry pipeline).
- [Pipe-flow Reynolds regime map](figures/reynolds_regime_map.svg)
  -- mirrors `cfd-suite` example `pipe_flow_validation` (laminar /
  transitional / turbulent band placement).
- [Richardson extrapolation log-log (2nd-order reference)](figures/richardson_loglog.svg)
  -- mirrors `cfd-validation` example `richardson_convergence`.
- [CFDrs layered architecture on the Atlas stack](figures/architecture_stack.svg)
  -- pure schematic, no example source.
