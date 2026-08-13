# cfd-core

Core abstractions for the [CFDrs](https://github.com/ryancinsight/CFDrs) simulation
suite. Every other `cfd-*` crate depends on this one; it depends on no other
workspace member.

## What it provides

- **`error`** — the canonical `Error` and `Result` types shared across the suite.
- **`physics`** — fluid models (`ConstantPropertyFluid`, blood rheology,
  non-Newtonian, temperature-dependent), boundary conditions, and dimensioned
  values (`Pressure`, `Velocity`, `Temperature`, `ReynoldsNumber`).
- **`abstractions`** — the `Problem` and `SimulationState` contracts.
- **`compute`** — `Solver` / `SolverConfig` traits, time integration, GPU dispatch.
- **`geometry`** — domain descriptions and mesh quality services.
- **`management`** — plugin registry, factories, and simulation aggregates.

Scalars are generic over `eunomia`'s element traits, so solvers built on these
abstractions are not fixed to `f64`.

## Example

```rust
use cfd_core::physics::fluid::{ConstantFluid, ConstantPropertyFluid};

fn main() -> Result<(), cfd_core::error::Error> {
    let water = ConstantPropertyFluid::<f64>::water_20c()?;
    assert!((water.density().into_base() - 998.2).abs() < 1e-9);
    Ok(())
}
```

`prelude` re-exports the commonly used subset:

```rust
use cfd_core::prelude::*;
```

## Features

| Feature | Default | Effect |
|---|---|---|
| `gpu` | yes | wgpu device acquisition via `hephaestus-wgpu` |
| `simd` | no | Architecture-conditional SIMD dispatch |
| `experimental` | no | Unstabilized APIs |

Build CPU-only with `--no-default-features`.

## Documentation

- [CFDrs book](https://ryancinsight.github.io/CFDrs/)
- [Workspace README](https://github.com/ryancinsight/CFDrs#readme)

## License

MIT OR Apache-2.0
