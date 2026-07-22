# Example: adaptive_time_stepping_demo

**Crate**: `cfd-suite` (workspace root)  
**Run**: `cargo run --example adaptive_time_stepping_demo`  
**Source**: [`examples/adaptive_time_stepping_demo.rs`](../../../examples/adaptive_time_stepping_demo.rs)

## What This Example Demonstrates

Three adaptive time-step strategies applied to the test ODE dy/dt = −2y (exact:
y(t) = exp(−2t)), showing accuracy and efficiency tradeoffs.

| Strategy | API | Basis |
|---|---|---|
| CFL-based | `AdaptationStrategy::Cfl` | Stability: Δt ≤ CFL·min(Δx/|u|) |
| Error-based | `AdaptationStrategy::Error` | Richardson extrapolation |
| Combined | `AdaptationStrategy::Combined` | min(CFL-limited, error-controlled) |

## Key Code Snippet

```rust
use cfd_2d::schemes::time::{
    AdaptationStrategy, AdaptiveController, AdaptiveTimeIntegrator, StateVector,
};

let y0 = StateVector::from_shape_vec([1], vec![1.0])?;

// CFL strategy: stable for advection-dominated problems
demonstrate_cfl_adaptation(&y0, t_final)?;

// Error strategy: targets local truncation error ε ≤ ε_target
demonstrate_error_adaptation(&y0, t_final)?;

// Combined: safest for production simulations
demonstrate_combined_adaptation(&y0, t_final)?;
```

## Physics Background

For explicit schemes, **CFL condition** bounds the time step:
```
CFL = |u| · Δt / Δx ≤ CFL_max
```
CFL > 1 causes the numerical domain of dependence to exceed the physical one,
leading to instability. **Richardson extrapolation** estimates the local
truncation error by comparing a full step with two half-steps:
```
ε ≈ |y_full − y_half| / (2^p − 1)
```
where p is the scheme order.

## Book Chapter

[← CFDrs Architecture and Problem Setup](../foundations.md)

