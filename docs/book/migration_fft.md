# Chapter 15 — Apollo: FFT and Spectral Methods

CFDrs's spectral Poisson solver, FFT-based filters, and boundary-layer
analyses migrate from `rustfft` (and `realfft` for real-to-real) to
**Apollo**'s forward FFT crate.  Apollo is *forward-only* (no inverse FFT
because every Atlas consumer can invert by conjugation) and is wired
for autodiff through [`coeus::Tensor`] — so a spectral solver that has
gradient information flows directly into optimization.

## The Apollo Surface

```rust
pub struct FftPlan {
    shape: Vec<usize>,
    inner: PlanInner,
}

pub trait FftEngine {
    fn plan(shape: &[usize]) -> FftPlan;
    fn forward_real(&self, x: &[f64], out: &mut [Complex64]);
    fn forward_complex(&self, x: &[Complex64], out: &mut [Complex64]);
    fn coeus_forward<T: CoeusScalar>(&self, x: &Tensor<T>) -> Tensor<T>;
}
```

A typical CFDrs spectral solve becomes:

```rust
let plan = FftPlan::new(&[nx, ny, nz]);
let mut spectrum = vec![Complex64::zero(); nx * ny * nz];
plan.forward_real(&pressure_field, &mut spectrum);

// Spectral solve: divide by eigenvalue of -Laplacian
for (k, s) in spectrum.iter_mut().enumerate() {
    *s /= 1.0 + eigenvalue(k, nx, ny, nz);
}

// IFFT via conjugation (Apollo is forward-only)
let recovered = plan.forward_complex(&spectrum, ...)
    .map(|c| c.conj());  // + normalization pass
```

## Migration Procedure

| Legacy | Atlas |
|---|---|
| `rustfft::FftPlanner::new().plan_fft(N)` | `FftPlan::new(&[N])` |
| `realfft::RealFftPlanner::new()` | `plan.forward_real(x, out)` |
| manual real-to-complex + zero-pad | `forward_real` (one call) |
| `num_complex::Complex64` | `eunomia::ComplexField` impl |
| `rustfft::inverse` | forward + conjugate (Apollo) |

The legacy `Inverse` path disappears — Atlas handles inverse via forward
+ complex-conjugate, which is mathematically identical and unifies the
solver's autodiff path.

## Autodiff Bridge

When a solver needs a gradient (e.g. for adjoint optimization), Apollo
returns a `Tensor<T>` from `coeus`: the FFT preserves the autodiff tape
so that a downstream loss yields a backward through the FFT.  This is the
critical feature that motivates Apollo's design over `rustfft`.

## Validation Examples

- [`spectral_3d_poisson`](examples/spectral_3d_poisson.md) — 3-D spectral
  Poisson solve, ported from `rustfft` to Apollo.
- [`spectral_performance`](examples/spectral_performance.md) — perf
  benchmark showing identical throughput, simpler integration.
- [`2d_heat_diffusion`](examples/2d_heat_diffusion.md) — repeated
  forward FFTs for an explicit diffusion time-stepper.
- [`matrix_free_demo`](examples/matrix_free_demo.md) — boundary-layer
  spectral preconditioner.

## Further Reading

- [`apollo` source](../../../apollo/)
- [Apollo ↔ Coeus integration notes](../../../coeus/)
