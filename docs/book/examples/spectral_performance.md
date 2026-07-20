# Example: spectral_performance

**Crate**: `cfd-suite` (workspace root)  
**Run**: `cargo run --example spectral_performance`  
**Source**: [`examples/spectral_performance.rs`](../../examples/spectral_performance.rs)

## What This Example Demonstrates

Apollo-backed FFT round-trip and spectral-derivative accuracy on a 32-point
sin + 0.5·cos signal, verifying that `cfd-3d`'s spectral wrapper correctly
delegates plan reuse to `apollo-fft`.

| API | Purpose |
|---|---|
| `FourierTransform::<f64>::new(n)` | Build reusable FFT plan |
| `transform.forward(&signal)` | Forward DFT |
| `transform.inverse(&spectrum)` | Inverse DFT |
| `SpectralDerivative::<f64>::new(n)` | Spectral derivative operator |
| `derivative.derivative(&signal, 1)` | First derivative in spectral space |

## Key Code Snippet

```rust
use cfd_3d::spectral::{FourierTransform, SpectralDerivative};
use leto::Array1;

let n = 32;
let transform = FourierTransform::<f64>::new(n)?;

let signal = Array1::from_shape_fn([n], |[i]| {
    let x = 2.0 * PI * i as f64 / n as f64;
    x.sin() + 0.5 * (2.0 * x).cos()
});

let spectrum  = transform.forward(&signal)?;
let recovered = transform.inverse(&spectrum)?;
let err = signal.iter().zip(recovered.iter())
    .map(|(a, b)| (a - b).abs()).fold(0.0_f64, f64::max);
assert!(err < 1e-12, "round-trip error = {err:.3e}");
```

## Atlas Integration

`apollo-fft` owns the reusable FFTW-style plans. `cfd-3d`'s `FourierTransform`
is a thin wrapper that routes into `apollo-fft` plan objects so CFDrs never
allocates planning buffers itself. This is the zero-copy, SSOT design for
spectral methods across the Atlas stack.

## Book Chapter

[← SIMD, GPU, and Backend Migration](../performance_and_atlas.md)
