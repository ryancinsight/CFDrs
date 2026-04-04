//! Demonstration of the Apollo-backed spectral solver
//!
//! This example exercises the CFDrs compatibility wrapper around Apollo's FFT
//! plans and reports the observed round-trip error.

use cfd_3d::spectral::{FourierTransform, SpectralDerivative};
use nalgebra::DVector;
use std::f64::consts::PI;

fn main() -> Result<(), Box<dyn std::error::Error>> {
    println!("=== Apollo-backed Spectral Transform Demo ===\n");

    let n = 32usize;
    let transform = FourierTransform::<f64>::new(n)?;
    let signal = DVector::from_iterator(
        n,
        (0..n).map(|i| {
            let x = 2.0 * PI * i as f64 / n as f64;
            x.sin() + 0.5 * (2.0 * x).cos()
        }),
    );

    let spectrum = transform.forward(&signal)?;
    let recovered = transform.inverse(&spectrum)?;
    let max_error = signal
        .iter()
        .zip(recovered.iter())
        .map(|(lhs, rhs)| (lhs - rhs).abs())
        .fold(0.0_f64, f64::max);

    println!("Round-trip max error: {max_error:.3e}");

    let derivative = SpectralDerivative::<f64>::new(n)?;
    let first_derivative = derivative.derivative(&signal, 1)?;
    println!("Sample spectral derivative value: {:.3e}", first_derivative[0]);
    println!();
    println!("Apollo owns the reusable FFT plans; CFDrs keeps the spectral wrapper.");

    Ok(())
}
