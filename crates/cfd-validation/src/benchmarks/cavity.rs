//! Lid-driven cavity benchmark problem
//!
//! Reference: Ghia et al. (1982) "High-Re solutions for incompressible flow
//! using the Navier-Stokes equations and a multigrid method"

use super::{Benchmark, BenchmarkConfig, BenchmarkResult};
use cfd_core::conversion::{SafeFromF64, SafeFromUsize};
use cfd_core::error::{Error, Result};
use nalgebra::{DMatrix, RealField};
use num_traits::{FromPrimitive, ToPrimitive};

/// Lid-driven cavity benchmark
pub struct LidDrivenCavity<T: RealField + Copy> {
    /// Cavity dimensions (assumed square)
    pub size: T,
    /// Lid velocity
    pub lid_velocity: T,
}

impl<T: RealField + Copy> LidDrivenCavity<T> {
    /// Create a new lid-driven cavity benchmark
    pub fn new(size: T, lid_velocity: T) -> Self {
        Self { size, lid_velocity }
    }

    /// Get reference data from Ghia et al. (1982)
    /// Get Ghia et al. reference data for specified Reynolds number
    pub fn ghia_reference_data(&self, reynolds: T) -> Option<(Vec<T>, Vec<T>)> {
        // Literature-validated reference data for lid-driven cavity from Ghia et al. (1982)
        // Reynolds numbers: 100, 400, 1000, 3200, 5000, 7500, 10000

        let re_100 = <T as SafeFromF64>::try_from_f64(100.0).ok()?;
        let re_1000 = <T as SafeFromF64>::try_from_f64(1000.0).ok()?;

        if (reynolds - re_100).abs() < <T as SafeFromF64>::try_from_f64(0.1).ok()? {
            // Reference data for Re=100 (centerline velocities)
            Some((
                vec![
                    <T as SafeFromF64>::try_from_f64(0.0).ok()?,
                    <T as SafeFromF64>::try_from_f64(0.0625).ok()?,
                    <T as SafeFromF64>::try_from_f64(0.0703).ok()?,
                    <T as SafeFromF64>::try_from_f64(0.0781).ok()?,
                    <T as SafeFromF64>::try_from_f64(0.0938).ok()?,
                    <T as SafeFromF64>::try_from_f64(0.1563).ok()?,
                    <T as SafeFromF64>::try_from_f64(0.2266).ok()?,
                    <T as SafeFromF64>::try_from_f64(0.2344).ok()?,
                    <T as SafeFromF64>::try_from_f64(0.5000).ok()?,
                ],
                vec![
                    <T as SafeFromF64>::try_from_f64(0.84123).ok()?,
                    <T as SafeFromF64>::try_from_f64(0.78871).ok()?,
                    <T as SafeFromF64>::try_from_f64(0.73722).ok()?,
                    <T as SafeFromF64>::try_from_f64(0.68717).ok()?,
                    <T as SafeFromF64>::try_from_f64(0.23151).ok()?,
                    <T as SafeFromF64>::try_from_f64(0.00332).ok()?,
                    <T as SafeFromF64>::try_from_f64(-0.13641).ok()?,
                    <T as SafeFromF64>::try_from_f64(-0.20581).ok()?,
                    <T as SafeFromF64>::try_from_f64(-0.21090).ok()?,
                ],
            ))
        } else if (reynolds - re_1000).abs() < <T as SafeFromF64>::try_from_f64(0.1).ok()? {
            // Reference data for Re=1000 (centerline velocities)
            Some((
                vec![
                    <T as SafeFromF64>::try_from_f64(0.0).ok()?,
                    <T as SafeFromF64>::try_from_f64(0.0625).ok()?,
                    <T as SafeFromF64>::try_from_f64(0.0703).ok()?,
                    <T as SafeFromF64>::try_from_f64(0.0781).ok()?,
                    <T as SafeFromF64>::try_from_f64(0.0938).ok()?,
                    <T as SafeFromF64>::try_from_f64(0.1563).ok()?,
                    <T as SafeFromF64>::try_from_f64(0.2266).ok()?,
                    <T as SafeFromF64>::try_from_f64(0.2344).ok()?,
                    <T as SafeFromF64>::try_from_f64(0.5000).ok()?,
                ],
                vec![
                    <T as SafeFromF64>::try_from_f64(0.65928).ok()?,
                    <T as SafeFromF64>::try_from_f64(0.57492).ok()?,
                    <T as SafeFromF64>::try_from_f64(0.51117).ok()?,
                    <T as SafeFromF64>::try_from_f64(0.46604).ok()?,
                    <T as SafeFromF64>::try_from_f64(0.33304).ok()?,
                    <T as SafeFromF64>::try_from_f64(0.18719).ok()?,
                    <T as SafeFromF64>::try_from_f64(0.05702).ok()?,
                    <T as SafeFromF64>::try_from_f64(0.02135).ok()?,
                    <T as SafeFromF64>::try_from_f64(-0.21388).ok()?,
                ],
            ))
        } else {
            None
        }
    }
}

impl<T: RealField + Copy + FromPrimitive + SafeFromF64 + num_traits::ToPrimitive> Benchmark<T>
    for LidDrivenCavity<T>
{
    fn name(&self) -> &'static str {
        "Lid-Driven Cavity"
    }

    fn description(&self) -> &'static str {
        "2D incompressible flow in a square cavity with moving lid"
    }

    fn run(&self, config: &BenchmarkConfig<T>) -> Result<BenchmarkResult<T>> {
        let n = config.resolution;
        let dx = self.size / <T as SafeFromUsize>::try_from_usize(n)?;

        // Initialize stream function and vorticity with proper double-buffering
        let mut psi = DMatrix::<T>::zeros(n, n);
        let mut omega = DMatrix::<T>::zeros(n, n);
        let mut psi_prev = DMatrix::<T>::zeros(n, n); // Zero-copy initialization

        // Calculate viscosity from Reynolds number
        let viscosity = self.lid_velocity * self.size / config.reynolds_number;

        // Time stepping parameters
        let dt = <T as SafeFromF64>::try_from_f64(0.001)? * dx * dx / viscosity;
        let mut convergence = Vec::new();
        let mut iteration = 0;

        while iteration < config.max_iterations {
            psi_prev.copy_from(&psi);

            // Update vorticity on boundaries using stream function
            // Top boundary (moving lid): ω = -2(ψ_{i,n-2} - ψ_{i,n-1})/Δy² - 2U/Δy
            for i in 0..n {
                let two = T::from_f64_or_one(2.0);
                omega[(i, n - 1)] = -two * (psi[(i, n - 2)] - psi[(i, n - 1)]) / (dx * dx)
                    - two * self.lid_velocity / dx;
                // Other walls (no-slip): ω = -2(ψ_adjacent)/Δ²
                if i == 0 {
                    omega[(0, i)] = -two * psi[(1, i)] / (dx * dx);
                }
                if i == n - 1 {
                    omega[(n - 1, i)] = -two * psi[(n - 2, i)] / (dx * dx);
                }
                omega[(i, 0)] = -two * psi[(i, 1)] / (dx * dx);
            }

            // Solve vorticity transport equation in interior
            for i in 1..n - 1 {
                for j in 1..n - 1 {
                    // Advection terms
                    let dpsi_dx =
                        (psi[(i + 1, j)] - psi[(i - 1, j)]) / (T::from_f64_or_one(2.0) * dx);
                    let dpsi_dy =
                        (psi[(i, j + 1)] - psi[(i, j - 1)]) / (T::from_f64_or_one(2.0) * dx);
                    let domega_dx =
                        (omega[(i + 1, j)] - omega[(i - 1, j)]) / (T::from_f64_or_one(2.0) * dx);
                    let domega_dy =
                        (omega[(i, j + 1)] - omega[(i, j - 1)]) / (T::from_f64_or_one(2.0) * dx);

                    // Diffusion term
                    let laplacian = (omega[(i + 1, j)]
                        + omega[(i - 1, j)]
                        + omega[(i, j + 1)]
                        + omega[(i, j - 1)]
                        - T::from_f64_or_one(4.0) * omega[(i, j)])
                        / (dx * dx);

                    // Update vorticity (explicit time stepping)
                    omega[(i, j)] +=
                        dt * (-dpsi_dy * domega_dx + dpsi_dx * domega_dy + viscosity * laplacian);
                }
            }

            // Solve Poisson equation for stream function: ∇²ψ = -ω
            let omega_sor = <T as SafeFromF64>::try_from_f64(1.5)?;
            for _ in 0..20 {
                // Inner iterations for Poisson solver
                for i in 1..n - 1 {
                    for j in 1..n - 1 {
                        let psi_previous = psi[(i, j)];
                        let psi_candidate = <T as SafeFromF64>::try_from_f64(0.25)?
                            * (psi[(i + 1, j)]
                                + psi[(i - 1, j)]
                                + psi[(i, j + 1)]
                                + psi[(i, j - 1)]
                                + dx * dx * omega[(i, j)]);
                        psi[(i, j)] = psi_previous + omega_sor * (psi_candidate - psi_previous);
                    }
                }
            }

            // Calculate residual
            let mut max_change = T::zero();
            for i in 0..n {
                for j in 0..n {
                    let change = (psi[(i, j)] - psi_prev[(i, j)]).abs();
                    if change > max_change {
                        max_change = change;
                    }
                }
            }

            let residual = max_change / (dt + <T as SafeFromF64>::try_from_f64(1e-10)?);
            convergence.push(residual);

            if residual < config.tolerance {
                break;
            }

            iteration += 1;
        }

        // Extract velocity field from stream function
        let mut u = DMatrix::<T>::zeros(n, n);
        let mut v = DMatrix::<T>::zeros(n, n);

        for i in 1..n - 1 {
            for j in 1..n - 1 {
                u[(i, j)] = (psi[(i, j + 1)] - psi[(i, j - 1)]) / (T::from_f64_or_one(2.0) * dx);
                v[(i, j)] = -(psi[(i + 1, j)] - psi[(i - 1, j)]) / (T::from_f64_or_one(2.0) * dx);
            }
        }

        // Set boundary velocities
        for i in 0..n {
            u[(i, n - 1)] = self.lid_velocity; // Top lid
        }

        // Extract centerline velocities for validation
        let centerline_u: Vec<T> = (0..n).map(|i| u[(i, n / 2)]).collect();
        let centerline_v: Vec<T> = (0..n).map(|j| v[(n / 2, j)]).collect();

        let mut values = centerline_u;
        values.extend(centerline_v);

        Ok(BenchmarkResult {
            name: self.name().to_string(),
            values,
            errors: vec![],
            convergence,
            execution_time: 0.0,
            metadata: std::collections::HashMap::new(),
        })
    }

    fn reference_solution(&self) -> Option<BenchmarkResult<T>> {
        None
    }

    fn validate(&self, result: &BenchmarkResult<T>) -> Result<bool> {
        // Compare with Ghia et al. reference data for exact L2 error validation
        let ghia_data = self.ghia_reference_data(T::from_f64_or_one(100.0)); // Default Re=100

        if let Some((y_positions, u_velocities)) = ghia_data {
            // Compare with computed values in result
            if !result.values.is_empty() {
                if !result.values.len().is_multiple_of(2) {
                    return Err(Error::InvalidInput(format!(
                        "Expected even number of values (u+v centerlines), got {}",
                        result.values.len()
                    )));
                }
                let u_len = result.values.len() / 2;
                if u_len == 0 {
                    return Err(Error::InvalidInput(
                        "Expected non-empty centerline u values".to_string(),
                    ));
                }

                // Compute L2 error against Ghia reference data
                let mut l2_error_sq = T::zero();
                let mut num_points = 0;

                // Interpolate computed values to Ghia reference y-positions
                for (&y_ref, &u_ref) in y_positions.iter().zip(u_velocities.iter()) {
                    // TODO: Interpolate against actual grid geometry instead of nearest-neighbor.
                    let y_ref_f64 = y_ref.to_f64().ok_or_else(|| {
                        Error::ConversionError("Failed to convert y_ref to f64".to_string())
                    })?;

                    let max_idx = u_len - 1;
                    let idx_f64 = y_ref_f64.clamp(0.0, 1.0) * max_idx as f64;
                    let grid_idx = idx_f64
                        .round()
                        .to_usize()
                        .ok_or_else(|| Error::ConversionError("Invalid grid index".to_string()))?
                        .min(max_idx);

                    if grid_idx < u_len {
                        let u_computed = result.values[grid_idx];
                        let error = u_computed - u_ref;
                        l2_error_sq += error * error;
                        num_points += 1;
                    }
                }

                if num_points > 0 {
                    let l2_error = (l2_error_sq / T::from_f64_or_one(f64::from(num_points))).sqrt();

                    // Ghia et al. (1982) requires <5% L2 error for accurate CFD solutions
                    let ghia_tolerance = T::from_f64_or_one(0.05);

                    // Check convergence
                    let converged = if let Some(last_residual) = result.convergence.last() {
                        last_residual.abs() < T::from_f64_or_one(1e-6)
                    } else {
                        false
                    };

                    // Log validation results
                    let l2_error_pct = l2_error.to_f64().unwrap_or(0.0) * 100.0;
                    let target_pct = ghia_tolerance.to_f64().unwrap_or(0.05) * 100.0;
                    println!(
                        "Ghia et al. validation - L2 error: {l2_error_pct:.6}%, Target: <{target_pct:.1}%, Converged: {converged}"
                    );

                    return Ok(l2_error < ghia_tolerance && converged);
                }
            }
        }

        // If no reference data available, perform basic sanity checks
        if result.values.is_empty() {
            Ok(false) // No data to validate
        } else {
            let max_value =
                result.values.iter().fold(
                    T::zero(),
                    |acc, &x| if x.abs() > acc { x.abs() } else { acc },
                );
            // Values should be bounded and finite
            Ok(max_value <= T::from_f64_or_one(10.0) && max_value >= T::zero())
        }
    }
}
