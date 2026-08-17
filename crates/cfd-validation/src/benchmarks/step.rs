//! Backward facing step benchmark problem
//!
//! Reference: Gartling (1990) "A test problem for outflow boundary conditions"

use super::{Benchmark, BenchmarkConfig, BenchmarkResult};
use crate::matrix::DMatrix;
use crate::scalar;
use cfd_core::error::{Error, Result};
use eunomia::{FloatElement, RealField};

const POISSON_SWEEPS_PER_VORTICITY_STEP: usize = 8;

/// Backward facing step benchmark
pub struct BackwardFacingStep<T: RealField + Copy> {
    /// Step height
    pub step_height: T,
    /// Channel height (after step)
    pub channel_height: T,
    /// Channel length
    pub channel_length: T,
    /// Inlet velocity
    pub inlet_velocity: T,
}

impl<T: RealField + Copy + FloatElement> BackwardFacingStep<T> {
    /// Create a new backward-facing-step benchmark.
    pub fn new(step_height: T, channel_height: T, channel_length: T, inlet_velocity: T) -> Self {
        Self {
            step_height,
            channel_height,
            channel_length,
            inlet_velocity,
        }
    }

    /// Return the inlet stream function at a vertical grid coordinate.
    fn inlet_streamfunction(&self, y: T, inlet_height: T, flow_rate: T) -> T {
        let zero = <T as FloatElement>::from_f64(0.0);
        if y <= self.step_height {
            return zero;
        }

        let eta = (y - self.step_height) / inlet_height;
        let two = <T as FloatElement>::from_f64(2.0);
        let three = <T as FloatElement>::from_f64(3.0);
        flow_rate * (three * eta * eta - two * eta * eta * eta)
    }

    /// Return the prescribed parabolic inlet velocity at a vertical coordinate.
    fn inlet_velocity_at(&self, y: T, inlet_height: T) -> T {
        let zero = <T as FloatElement>::from_f64(0.0);
        if y <= self.step_height {
            return zero;
        }

        let eta = (y - self.step_height) / inlet_height;
        let one = <T as FloatElement>::from_f64(1.0);
        let six = <T as FloatElement>::from_f64(6.0);
        six * self.inlet_velocity * eta * (one - eta)
    }

    /// Apply the Dirichlet stream-function boundaries for the step domain.
    ///
    /// The lower wall and the vertical step have `ψ = 0`. The upper wall has
    /// the inlet flow rate, the inlet has the integrated parabolic profile,
    /// and the outlet uses a zero stream-function gradient. This makes the
    /// solid step and the inlet/outlet contract explicit before each Poisson
    /// relaxation sweep.
    fn apply_streamfunction_boundaries(
        &self,
        psi: &mut [T],
        nx: usize,
        ny: usize,
        dy: T,
        inlet_height: T,
        flow_rate: T,
    ) {
        let zero = <T as FloatElement>::from_f64(0.0);
        let top_row = (ny - 1) * nx;

        for i in 0..nx {
            psi[i] = zero;
            psi[top_row + i] = flow_rate;
        }

        for j in 0..ny {
            let y = scalar::from_usize::<T>(j) * dy;
            psi[j * nx] = self.inlet_streamfunction(y, inlet_height, flow_rate);
        }

        for j in 1..ny - 1 {
            let row = j * nx;
            psi[row + nx - 1] = psi[row + nx - 2];
        }

        psi[top_row] = flow_rate;
    }

    /// Apply vorticity values implied by the no-slip and inlet boundaries.
    fn apply_vorticity_boundaries(
        &self,
        omega: &mut [T],
        psi: &[T],
        nx: usize,
        ny: usize,
        dx_squared: T,
        dy_squared: T,
        dy: T,
        flow_rate: T,
    ) {
        let zero = <T as FloatElement>::from_f64(0.0);
        let two = <T as FloatElement>::from_f64(2.0);
        let top = (ny - 1) * nx;

        for i in 1..nx {
            omega[i] = zero - two * psi[nx + i] / dy_squared;
            omega[top + i] = zero - two * (psi[top - nx + i] - flow_rate) / dy_squared;
        }

        for j in 1..ny - 1 {
            let row = j * nx;
            let y = scalar::from_usize::<T>(j) * dy;
            if y < self.step_height {
                omega[row] = zero - two * psi[row + 1] / dx_squared;
            } else {
                omega[row] = zero - (psi[row + nx] - two * psi[row] + psi[row - nx]) / dy_squared;
            }
        }

        for j in 1..ny - 1 {
            let row = j * nx;
            omega[row + nx - 1] = omega[row + nx - 2];
        }

        omega[top] = zero;
        omega[0] = zero;
    }

    /// Compute velocities from the divergence-free stream function.
    fn velocities_from_streamfunction(
        &self,
        psi: &[T],
        u: &mut [T],
        v: &mut [T],
        nx: usize,
        ny: usize,
        dx: T,
        dy: T,
        inlet_height: T,
    ) {
        let zero = <T as FloatElement>::from_f64(0.0);
        let two = <T as FloatElement>::from_f64(2.0);

        for j in 1..ny - 1 {
            let row = j * nx;
            for i in 1..nx - 1 {
                let index = row + i;
                u[index] = (psi[index + nx] - psi[index - nx]) / (two * dy);
                v[index] = zero - (psi[index + 1] - psi[index - 1]) / (two * dx);
            }
        }

        for i in 0..nx {
            u[i] = zero;
            v[i] = zero;
            u[(ny - 1) * nx + i] = zero;
            v[(ny - 1) * nx + i] = zero;
        }

        for j in 1..ny - 1 {
            let row = j * nx;
            let y = scalar::from_usize::<T>(j) * dy;
            u[row] = self.inlet_velocity_at(y, inlet_height);
            v[row] = zero;
            u[row + nx - 1] = u[row + nx - 2];
            v[row + nx - 1] = v[row + nx - 2];
        }

        u[0] = zero;
        v[0] = zero;
        u[(ny - 1) * nx] = zero;
        v[(ny - 1) * nx] = zero;
    }

    /// Calculate reattachment from the first downstream wall-shear sign change.
    ///
    /// The bottom wall satisfies no-slip, so the first off-wall velocity gives
    /// the one-sided shear-rate approximation `du/dy`. Recirculation has
    /// negative streamwise wall shear; reattachment is the first negative to
    /// non-negative crossing. Linear interpolation between adjacent grid
    /// values avoids snapping the result to a cell center.
    fn calculate_reattachment_length(u_field: &DMatrix<T>, dx: T, dy: T) -> Result<T> {
        let (ny, nx) = u_field.shape();
        if ny < 3 || nx < 3 {
            return Err(Error::InvalidInput(
                "backward-facing-step wall shear requires at least a 3x3 field".to_string(),
            ));
        }

        let zero = <T as FloatElement>::from_f64(0.0);
        let mut previous_shear = None;
        for i in 1..nx {
            let shear = (u_field[(1, i)] - u_field[(0, i)]) / dy;
            if !shear.is_finite() {
                return Err(Error::InvalidInput(
                    "backward-facing-step wall shear contains a non-finite value".to_string(),
                ));
            }

            if let Some(previous) = previous_shear {
                if previous < zero && shear >= zero {
                    let fraction = (zero - previous) / (shear - previous);
                    return Ok((scalar::from_usize::<T>(i - 1) + fraction) * dx);
                }
            }
            previous_shear = Some(shear);
        }

        Err(Error::InvalidInput(
            "backward-facing-step field has no downstream wall-shear reattachment crossing"
                .to_string(),
        ))
    }
}

impl<T: RealField + Copy + FloatElement> Benchmark<T> for BackwardFacingStep<T> {
    fn name(&self) -> &'static str {
        "Backward Facing Step"
    }

    fn description(&self) -> &'static str {
        "2D laminar flow over a backward-facing step"
    }

    fn run(&self, config: &BenchmarkConfig<T>) -> Result<BenchmarkResult<T>> {
        let zero = <T as FloatElement>::from_f64(0.0);
        if config.resolution < 4 || config.max_iterations == 0 {
            return Err(Error::InvalidConfiguration(
                "backward-facing-step resolution must be at least four and iterations positive"
                    .to_string(),
            ));
        }
        if !self.step_height.is_finite()
            || !self.channel_height.is_finite()
            || !self.channel_length.is_finite()
            || !self.inlet_velocity.is_finite()
            || !config.reynolds_number.is_finite()
            || self.step_height <= zero
            || self.channel_height <= self.step_height
            || self.channel_length <= zero
            || self.inlet_velocity <= zero
            || config.reynolds_number <= zero
        {
            return Err(Error::InvalidConfiguration(
                "backward-facing-step dimensions, inlet velocity, and Reynolds number must be finite and positive"
                    .to_string(),
            ));
        }

        let nx = config.resolution * 3; // Longer domain
        let ny = config.resolution;
        let inlet_height = self.channel_height - self.step_height;
        let flow_rate = self.inlet_velocity * inlet_height;

        // The streamfunction makes incompressibility exact on the discrete
        // grid; velocity is recovered from its curl after each relaxation.
        let mut psi = DMatrix::<T>::zeros(ny, nx);
        let mut omega = DMatrix::<T>::zeros(ny, nx);
        let mut next_omega = DMatrix::<T>::zeros(ny, nx);
        let mut u = DMatrix::<T>::zeros(ny, nx);
        let mut v = DMatrix::<T>::zeros(ny, nx);

        let dx = self.channel_length / scalar::from_usize::<T>(nx);
        let dy = self.channel_height / scalar::from_usize::<T>(ny);
        let nu = self.inlet_velocity * self.step_height / config.reynolds_number;
        let min_spacing = scalar::min(dx, dy);
        let diffusion_limit = <T as FloatElement>::from_f64(0.25)
            / (nu * (dx.recip() * dx.recip() + dy.recip() * dy.recip()));
        let advection_limit =
            <T as FloatElement>::from_f64(0.5) * min_spacing / self.inlet_velocity;
        let dt_limit = scalar::min(diffusion_limit, advection_limit);
        let dt = config
            .time_step
            .map_or(dt_limit, |requested| scalar::min(requested, dt_limit));
        let dx_squared = dx * dx;
        let dy_squared = dy * dy;

        let psi_data = psi
            .as_slice_mut()
            .expect("invariant: zero-initialized streamfunction is contiguous");
        self.apply_streamfunction_boundaries(psi_data, nx, ny, dy, inlet_height, flow_rate);
        let omega_data = omega
            .as_slice_mut()
            .expect("invariant: zero-initialized vorticity is contiguous");
        self.apply_vorticity_boundaries(
            omega_data, psi_data, nx, ny, dx_squared, dy_squared, dy, flow_rate,
        );

        let mut convergence = Vec::new();
        let two = <T as FloatElement>::from_f64(2.0);
        let poisson_denominator = two * (dx_squared + dy_squared);
        let poisson_source = dx_squared * dy_squared;

        for _ in 0..config.max_iterations {
            let psi_data = psi
                .as_slice_mut()
                .expect("invariant: zero-initialized streamfunction is contiguous");
            let omega_data = omega
                .as_slice()
                .expect("invariant: zero-initialized vorticity is contiguous");

            let mut poisson_residual = zero;
            for _ in 0..POISSON_SWEEPS_PER_VORTICITY_STEP {
                for j in 1..ny - 1 {
                    let row = j * nx;
                    for i in 1..nx - 1 {
                        let index = row + i;
                        let updated = (dy_squared * (psi_data[index + 1] + psi_data[index - 1])
                            + dx_squared * (psi_data[index + nx] + psi_data[index - nx])
                            + poisson_source * omega_data[index])
                            / poisson_denominator;
                        let change = scalar::abs(updated - psi_data[index]);
                        if change > poisson_residual {
                            poisson_residual = change;
                        }
                        psi_data[index] = updated;
                    }
                }
                self.apply_streamfunction_boundaries(psi_data, nx, ny, dy, inlet_height, flow_rate);
            }

            let u_data = u
                .as_slice_mut()
                .expect("invariant: zero-initialized velocity field is contiguous");
            let v_data = v
                .as_slice_mut()
                .expect("invariant: zero-initialized velocity field is contiguous");
            self.velocities_from_streamfunction(
                psi_data,
                u_data,
                v_data,
                nx,
                ny,
                dx,
                dy,
                inlet_height,
            );

            let next_omega_data = next_omega
                .as_slice_mut()
                .expect("invariant: zero-initialized vorticity workspace is contiguous");
            let mut vorticity_residual = poisson_residual;
            for j in 1..ny - 1 {
                let row = j * nx;
                for i in 1..nx - 1 {
                    let index = row + i;
                    let laplacian = (omega_data[index + 1] - two * omega_data[index]
                        + omega_data[index - 1])
                        / dx_squared
                        + (omega_data[index + nx] - two * omega_data[index]
                            + omega_data[index - nx])
                            / dy_squared;
                    let advection = u_data[index] * (omega_data[index + 1] - omega_data[index - 1])
                        / (two * dx)
                        + v_data[index] * (omega_data[index + nx] - omega_data[index - nx])
                            / (two * dy);
                    let updated = omega_data[index] + dt * (nu * laplacian - advection);
                    let change = scalar::abs(updated - omega_data[index]);
                    if change > vorticity_residual {
                        vorticity_residual = change;
                    }
                    next_omega_data[index] = updated;
                }
            }

            std::mem::swap(&mut omega, &mut next_omega);
            let omega_data = omega
                .as_slice_mut()
                .expect("invariant: zero-initialized vorticity is contiguous");
            let psi_data = psi
                .as_slice()
                .expect("invariant: zero-initialized streamfunction is contiguous");
            self.apply_vorticity_boundaries(
                omega_data, psi_data, nx, ny, dx_squared, dy_squared, dy, flow_rate,
            );
            convergence.push(vorticity_residual);

            if vorticity_residual < config.tolerance {
                break;
            }
        }

        let u_data = u
            .as_slice_mut()
            .expect("invariant: zero-initialized velocity field is contiguous");
        let v_data = v
            .as_slice_mut()
            .expect("invariant: zero-initialized velocity field is contiguous");
        let psi_data = psi
            .as_slice()
            .expect("invariant: zero-initialized streamfunction is contiguous");
        self.velocities_from_streamfunction(psi_data, u_data, v_data, nx, ny, dx, dy, inlet_height);
        let reattachment = Self::calculate_reattachment_length(&u, dx, dy)?;

        Ok(BenchmarkResult {
            name: self.name().to_string(),
            values: vec![reattachment],
            errors: vec![],
            convergence,
            execution_time: 0.0,
            metrics: std::collections::HashMap::new(),
            metadata: std::collections::HashMap::new(),
        })
    }

    fn reference_solution(&self) -> Option<BenchmarkResult<T>> {
        // Reference reattachment lengths from Gartling (1990)
        // "A test problem for outflow boundary conditions"
        // For expansion ratio ER = 2.0 (typical backward-facing step):
        // Re = 100:  x_r/h ≈ 3.0
        // Re = 200:  x_r/h ≈ 5.0
        // Re = 400:  x_r/h ≈ 6.8
        // Re = 800:  x_r/h ≈ 10.0
        //
        // For this implementation, we use the generic correlation
        // from experimental data: x_r/h ≈ 0.1 * Re^0.6 (Armaly et al. 1983)
        //
        // At Re=100: x_r/h ≈ 1.58 * h (conservative lower bound)
        // At Re=200: x_r/h ≈ 2.51 * h
        // At Re=400: x_r/h ≈ 3.98 * h
        // At Re=800: x_r/h ≈ 6.31 * h

        // Default reference for moderate Reynolds number (Re~200-400)
        let reference_reattachment = <T as FloatElement>::from_f64(6.0) * self.step_height;

        Some(BenchmarkResult {
            name: "Backward Facing Step (Reference)".to_string(),
            values: vec![reference_reattachment],
            errors: vec![],
            convergence: vec![],
            execution_time: 0.0,
            metrics: std::collections::HashMap::new(),
            metadata: std::collections::HashMap::new(),
        })
    }

    fn validate(&self, result: &BenchmarkResult<T>) -> Result<bool> {
        // Validate reattachment length against reference data
        if result.values.is_empty() {
            return Ok(false);
        }

        let computed_reattachment = result.values[0];

        // Get reference solution
        if let Some(reference) = self.reference_solution() {
            let reference_reattachment = reference.values[0];

            // Allow 30% tolerance due to:
            // 1. Different numerical schemes yield different reattachment points
            // 2. Grid resolution effects (coarse grids over-predict)
            // 3. Reference data variability across different studies
            // This is consistent with literature comparisons (Gartling 1990, Armaly et al. 1983)
            let tolerance = <T as FloatElement>::from_f64(0.30);
            let relative_error = scalar::abs(computed_reattachment - reference_reattachment)
                / reference_reattachment;

            // Validation passes if within 30% of reference
            let within_tolerance = relative_error <= tolerance;

            // Additional sanity checks
            let physically_reasonable = computed_reattachment > <T as FloatElement>::from_f64(0.0)
                && computed_reattachment < <T as FloatElement>::from_f64(20.0) * self.step_height;

            // Check convergence occurred
            let converged = if let Some(last_residual) = result.convergence.last() {
                scalar::abs(*last_residual) < <T as FloatElement>::from_f64(1e-4)
            } else {
                false
            };

            return Ok(within_tolerance && physically_reasonable && converged);
        }

        // Fallback: basic sanity checks without reference
        let physically_reasonable = computed_reattachment > <T as FloatElement>::from_f64(0.0)
            && computed_reattachment < <T as FloatElement>::from_f64(20.0) * self.step_height;

        Ok(physically_reasonable)
    }
}

#[cfg(test)]
mod tests {
    use super::{BackwardFacingStep, DMatrix};

    #[test]
    fn reattachment_interpolates_the_wall_shear_crossing() {
        let samples = [-1.0, -1.0, 1.0, 1.0, 1.0];
        let field = DMatrix::from_fn(4, samples.len(), |row, column| {
            if row == 1 {
                samples[column]
            } else {
                0.0
            }
        });

        let reattachment =
            BackwardFacingStep::<f64>::calculate_reattachment_length(&field, 1.0, 1.0)
                .expect("synthetic field has one valid crossing");
        assert!((reattachment - 1.5).abs() < f64::EPSILON);
    }

    #[test]
    fn reattachment_rejects_a_field_without_a_crossing() {
        let field = DMatrix::from_fn(4, 5, |row, _| if row == 1 { 1.0 } else { 0.0 });
        let error = BackwardFacingStep::<f64>::calculate_reattachment_length(&field, 1.0, 1.0)
            .expect_err("positive wall shear cannot identify reattachment");
        assert_eq!(
            error.to_string(),
            "Invalid input: backward-facing-step field has no downstream wall-shear reattachment crossing"
        );
    }

    #[test]
    fn streamfunction_boundaries_encode_step_and_channel_flow() {
        let benchmark = BackwardFacingStep::new(1.0_f64, 2.0, 10.0, 1.0);
        let mut psi = vec![0.0; 6 * 4];
        benchmark.apply_streamfunction_boundaries(&mut psi, 4, 6, 0.4, 1.0, 1.0);

        assert!(psi[5 * 4..].iter().all(|value| *value == 1.0));
        assert!(psi[0..4].iter().all(|value| *value == 0.0));
        assert_eq!(psi[4], 0.0);
        assert!(psi[12] > 0.0);
        assert_eq!(benchmark.inlet_velocity_at(0.8, 1.0), 0.0);
    }
}
