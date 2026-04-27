//! Periodic pseudospectral DNS for incompressible 3D Navier-Stokes.
//!
//! This module provides a small Apollo-backed stepping kernel for periodic
//! domains. It evaluates the nonlinear convective term in physical space,
//! applies a 2/3 de-aliasing filter in Fourier space, and advances the flow
//! with a projected semi-implicit Euler step. An optional forcing field can
//! be injected through the same Fourier-space projection path.
//!
//! # Theorem - 2/3 Dealiasing for Quadratic Nonlinearity
//!
//! If all Fourier modes with index magnitude greater than `N/3` are removed
//! before forming quadratic products, then the resulting pseudospectral
//! evaluation is alias-free for a quadratic nonlinearity on an `N`-point grid.
//! This is the standard Orszag 2/3 rule used here to keep the nonlinear term
//! free of wraparound contamination.

use apollofft::{Complex64, FftPlan3D, Shape3D};
use cfd_core::error::{Error, Result};
use cfd_core::physics::fluid_dynamics::VelocityField;
use nalgebra::Vector3;
use ndarray::Array3;
use std::sync::Arc;

/// Configuration for a periodic pseudospectral DNS stepper.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct PeriodicPseudospectralDnsConfig {
    /// Grid dimensions `(nx, ny, nz)`.
    pub dimensions: (usize, usize, usize),
    /// Physical domain lengths `(lx, ly, lz)`.
    pub lengths: (f64, f64, f64),
    /// Kinematic viscosity `nu`.
    pub kinematic_viscosity: f64,
    /// Time step `dt`.
    pub time_step: f64,
}

impl PeriodicPseudospectralDnsConfig {
    /// Create a validated configuration.
    pub fn new(
        dimensions: (usize, usize, usize),
        lengths: (f64, f64, f64),
        kinematic_viscosity: f64,
        time_step: f64,
    ) -> Result<Self> {
        let (nx, ny, nz) = dimensions;
        let (lx, ly, lz) = lengths;

        if nx == 0 || ny == 0 || nz == 0 {
            return Err(Error::InvalidConfiguration(
                "PeriodicPseudospectralDnsConfig: dimensions must be greater than zero".into(),
            ));
        }
        if lx <= 0.0 || ly <= 0.0 || lz <= 0.0 {
            return Err(Error::InvalidConfiguration(
                "PeriodicPseudospectralDnsConfig: domain lengths must be positive".into(),
            ));
        }
        if !kinematic_viscosity.is_finite() || kinematic_viscosity < 0.0 {
            return Err(Error::InvalidConfiguration(
                "PeriodicPseudospectralDnsConfig: kinematic viscosity must be finite and non-negative".into(),
            ));
        }
        if !time_step.is_finite() || time_step <= 0.0 {
            return Err(Error::InvalidConfiguration(
                "PeriodicPseudospectralDnsConfig: time step must be finite and positive".into(),
            ));
        }

        Ok(Self {
            dimensions,
            lengths,
            kinematic_viscosity,
            time_step,
        })
    }
}

/// Apollo-backed periodic pseudospectral Navier-Stokes stepper.
#[derive(Debug, Clone)]
pub struct PeriodicPseudospectralDns3D {
    config: PeriodicPseudospectralDnsConfig,
    fft_plan: Arc<FftPlan3D>,
    wavenumbers_x: Vec<f64>,
    wavenumbers_y: Vec<f64>,
    wavenumbers_z: Vec<f64>,
    cutoff_x: usize,
    cutoff_y: usize,
    cutoff_z: usize,
}

impl PeriodicPseudospectralDns3D {
    /// Create a new periodic DNS stepper.
    pub fn new(config: PeriodicPseudospectralDnsConfig) -> Result<Self> {
        let (nx, ny, nz) = config.dimensions;
        let (lx, ly, lz) = config.lengths;
        let shape = Shape3D::new(nx, ny, nz).map_err(|error| {
            Error::InvalidConfiguration(format!(
                "PeriodicPseudospectralDns3D: invalid Apollo FFT shape: {error}"
            ))
        })?;

        Ok(Self {
            fft_plan: Arc::new(FftPlan3D::new(shape)),
            wavenumbers_x: Self::wavenumber_axis(nx, lx),
            wavenumbers_y: Self::wavenumber_axis(ny, ly),
            wavenumbers_z: Self::wavenumber_axis(nz, lz),
            cutoff_x: nx / 3,
            cutoff_y: ny / 3,
            cutoff_z: nz / 3,
            config,
        })
    }

    /// Return the current configuration.
    #[must_use]
    pub fn config(&self) -> PeriodicPseudospectralDnsConfig {
        self.config
    }

    /// Compute the nonlinear advection term `u · ∇u` in physical space.
    pub fn nonlinear_advection(&self, velocity: &VelocityField<f64>) -> Result<VelocityField<f64>> {
        let spectra = self.velocity_spectrum(velocity)?;
        self.nonlinear_advection_from_spectrum(&spectra)
    }

    /// Advance one explicit semi-implicit Euler step.
    ///
    /// The nonlinear term is evaluated pseudospectrally and projected back to a
    /// divergence-free subspace before the viscous update.
    pub fn advance_euler(&self, velocity: &VelocityField<f64>) -> Result<VelocityField<f64>> {
        let spectra = self.velocity_spectrum(velocity)?;
        let nonlinear_hat = self.nonlinear_advection_hat(&spectra)?;
        self.advance_from_spectra(&spectra, &nonlinear_hat, None)
    }

    /// Advance one step with an externally supplied physical-space forcing field.
    pub fn advance_euler_with_forcing(
        &self,
        velocity: &VelocityField<f64>,
        forcing: &VelocityField<f64>,
    ) -> Result<VelocityField<f64>> {
        let spectra = self.velocity_spectrum(velocity)?;
        let nonlinear_hat = self.nonlinear_advection_hat(&spectra)?;
        let forcing_hat = self.velocity_spectrum(forcing)?;
        self.advance_from_spectra(&spectra, &nonlinear_hat, Some(&forcing_hat))
    }

    fn velocity_spectrum(&self, velocity: &VelocityField<f64>) -> Result<[Array3<Complex64>; 3]> {
        let (nx, ny, nz) = self.config.dimensions;
        if velocity.dimensions != self.config.dimensions {
            return Err(Error::InvalidConfiguration(format!(
                "PeriodicPseudospectralDns3D: expected dimensions ({nx}, {ny}, {nz}), found ({}, {}, {})",
                velocity.dimensions.0, velocity.dimensions.1, velocity.dimensions.2
            )));
        }
        if velocity.components.len() != nx * ny * nz {
            return Err(Error::InvalidConfiguration(format!(
                "PeriodicPseudospectralDns3D: expected {} samples, found {}",
                nx * ny * nz,
                velocity.components.len()
            )));
        }

        let (u, v, w) = self.velocity_components_to_arrays(velocity);
        let mut u_hat = self.fft_plan.forward_real_to_complex(&u);
        let mut v_hat = self.fft_plan.forward_real_to_complex(&v);
        let mut w_hat = self.fft_plan.forward_real_to_complex(&w);
        self.apply_dealiasing_filter(&mut u_hat);
        self.apply_dealiasing_filter(&mut v_hat);
        self.apply_dealiasing_filter(&mut w_hat);
        Ok([u_hat, v_hat, w_hat])
    }

    fn nonlinear_advection_from_spectrum(
        &self,
        spectra: &[Array3<Complex64>; 3],
    ) -> Result<VelocityField<f64>> {
        let [u_hat, v_hat, w_hat] = spectra;
        let (u, v, w) = (
            self.fft_plan.inverse_complex_to_real(u_hat),
            self.fft_plan.inverse_complex_to_real(v_hat),
            self.fft_plan.inverse_complex_to_real(w_hat),
        );

        let du_dx = self.spectral_derivative(u_hat, 0);
        let du_dy = self.spectral_derivative(u_hat, 1);
        let du_dz = self.spectral_derivative(u_hat, 2);
        let dv_dx = self.spectral_derivative(v_hat, 0);
        let dv_dy = self.spectral_derivative(v_hat, 1);
        let dv_dz = self.spectral_derivative(v_hat, 2);
        let dw_dx = self.spectral_derivative(w_hat, 0);
        let dw_dy = self.spectral_derivative(w_hat, 1);
        let dw_dz = self.spectral_derivative(w_hat, 2);

        let mut adv_x = Array3::<f64>::zeros(self.config.dimensions);
        let mut adv_y = Array3::<f64>::zeros(self.config.dimensions);
        let mut adv_z = Array3::<f64>::zeros(self.config.dimensions);

        for k in 0..self.config.dimensions.2 {
            for j in 0..self.config.dimensions.1 {
                for i in 0..self.config.dimensions.0 {
                    let velocity_x = u[[i, j, k]];
                    let velocity_y = v[[i, j, k]];
                    let velocity_z = w[[i, j, k]];
                    adv_x[[i, j, k]] = velocity_x * du_dx[[i, j, k]]
                        + velocity_y * du_dy[[i, j, k]]
                        + velocity_z * du_dz[[i, j, k]];
                    adv_y[[i, j, k]] = velocity_x * dv_dx[[i, j, k]]
                        + velocity_y * dv_dy[[i, j, k]]
                        + velocity_z * dv_dz[[i, j, k]];
                    adv_z[[i, j, k]] = velocity_x * dw_dx[[i, j, k]]
                        + velocity_y * dw_dy[[i, j, k]]
                        + velocity_z * dw_dz[[i, j, k]];
                }
            }
        }

        Ok(self.arrays_to_velocity(adv_x, adv_y, adv_z))
    }

    fn nonlinear_advection_hat(
        &self,
        spectra: &[Array3<Complex64>; 3],
    ) -> Result<[Array3<Complex64>; 3]> {
        let adv = self.nonlinear_advection_from_spectrum(spectra)?;
        let (adv_x, adv_y, adv_z) = self.velocity_components_to_arrays(&adv);
        let mut adv_hat_x = self.fft_plan.forward_real_to_complex(&adv_x);
        let mut adv_hat_y = self.fft_plan.forward_real_to_complex(&adv_y);
        let mut adv_hat_z = self.fft_plan.forward_real_to_complex(&adv_z);
        self.apply_dealiasing_filter(&mut adv_hat_x);
        self.apply_dealiasing_filter(&mut adv_hat_y);
        self.apply_dealiasing_filter(&mut adv_hat_z);
        Ok([adv_hat_x, adv_hat_y, adv_hat_z])
    }

    fn advance_from_spectra(
        &self,
        spectra: &[Array3<Complex64>; 3],
        nonlinear_hat: &[Array3<Complex64>; 3],
        forcing_hat: Option<&[Array3<Complex64>; 3]>,
    ) -> Result<VelocityField<f64>> {
        let (nx, ny, nz) = self.config.dimensions;
        let mut next_u_hat = Array3::<Complex64>::zeros((nx, ny, nz));
        let mut next_v_hat = Array3::<Complex64>::zeros((nx, ny, nz));
        let mut next_w_hat = Array3::<Complex64>::zeros((nx, ny, nz));

        let [u_hat, v_hat, w_hat] = spectra;
        let [adv_hat_x, adv_hat_y, adv_hat_z] = nonlinear_hat;

        for k in 0..nz {
            for j in 0..ny {
                for i in 0..nx {
                    let kx = self.wavenumbers_x[i];
                    let ky = self.wavenumbers_y[j];
                    let kz = self.wavenumbers_z[k];
                    let k2 = kx * kx + ky * ky + kz * kz;
                    let dt = self.config.time_step;
                    let nu = self.config.kinematic_viscosity;

                    let mut rhs_x = -adv_hat_x[[i, j, k]];
                    let mut rhs_y = -adv_hat_y[[i, j, k]];
                    let mut rhs_z = -adv_hat_z[[i, j, k]];

                    if let Some(forcing_hat) = forcing_hat {
                        rhs_x += forcing_hat[0][[i, j, k]];
                        rhs_y += forcing_hat[1][[i, j, k]];
                        rhs_z += forcing_hat[2][[i, j, k]];
                    }

                    if k2 > 0.0 {
                        let dot = rhs_x * kx + rhs_y * ky + rhs_z * kz;
                        let k2_inv = 1.0 / k2;
                        rhs_x -= dot * kx * k2_inv;
                        rhs_y -= dot * ky * k2_inv;
                        rhs_z -= dot * kz * k2_inv;
                    } else {
                        rhs_x = Complex64::new(0.0, 0.0);
                        rhs_y = Complex64::new(0.0, 0.0);
                        rhs_z = Complex64::new(0.0, 0.0);
                    }

                    let viscous_denominator = 1.0 + nu * k2 * dt;
                    next_u_hat[[i, j, k]] = (u_hat[[i, j, k]] + dt * rhs_x) / viscous_denominator;
                    next_v_hat[[i, j, k]] = (v_hat[[i, j, k]] + dt * rhs_y) / viscous_denominator;
                    next_w_hat[[i, j, k]] = (w_hat[[i, j, k]] + dt * rhs_z) / viscous_denominator;
                }
            }
        }

        self.apply_dealiasing_filter(&mut next_u_hat);
        self.apply_dealiasing_filter(&mut next_v_hat);
        self.apply_dealiasing_filter(&mut next_w_hat);

        let next_u = self.fft_plan.inverse_complex_to_real(&next_u_hat);
        let next_v = self.fft_plan.inverse_complex_to_real(&next_v_hat);
        let next_w = self.fft_plan.inverse_complex_to_real(&next_w_hat);
        Ok(self.arrays_to_velocity(next_u, next_v, next_w))
    }

    fn spectral_derivative(&self, spectrum: &Array3<Complex64>, axis: usize) -> Array3<f64> {
        let mut derivative_hat = spectrum.clone();
        let (nx, ny, nz) = self.config.dimensions;
        for k in 0..nz {
            for j in 0..ny {
                for i in 0..nx {
                    let wave_number = match axis {
                        0 => self.wavenumbers_x[i],
                        1 => self.wavenumbers_y[j],
                        _ => self.wavenumbers_z[k],
                    };
                    derivative_hat[[i, j, k]] *= Complex64::new(0.0, wave_number);
                }
            }
        }
        self.fft_plan.inverse_complex_to_real(&derivative_hat)
    }

    fn apply_dealiasing_filter(&self, spectrum: &mut Array3<Complex64>) {
        let (nx, ny, nz) = self.config.dimensions;
        for k in 0..nz {
            let keep_z = self.mode_is_kept(k, nz, self.cutoff_z);
            for j in 0..ny {
                let keep_y = self.mode_is_kept(j, ny, self.cutoff_y);
                for i in 0..nx {
                    if !(self.mode_is_kept(i, nx, self.cutoff_x) && keep_y && keep_z) {
                        spectrum[[i, j, k]] = Complex64::new(0.0, 0.0);
                    }
                }
            }
        }
    }

    fn mode_is_kept(&self, index: usize, size: usize, cutoff: usize) -> bool {
        self.signed_mode(index, size).unsigned_abs() <= cutoff
    }

    fn signed_mode(&self, index: usize, size: usize) -> isize {
        let half = size / 2;
        if index <= half {
            index as isize
        } else {
            index as isize - size as isize
        }
    }

    fn wavenumber_axis(size: usize, length: f64) -> Vec<f64> {
        let factor = 2.0 * std::f64::consts::PI / length;
        (0..size)
            .map(|index| Self::signed_mode_static(index, size) as f64 * factor)
            .collect()
    }

    fn signed_mode_static(index: usize, size: usize) -> isize {
        let half = size / 2;
        if index <= half {
            index as isize
        } else {
            index as isize - size as isize
        }
    }

    fn velocity_components_to_arrays(
        &self,
        velocity: &VelocityField<f64>,
    ) -> (Array3<f64>, Array3<f64>, Array3<f64>) {
        let (nx, ny, nz) = self.config.dimensions;
        let mut u = Array3::<f64>::zeros((nx, ny, nz));
        let mut v = Array3::<f64>::zeros((nx, ny, nz));
        let mut w = Array3::<f64>::zeros((nx, ny, nz));

        for k in 0..nz {
            for j in 0..ny {
                for i in 0..nx {
                    let idx = k * nx * ny + j * nx + i;
                    let sample = velocity.components[idx];
                    u[[i, j, k]] = sample.x;
                    v[[i, j, k]] = sample.y;
                    w[[i, j, k]] = sample.z;
                }
            }
        }

        (u, v, w)
    }

    fn arrays_to_velocity(
        &self,
        u: Array3<f64>,
        v: Array3<f64>,
        w: Array3<f64>,
    ) -> VelocityField<f64> {
        let (nx, ny, nz) = self.config.dimensions;
        let mut components = Vec::with_capacity(nx * ny * nz);

        for k in 0..nz {
            for j in 0..ny {
                for i in 0..nx {
                    components.push(Vector3::new(u[[i, j, k]], v[[i, j, k]], w[[i, j, k]]));
                }
            }
        }

        VelocityField {
            components,
            dimensions: self.config.dimensions,
        }
    }
}

#[cfg(test)]
mod tests {
    use super::{PeriodicPseudospectralDns3D, PeriodicPseudospectralDnsConfig};
    use crate::spectral::{BandLimitedRandomPhaseForcing3D, BandLimitedRandomPhaseForcingConfig};
    use apollofft::Complex64;
    use cfd_core::physics::fluid_dynamics::VelocityField;
    use nalgebra::Vector3;

    fn constant_velocity_field(
        dimensions: (usize, usize, usize),
        velocity: Vector3<f64>,
    ) -> VelocityField<f64> {
        let total_points = dimensions.0 * dimensions.1 * dimensions.2;
        VelocityField {
            components: vec![velocity; total_points],
            dimensions,
        }
    }

    #[test]
    fn constant_field_is_a_fixed_point() {
        let config = PeriodicPseudospectralDnsConfig::new((4, 4, 4), (1.0, 1.0, 1.0), 0.01, 0.05)
            .expect("config should be valid");
        let stepper = PeriodicPseudospectralDns3D::new(config).expect("stepper should be valid");
        let velocity = constant_velocity_field((4, 4, 4), Vector3::new(1.0, -0.5, 0.25));

        let nonlinear = stepper
            .nonlinear_advection(&velocity)
            .expect("nonlinear term should compute");
        let advanced = stepper
            .advance_euler(&velocity)
            .expect("advance should succeed");

        for sample in nonlinear.components {
            assert!(sample.x.abs() < 1e-12);
            assert!(sample.y.abs() < 1e-12);
            assert!(sample.z.abs() < 1e-12);
        }

        for sample in advanced.components {
            assert!((sample.x - 1.0).abs() < 1e-12);
            assert!((sample.y + 0.5).abs() < 1e-12);
            assert!((sample.z - 0.25).abs() < 1e-12);
        }
    }

    #[test]
    fn dealiasing_filter_removes_high_modes() {
        let config = PeriodicPseudospectralDnsConfig::new((6, 6, 6), (1.0, 1.0, 1.0), 0.0, 0.05)
            .expect("config should be valid");
        let stepper = PeriodicPseudospectralDns3D::new(config).expect("stepper should be valid");
        let mut spectrum = ndarray::Array3::from_elem((6, 6, 6), Complex64::new(1.0, 0.0));

        stepper.apply_dealiasing_filter(&mut spectrum);

        let cutoff = 6 / 3;
        for k in 0..6 {
            let kz = if k <= 3 { k as isize } else { k as isize - 6 };
            for j in 0..6 {
                let ky = if j <= 3 { j as isize } else { j as isize - 6 };
                for i in 0..6 {
                    let kx = if i <= 3 { i as isize } else { i as isize - 6 };
                    let kept = kx.abs() as usize <= cutoff
                        && ky.abs() as usize <= cutoff
                        && kz.abs() as usize <= cutoff;
                    let value = spectrum[[i, j, k]];
                    if kept {
                        assert!((value.re - 1.0).abs() < 1e-12);
                    } else {
                        assert!(
                            value.norm() < 1e-12,
                            "mode ({i}, {j}, {k}) should be filtered"
                        );
                    }
                }
            }
        }
    }

    #[test]
    fn forcing_drives_nontrivial_update() {
        let dns_config = PeriodicPseudospectralDnsConfig::new((4, 4, 4), (1.0, 1.0, 1.0), 0.0, 0.1)
            .expect("config should be valid");
        let stepper =
            PeriodicPseudospectralDns3D::new(dns_config).expect("stepper should be valid");
        let forcing_config =
            BandLimitedRandomPhaseForcingConfig::new((4, 4, 4), (1.0, 1.0, 1.0), 1, 0, 0.5, 123)
                .expect("forcing config should be valid");
        let forcing = BandLimitedRandomPhaseForcing3D::new(forcing_config)
            .expect("forcing generator should be valid");
        let forcing_field = forcing
            .sample_physical_forcing()
            .expect("forcing field should sample");
        let velocity = constant_velocity_field((4, 4, 4), Vector3::zeros());

        let advanced = stepper
            .advance_euler_with_forcing(&velocity, &forcing_field)
            .expect("forced advance should succeed");

        assert!(
            advanced
                .components
                .iter()
                .any(|sample| sample.norm() > 1e-10),
            "forcing must drive a nontrivial update"
        );
    }
}
