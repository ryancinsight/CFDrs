//! Fourier-space forcing and initialization utilities for periodic DNS.
//!
//! The helpers in this module generate seeded, band-limited, divergence-free
//! velocity fields in spectral space and convert them back to physical space
//! using Apollo FFTs. This is a practical first step toward seeded turbulence
//! initialization and stationary forcing.
//!
//! # Theorem - Solenoidal Projection
//!
//! For any nonzero wavevector `k`, the projection
//!
//! ```text
//! f_⊥ = f - k (k · f) / |k|²
//! ```
//!
//! removes the component parallel to `k`, so the resulting Fourier-space
//! forcing satisfies `k · f_⊥ = 0` and is divergence-free.

use apollofft::{ifft_3d_array, Complex64};
use cfd_core::error::{Error, Result};
use cfd_core::physics::fluid_dynamics::VelocityField;
use nalgebra::Vector3;
use ndarray::Array3;
use serde::{Deserialize, Serialize};

/// Seeded random-phase forcing configuration for periodic DNS.
#[derive(Debug, Clone, Copy, PartialEq, Serialize, Deserialize)]
pub struct BandLimitedRandomPhaseForcingConfig {
    /// Grid dimensions `(nx, ny, nz)`.
    pub dimensions: (usize, usize, usize),
    /// Physical domain lengths `(lx, ly, lz)`.
    pub lengths: (f64, f64, f64),
    /// Central shell index targeted by the forcing envelope.
    pub target_shell: usize,
    /// Half-width of the shell band around `target_shell`.
    pub shell_bandwidth: usize,
    /// Forcing amplitude.
    pub amplitude: f64,
    /// Deterministic seed used to synthesize pseudo-random phases.
    pub seed: u64,
}

impl BandLimitedRandomPhaseForcingConfig {
    /// Create a validated forcing configuration.
    pub fn new(
        dimensions: (usize, usize, usize),
        lengths: (f64, f64, f64),
        target_shell: usize,
        shell_bandwidth: usize,
        amplitude: f64,
        seed: u64,
    ) -> Result<Self> {
        let (nx, ny, nz) = dimensions;
        let (lx, ly, lz) = lengths;

        if nx == 0 || ny == 0 || nz == 0 {
            return Err(Error::InvalidConfiguration(
                "BandLimitedRandomPhaseForcingConfig: dimensions must be greater than zero".into(),
            ));
        }
        if lx <= 0.0 || ly <= 0.0 || lz <= 0.0 {
            return Err(Error::InvalidConfiguration(
                "BandLimitedRandomPhaseForcingConfig: domain lengths must be positive".into(),
            ));
        }
        if target_shell == 0 {
            return Err(Error::InvalidConfiguration(
                "BandLimitedRandomPhaseForcingConfig: target shell must be greater than zero"
                    .into(),
            ));
        }
        if amplitude < 0.0 || !amplitude.is_finite() {
            return Err(Error::InvalidConfiguration(
                "BandLimitedRandomPhaseForcingConfig: amplitude must be finite and non-negative"
                    .into(),
            ));
        }

        Ok(Self {
            dimensions,
            lengths,
            target_shell,
            shell_bandwidth,
            amplitude,
            seed,
        })
    }
}

/// Deterministic seeded forcing generator in Fourier space.
#[derive(Debug, Clone)]
pub struct BandLimitedRandomPhaseForcing3D {
    config: BandLimitedRandomPhaseForcingConfig,
}

impl BandLimitedRandomPhaseForcing3D {
    /// Create a new seeded band-limited forcing generator.
    pub fn new(config: BandLimitedRandomPhaseForcingConfig) -> Result<Self> {
        Ok(Self { config })
    }

    /// Return the current configuration.
    #[must_use]
    pub fn config(&self) -> BandLimitedRandomPhaseForcingConfig {
        self.config
    }

    /// Sample the forcing directly in Fourier space.
    pub fn sample_spectral_forcing(&self) -> Result<[Array3<Complex64>; 3]> {
        let (nx, ny, nz) = self.config.dimensions;
        let mut force_x = Array3::<Complex64>::zeros((nx, ny, nz));
        let mut force_y = Array3::<Complex64>::zeros((nx, ny, nz));
        let mut force_z = Array3::<Complex64>::zeros((nx, ny, nz));

        let lower_shell = self
            .config
            .target_shell
            .saturating_sub(self.config.shell_bandwidth);
        let upper_shell = self.config.target_shell + self.config.shell_bandwidth;

        for k in 0..nz {
            for j in 0..ny {
                for i in 0..nx {
                    let signed_kx = Self::signed_mode(i, nx);
                    let signed_ky = Self::signed_mode(j, ny);
                    let signed_kz = Self::signed_mode(k, nz);
                    let shell = Self::mode_shell(signed_kx, signed_ky, signed_kz);
                    if shell < lower_shell || shell > upper_shell || shell == 0 {
                        continue;
                    }

                    if !Self::is_canonical_mode(signed_kx, signed_ky, signed_kz) {
                        continue;
                    }

                    let envelope = self.shell_envelope(shell);
                    let mode = Vector3::new(signed_kx as f64, signed_ky as f64, signed_kz as f64);
                    let mut coeff = self.random_mode_vector(i, j, k, envelope);
                    self.project_divergence_free(&mode, &mut coeff);

                    if self.is_self_conjugate(i, j, k) {
                        coeff.x.im = 0.0;
                        coeff.y.im = 0.0;
                        coeff.z.im = 0.0;
                    }

                    force_x[[i, j, k]] = coeff.x;
                    force_y[[i, j, k]] = coeff.y;
                    force_z[[i, j, k]] = coeff.z;

                    let (mi, mj, mk) = self.mirror_index(i, j, k);
                    if (mi, mj, mk) != (i, j, k) {
                        force_x[[mi, mj, mk]] = coeff.x.conj();
                        force_y[[mi, mj, mk]] = coeff.y.conj();
                        force_z[[mi, mj, mk]] = coeff.z.conj();
                    }
                }
            }
        }

        Ok([force_x, force_y, force_z])
    }

    /// Sample the forcing in physical space.
    pub fn sample_physical_forcing(&self) -> Result<VelocityField<f64>> {
        let spectra = self.sample_spectral_forcing()?;
        let force_x = ifft_3d_array(&spectra[0]);
        let force_y = ifft_3d_array(&spectra[1]);
        let force_z = ifft_3d_array(&spectra[2]);
        Ok(self.arrays_to_velocity(force_x, force_y, force_z))
    }

    fn random_mode_vector(
        &self,
        i: usize,
        j: usize,
        k: usize,
        envelope: f64,
    ) -> Vector3<Complex64> {
        let amplitude = self.config.amplitude * envelope;
        let phase_x = self.seeded_phase(i, j, k, 0);
        let phase_y = self.seeded_phase(i, j, k, 1);
        let phase_z = self.seeded_phase(i, j, k, 2);

        Vector3::new(
            Complex64::from_polar(amplitude, phase_x),
            Complex64::from_polar(amplitude, phase_y),
            Complex64::from_polar(amplitude, phase_z),
        )
    }

    fn project_divergence_free(&self, mode: &Vector3<f64>, coeff: &mut Vector3<Complex64>) {
        let k2 = mode.dot(mode);
        if k2 <= 0.0 {
            return;
        }

        let k_dot_f = coeff.x * mode.x + coeff.y * mode.y + coeff.z * mode.z;
        let inv_k2 = 1.0 / k2;
        coeff.x -= k_dot_f * mode.x * inv_k2;
        coeff.y -= k_dot_f * mode.y * inv_k2;
        coeff.z -= k_dot_f * mode.z * inv_k2;
    }

    fn shell_envelope(&self, shell: usize) -> f64 {
        let distance = shell.abs_diff(self.config.target_shell) as f64;
        let width = self.config.shell_bandwidth.max(1) as f64;
        (-0.5 * (distance / width).powi(2)).exp()
    }

    fn seed_from_mode(&self, i: usize, j: usize, k: usize, component: u64) -> u64 {
        let mut state = self.config.seed
            ^ ((i as u64).wrapping_mul(0x9E37_79B9_7F4A_7C15))
            ^ ((j as u64).wrapping_mul(0xBF58_476D_1CE4_E5B9))
            ^ ((k as u64).wrapping_mul(0x94D0_49BB_1331_11EB))
            ^ component.rotate_left(17);
        state = state.wrapping_add(0x9E37_79B9_7F4A_7C15);
        state ^= state >> 30;
        state = state.wrapping_mul(0xBF58_476D_1CE4_E5B9);
        state ^= state >> 27;
        state = state.wrapping_mul(0x94D0_49BB_1331_11EB);
        state ^ (state >> 31)
    }

    fn seeded_phase(&self, i: usize, j: usize, k: usize, component: u64) -> f64 {
        let hash = self.seed_from_mode(i, j, k, component);
        let unit = (hash as f64) / (u64::MAX as f64);
        unit * std::f64::consts::TAU
    }

    fn signed_mode(index: usize, size: usize) -> isize {
        let half = size / 2;
        if index <= half {
            index as isize
        } else {
            index as isize - size as isize
        }
    }

    fn mirror_index(&self, i: usize, j: usize, k: usize) -> (usize, usize, usize) {
        let (nx, ny, nz) = self.config.dimensions;
        ((nx - i) % nx, (ny - j) % ny, (nz - k) % nz)
    }

    fn is_self_conjugate(&self, i: usize, j: usize, k: usize) -> bool {
        self.mirror_index(i, j, k) == (i, j, k)
    }

    fn is_canonical_mode(signed_kx: isize, signed_ky: isize, signed_kz: isize) -> bool {
        signed_kx > 0
            || (signed_kx == 0 && signed_ky > 0)
            || (signed_kx == 0 && signed_ky == 0 && signed_kz >= 0)
    }

    fn mode_shell(signed_kx: isize, signed_ky: isize, signed_kz: isize) -> usize {
        let radius =
            ((signed_kx * signed_kx + signed_ky * signed_ky + signed_kz * signed_kz) as f64).sqrt();
        radius.floor() as usize
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

/// Configuration for a time-resampled random-phase forcing schedule.
#[derive(Debug, Clone, Copy, PartialEq, Serialize, Deserialize)]
pub struct TimeResampledBandLimitedForcingConfig {
    /// Base band-limited forcing configuration.
    pub forcing: BandLimitedRandomPhaseForcingConfig,
    /// Number of DNS steps to reuse each forcing realization.
    pub resample_stride: usize,
}

impl TimeResampledBandLimitedForcingConfig {
    /// Create a validated time-resampled forcing configuration.
    pub fn new(
        forcing: BandLimitedRandomPhaseForcingConfig,
        resample_stride: usize,
    ) -> Result<Self> {
        if resample_stride == 0 {
            return Err(Error::InvalidConfiguration(
                "TimeResampledBandLimitedForcingConfig: resample stride must be greater than zero"
                    .into(),
            ));
        }

        Ok(Self {
            forcing,
            resample_stride,
        })
    }
}

/// Deterministic time-resampled forcing schedule for periodic DNS.
#[derive(Debug, Clone)]
pub struct TimeResampledBandLimitedForcing3D {
    config: TimeResampledBandLimitedForcingConfig,
}

impl TimeResampledBandLimitedForcing3D {
    /// Create a new time-resampled forcing schedule.
    pub fn new(config: TimeResampledBandLimitedForcingConfig) -> Result<Self> {
        Ok(Self { config })
    }

    /// Return the current schedule configuration.
    #[must_use]
    pub fn config(&self) -> TimeResampledBandLimitedForcingConfig {
        self.config
    }

    /// Sample the forcing field for a DNS step.
    pub fn sample_physical_forcing_for_step(&self, step: usize) -> Result<VelocityField<f64>> {
        self.forcing_for_step(step)?.sample_physical_forcing()
    }

    /// Sample the spectral forcing field for a DNS step.
    pub fn sample_spectral_forcing_for_step(&self, step: usize) -> Result<[Array3<Complex64>; 3]> {
        self.forcing_for_step(step)?.sample_spectral_forcing()
    }

    fn forcing_for_step(&self, step: usize) -> Result<BandLimitedRandomPhaseForcing3D> {
        let epoch = step / self.config.resample_stride;
        let seed = Self::mixed_seed(self.config.forcing.seed, epoch as u64);
        let forcing_config = BandLimitedRandomPhaseForcingConfig {
            seed,
            ..self.config.forcing
        };

        BandLimitedRandomPhaseForcing3D::new(forcing_config)
    }

    fn mixed_seed(base_seed: u64, epoch: u64) -> u64 {
        let mut state = base_seed ^ epoch.wrapping_mul(0x9E37_79B9_7F4A_7C15);
        state ^= state >> 30;
        state = state.wrapping_mul(0xBF58_476D_1CE4_E5B9);
        state ^= state >> 27;
        state = state.wrapping_mul(0x94D0_49BB_1331_11EB);
        state ^ (state >> 31)
    }
}

#[cfg(test)]
mod tests {
    use super::{
        BandLimitedRandomPhaseForcing3D, BandLimitedRandomPhaseForcingConfig,
        TimeResampledBandLimitedForcing3D, TimeResampledBandLimitedForcingConfig,
    };

    #[test]
    fn forcing_is_band_limited_and_divergence_free() {
        let config =
            BandLimitedRandomPhaseForcingConfig::new((8, 8, 8), (1.0, 1.0, 1.0), 2, 1, 0.5, 42)
                .expect("forcing config should be valid");
        let forcing = BandLimitedRandomPhaseForcing3D::new(config)
            .expect("forcing generator should be valid");
        let spectra = forcing
            .sample_spectral_forcing()
            .expect("spectral forcing should sample");

        let lower = config.target_shell.saturating_sub(config.shell_bandwidth);
        let upper = config.target_shell + config.shell_bandwidth;
        for k in 0..config.dimensions.2 {
            let signed_kz = if k <= config.dimensions.2 / 2 {
                k as isize
            } else {
                k as isize - config.dimensions.2 as isize
            };
            for j in 0..config.dimensions.1 {
                let signed_ky = if j <= config.dimensions.1 / 2 {
                    j as isize
                } else {
                    j as isize - config.dimensions.1 as isize
                };
                for i in 0..config.dimensions.0 {
                    let signed_kx = if i <= config.dimensions.0 / 2 {
                        i as isize
                    } else {
                        i as isize - config.dimensions.0 as isize
                    };
                    let shell = ((signed_kx * signed_kx
                        + signed_ky * signed_ky
                        + signed_kz * signed_kz) as f64)
                        .sqrt()
                        .floor() as usize;
                    let coeff_x = spectra[0][[i, j, k]];
                    let coeff_y = spectra[1][[i, j, k]];
                    let coeff_z = spectra[2][[i, j, k]];

                    if shell < lower || shell > upper || shell == 0 {
                        assert!(coeff_x.norm() < 1e-12);
                        assert!(coeff_y.norm() < 1e-12);
                        assert!(coeff_z.norm() < 1e-12);
                    } else {
                        let divergence = coeff_x * signed_kx as f64
                            + coeff_y * signed_ky as f64
                            + coeff_z * signed_kz as f64;
                        assert!(
                            divergence.norm() < 1e-10,
                            "mode ({i}, {j}, {k}) should be solenoidal"
                        );
                    }
                }
            }
        }
    }

    #[test]
    fn forcing_initialization_produces_real_nontrivial_velocity_field() {
        let config =
            BandLimitedRandomPhaseForcingConfig::new((6, 6, 6), (1.0, 1.0, 1.0), 1, 0, 0.25, 7)
                .expect("forcing config should be valid");
        let forcing = BandLimitedRandomPhaseForcing3D::new(config)
            .expect("forcing generator should be valid");
        let velocity = forcing
            .sample_physical_forcing()
            .expect("physical field should sample");

        let mut energy = 0.0;
        for sample in &velocity.components {
            assert!(sample.x.is_finite() && sample.y.is_finite() && sample.z.is_finite());
            energy += sample.norm_squared();
        }

        assert!(energy > 0.0, "forcing field should be nontrivial");
        assert_eq!(velocity.dimensions, (6, 6, 6));
    }

    #[test]
    fn time_resampled_schedule_reuses_epochs_and_changes_across_boundaries() {
        let forcing =
            BandLimitedRandomPhaseForcingConfig::new((8, 8, 8), (1.0, 1.0, 1.0), 2, 1, 0.5, 99)
                .expect("forcing config should be valid");
        let schedule = TimeResampledBandLimitedForcing3D::new(
            TimeResampledBandLimitedForcingConfig::new(forcing, 4)
                .expect("schedule config should be valid"),
        )
        .expect("schedule should be valid");

        let step_one = schedule
            .sample_physical_forcing_for_step(1)
            .expect("step-one forcing should sample");
        let step_three = schedule
            .sample_physical_forcing_for_step(3)
            .expect("step-three forcing should sample");
        let step_four = schedule
            .sample_physical_forcing_for_step(4)
            .expect("step-four forcing should sample");

        assert_eq!(step_one.components, step_three.components);
        assert_ne!(step_one.components, step_four.components);
    }
}
