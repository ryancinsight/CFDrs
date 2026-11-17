//! Spectral diagnostics for velocity fields.
//!
//! This module provides FFT-backed post-processing helpers that are useful for
//! DNS validation and turbulence analysis, including isotropic kinetic energy
//! spectra, enstrophy spectra, probe-signal power spectra, and normalized
//! temporal autocorrelation series.

mod signals;
mod spectra;

pub use signals::{ProbeSignalSpectrum, TemporalAutocorrelation, probe_signal_spectrum, temporal_autocorrelation};
pub use spectra::{EnstrophySpectrum, KineticEnergySpectrum, enstrophy_spectrum, kinetic_energy_spectrum};

#[cfg(test)]
mod tests {
    use super::*;
    use cfd_core::physics::fluid_dynamics::VelocityField;
    use leto::geometry::Vector3;

    #[test]
    fn constant_velocity_energy_lives_in_zero_shell() {
        let mut velocity = VelocityField {
            components: vec![Vector3::new(2.0, -1.0, 0.5); 4 * 4 * 4],
            dimensions: (4, 4, 4),
        };

        let spectrum = kinetic_energy_spectrum(&velocity).expect("spectrum should compute");
        let expected_energy = 0.5 * (2.0_f64 * 2.0 + 1.0_f64 + 0.5_f64 * 0.5) * 64.0;

        assert!((spectrum.total_energy - expected_energy).abs() < 1e-10);
        assert!((spectrum.shell_energy[0] - expected_energy).abs() < 1e-10);
        assert!(
            spectrum
                .shell_energy
                .iter()
                .skip(1)
                .all(|energy| energy.abs() < 1e-10)
        );

        // Ensure the input is left untouched.
        velocity.components[0] = Vector3::new(0.0, 0.0, 0.0);
        assert_eq!(velocity.components.len(), 64);
    }

    #[test]
    fn single_mode_energy_peaks_in_first_shell() {
        let (nx, ny, nz) = (8, 4, 4);
        let mut velocity = VelocityField {
            components: Vec::with_capacity(nx * ny * nz),
            dimensions: (nx, ny, nz),
        };

        for _k in 0..nz {
            for _j in 0..ny {
                for i in 0..nx {
                    let phase = 2.0 * std::f64::consts::PI * i as f64 / nx as f64;
                    velocity
                        .components
                        .push(Vector3::new(phase.sin(), 0.0, 0.0));
                }
            }
        }

        let spectrum = kinetic_energy_spectrum(&velocity).expect("spectrum should compute");
        let expected_energy = 0.25 * (nx * ny * nz) as f64;
        let dominant_shell = spectrum
            .shell_energy
            .iter()
            .enumerate()
            .max_by(|lhs, rhs| {
                lhs.1
                    .partial_cmp(rhs.1)
                    .expect("energy should be comparable")
            })
            .map(|(index, _)| index)
            .expect("spectrum should contain at least one shell");

        assert!((spectrum.total_energy - expected_energy).abs() < 1e-10);
        assert!(spectrum.shell_energy[0] < 1e-10);
        assert_eq!(dominant_shell, 1);
        assert!(spectrum.shell_energy[1] > 0.99 * spectrum.total_energy);
    }

    #[test]
    fn constant_velocity_has_zero_enstrophy() {
        let velocity = VelocityField {
            components: vec![Vector3::new(1.0, 2.0, -1.5); 4 * 4 * 4],
            dimensions: (4, 4, 4),
        };

        let spectrum = enstrophy_spectrum(&velocity, (1.0, 1.0, 1.0))
            .expect("enstrophy spectrum should compute");

        assert!(spectrum.total_enstrophy.abs() < 1e-10);
        assert!(
            spectrum
                .shell_enstrophy
                .iter()
                .all(|value| value.abs() < 1e-10)
        );
    }

    #[test]
    fn single_mode_enstrophy_peaks_in_first_shell() {
        let (nx, ny, nz) = (8, 8, 8);
        let mut velocity = VelocityField {
            components: Vec::with_capacity(nx * ny * nz),
            dimensions: (nx, ny, nz),
        };

        for _k in 0..nz {
            for j in 0..ny {
                let phase = 2.0 * std::f64::consts::PI * j as f64 / ny as f64;
                for _i in 0..nx {
                    velocity
                        .components
                        .push(Vector3::new(phase.sin(), 0.0, 0.0));
                }
            }
        }

        let spectrum = enstrophy_spectrum(&velocity, (1.0, 1.0, 1.0))
            .expect("enstrophy spectrum should compute");

        let dominant_shell = spectrum
            .shell_enstrophy
            .iter()
            .enumerate()
            .max_by(|lhs, rhs| {
                lhs.1
                    .partial_cmp(rhs.1)
                    .expect("enstrophy should be comparable")
            })
            .map(|(index, _)| index)
            .expect("spectrum should contain at least one shell");

        assert!(spectrum.total_enstrophy > 0.0);
        assert_eq!(dominant_shell, 1);
        assert!(spectrum.shell_enstrophy[1] > 0.99 * spectrum.total_enstrophy);
    }

    #[test]
    fn alternating_signal_has_predictable_autocorrelation_and_probe_peak() {
        let samples: Vec<f64> = (0..16)
            .map(|index| if index % 2 == 0 { 1.0 } else { -1.0 })
            .collect();

        let autocorrelation =
            temporal_autocorrelation(&samples, 0.25, 4).expect("autocorrelation should compute");
        assert_eq!(autocorrelation.lags, vec![0, 1, 2, 3, 4]);
        assert!((autocorrelation.values[0] - 1.0).abs() < 1e-10);
        assert!((autocorrelation.values[1] + 1.0).abs() < 1e-10);
        assert!((autocorrelation.values[2] - 1.0).abs() < 1e-10);

        let spectrum =
            probe_signal_spectrum(&samples, 0.25).expect("probe signal spectrum should compute");
        let dominant_bin = spectrum
            .spectral_energy
            .iter()
            .enumerate()
            .max_by(|lhs, rhs| {
                lhs.1
                    .partial_cmp(rhs.1)
                    .expect("power should be comparable")
            })
            .map(|(index, _)| index)
            .expect("spectrum should contain at least one bin");

        assert_eq!(dominant_bin, 8);
        assert!((spectrum.frequencies_hz[dominant_bin] - 2.0).abs() < 1e-10);
        assert!(spectrum.spectral_energy[dominant_bin] > 0.99 * spectrum.total_fluctuation_energy);
        assert!((spectrum.variance - 1.0).abs() < 1e-10);
    }
}
