//! 3D stationary forced isotropic turbulence benchmark.
//!
//! This benchmark exercises the Apollo-backed periodic DNS stepper with a
//! deterministic, time-resampled band-limited forcing schedule. It starts from
//! quiescent flow, drives energy into a small shell band, and records the
//! kinetic-energy spectrum over time so the validation layer can inspect both
//! the forcing cadence and the spectral response.

use super::super::{Benchmark, BenchmarkConfig, BenchmarkResult};
use cfd_3d::spectral::{
    enstrophy_spectrum, kinetic_energy_spectrum, probe_signal_spectrum,
    temporal_autocorrelation, BandLimitedRandomPhaseForcingConfig, EnstrophySpectrum,
    KineticEnergySpectrum, PeriodicPseudospectralDns3D, PeriodicPseudospectralDnsConfig,
    TimeResampledBandLimitedForcing3D, TimeResampledBandLimitedForcingConfig,
};
use cfd_core::error::{Error, Result};
use cfd_core::physics::fluid_dynamics::VelocityField;
use nalgebra::Vector3;
use serde::{Deserialize, Serialize};
use std::time::Instant;

/// Configuration for the stationary forced turbulence benchmark.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct ForcedTurbulenceBenchmarkConfig {
    /// Grid dimensions `(nx, ny, nz)`.
    pub dimensions: (usize, usize, usize),
    /// Characteristic length scale `L` used to build the periodic domain.
    pub length_scale: f64,
    /// Kinematic viscosity `nu`.
    pub kinematic_viscosity: f64,
    /// Fluid density used for kinetic-energy scaling.
    pub density: f64,
    /// DNS time step.
    pub time_step: f64,
    /// Final physical time for the benchmark run.
    pub final_time: f64,
    /// Target shell index for the forcing envelope.
    pub forcing_shell: usize,
    /// Half-width of the forcing shell band.
    pub forcing_bandwidth: usize,
    /// Forcing amplitude.
    pub forcing_amplitude: f64,
    /// Number of DNS steps to reuse each forcing realization.
    pub forcing_resample_stride: usize,
    /// Sample the spectrum every `spectrum_stride` steps.
    pub spectrum_stride: usize,
    /// Seed used by the deterministic forcing generator.
    pub seed: u64,
}

impl ForcedTurbulenceBenchmarkConfig {
    /// Create a validated benchmark configuration.
    pub fn new(
        dimensions: (usize, usize, usize),
        length_scale: f64,
        kinematic_viscosity: f64,
        density: f64,
        time_step: f64,
        final_time: f64,
        forcing_shell: usize,
        forcing_bandwidth: usize,
        forcing_amplitude: f64,
        forcing_resample_stride: usize,
        spectrum_stride: usize,
        seed: u64,
    ) -> Result<Self> {
        let (nx, ny, nz) = dimensions;
        if nx == 0 || ny == 0 || nz == 0 {
            return Err(Error::InvalidConfiguration(
                "ForcedTurbulenceBenchmarkConfig: dimensions must be greater than zero".into(),
            ));
        }
        if !length_scale.is_finite() || length_scale <= 0.0 {
            return Err(Error::InvalidConfiguration(
                "ForcedTurbulenceBenchmarkConfig: length_scale must be finite and positive".into(),
            ));
        }
        if !kinematic_viscosity.is_finite() || kinematic_viscosity < 0.0 {
            return Err(Error::InvalidConfiguration(
                "ForcedTurbulenceBenchmarkConfig: kinematic_viscosity must be finite and non-negative".into(),
            ));
        }
        if !density.is_finite() || density <= 0.0 {
            return Err(Error::InvalidConfiguration(
                "ForcedTurbulenceBenchmarkConfig: density must be finite and positive".into(),
            ));
        }
        if !time_step.is_finite() || time_step <= 0.0 {
            return Err(Error::InvalidConfiguration(
                "ForcedTurbulenceBenchmarkConfig: time_step must be finite and positive".into(),
            ));
        }
        if !final_time.is_finite() || final_time <= 0.0 {
            return Err(Error::InvalidConfiguration(
                "ForcedTurbulenceBenchmarkConfig: final_time must be finite and positive".into(),
            ));
        }
        if forcing_shell == 0 {
            return Err(Error::InvalidConfiguration(
                "ForcedTurbulenceBenchmarkConfig: forcing_shell must be greater than zero".into(),
            ));
        }
        if !forcing_amplitude.is_finite() || forcing_amplitude < 0.0 {
            return Err(Error::InvalidConfiguration(
                "ForcedTurbulenceBenchmarkConfig: forcing_amplitude must be finite and non-negative".into(),
            ));
        }
        if forcing_resample_stride == 0 {
            return Err(Error::InvalidConfiguration(
                "ForcedTurbulenceBenchmarkConfig: forcing_resample_stride must be greater than zero".into(),
            ));
        }
        if spectrum_stride == 0 {
            return Err(Error::InvalidConfiguration(
                "ForcedTurbulenceBenchmarkConfig: spectrum_stride must be greater than zero".into(),
            ));
        }

        Ok(Self {
            dimensions,
            length_scale,
            kinematic_viscosity,
            density,
            time_step,
            final_time,
            forcing_shell,
            forcing_bandwidth,
            forcing_amplitude,
            forcing_resample_stride,
            spectrum_stride,
            seed,
        })
    }

    /// Physical domain lengths used by the periodic DNS solver.
    #[must_use]
    pub fn domain_lengths(&self) -> (f64, f64, f64) {
        let domain_length = 2.0 * std::f64::consts::PI * self.length_scale;
        (domain_length, domain_length, domain_length)
    }

    /// Number of DNS steps implied by the benchmark horizon.
    #[must_use]
    pub fn step_count(&self) -> usize {
        let steps = (self.final_time / self.time_step).round();
        steps.max(1.0) as usize
    }

    fn with_runtime_overrides(&self, runtime: &BenchmarkConfig<f64>) -> Self {
        let mut resolved = self.clone();

        if runtime.resolution > 0 {
            resolved.dimensions = (runtime.resolution, runtime.resolution, runtime.resolution);
        }

        if let Some(time_step) = runtime.time_step {
            if time_step.is_finite() && time_step > 0.0 {
                resolved.time_step = time_step;
            }
        }

        let max_steps = runtime.max_iterations.max(1);
        let configured_steps = resolved.step_count();
        let effective_steps = configured_steps.min(max_steps);
        resolved.final_time = effective_steps as f64 * resolved.time_step;

        resolved
    }
}

impl Default for ForcedTurbulenceBenchmarkConfig {
    fn default() -> Self {
        Self {
            dimensions: (16, 16, 16),
            length_scale: 1.0,
            kinematic_viscosity: 0.01,
            density: 1.0,
            time_step: 0.005,
            final_time: 0.05,
            forcing_shell: 2,
            forcing_bandwidth: 1,
            forcing_amplitude: 0.35,
            forcing_resample_stride: 4,
            spectrum_stride: 1,
            seed: 11,
        }
    }
}

/// One stationary forcing checkpoint.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct ForcedTurbulenceCheckpoint {
    /// Step index for the checkpoint.
    pub step: usize,
    /// Physical time of the checkpoint.
    pub time: f64,
    /// Numerically measured kinetic energy.
    pub kinetic_energy: f64,
    /// Absolute energy change relative to the previous checkpoint.
    pub energy_delta: f64,
    /// Dominant shell in the isotropic spectrum.
    pub dominant_shell: usize,
    /// Total enstrophy measured at the checkpoint.
    pub enstrophy: f64,
    /// Forcing epoch that produced the current forcing realization.
    pub forcing_epoch: u64,
    /// Full kinetic-energy spectrum snapshot.
    pub spectrum: KineticEnergySpectrum,
    /// Full enstrophy spectrum snapshot.
    pub enstrophy_spectrum: EnstrophySpectrum,
}

/// Collection of forced-turbulence checkpoint data.
#[derive(Debug, Clone, Default, Serialize, Deserialize)]
pub struct ForcedTurbulenceBenchmarkHistory {
    /// Ordered checkpoints from initial to final state.
    pub checkpoints: Vec<ForcedTurbulenceCheckpoint>,
}

impl ForcedTurbulenceBenchmarkHistory {
    fn push(&mut self, checkpoint: ForcedTurbulenceCheckpoint) {
        self.checkpoints.push(checkpoint);
    }

    /// Return the last checkpoint, if any.
    #[must_use]
    pub fn final_checkpoint(&self) -> Option<&ForcedTurbulenceCheckpoint> {
        self.checkpoints.last()
    }
}

/// Combined benchmark report containing the generic result plus rich history.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct ForcedTurbulenceBenchmarkReport {
    /// Generic validation result expected by the benchmark framework.
    pub result: BenchmarkResult<f64>,
    /// Forced turbulence checkpoint history.
    pub history: ForcedTurbulenceBenchmarkHistory,
}

/// Stationary forced isotropic turbulence benchmark driver.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct ForcedTurbulenceBenchmark3D {
    /// Static benchmark configuration.
    pub config: ForcedTurbulenceBenchmarkConfig,
}

impl ForcedTurbulenceBenchmark3D {
    /// Create a new stationary forced turbulence benchmark.
    pub fn new(config: ForcedTurbulenceBenchmarkConfig) -> Self {
        Self { config }
    }

    /// Run the benchmark and return the rich report, including checkpoint history.
    pub fn run_with_report(
        &self,
        runtime: &BenchmarkConfig<f64>,
    ) -> Result<ForcedTurbulenceBenchmarkReport> {
        let config = self.config.with_runtime_overrides(runtime);
        let dns = PeriodicPseudospectralDns3D::new(PeriodicPseudospectralDnsConfig::new(
            config.dimensions,
            config.domain_lengths(),
            config.kinematic_viscosity,
            config.time_step,
        )?)?;

        let forcing_config = BandLimitedRandomPhaseForcingConfig::new(
            config.dimensions,
            config.domain_lengths(),
            config.forcing_shell,
            config.forcing_bandwidth,
            config.forcing_amplitude,
            config.seed,
        )?;
        let forcing_schedule = TimeResampledBandLimitedForcing3D::new(
            TimeResampledBandLimitedForcingConfig::new(
                forcing_config,
                config.forcing_resample_stride,
            )?,
        )?;

        let mut velocity = Self::zero_velocity_field(config.dimensions);
        let mut history = ForcedTurbulenceBenchmarkHistory::default();
        let start = Instant::now();
        let step_count = config.step_count();
        let probe_index = Self::probe_index(config.dimensions);
        let mut probe_samples = Vec::with_capacity(step_count + 1);
        let mut current_time = 0.0_f64;

        for step in 0..=step_count {
            if step == 0 || step % config.spectrum_stride == 0 || step == step_count {
                let checkpoint = Self::checkpoint(step, current_time, &velocity, &config)?;
                history.push(checkpoint);
                probe_samples.push(Self::probe_sample(&velocity, probe_index));
            }

            if step < step_count {
                let forcing_field = forcing_schedule.sample_physical_forcing_for_step(step)?;
                velocity = dns.advance_euler_with_forcing(&velocity, &forcing_field)?;
                current_time += config.time_step;
            }
        }

        let mut previous_energy: Option<f64> = None;
        for checkpoint in &mut history.checkpoints {
            checkpoint.energy_delta = previous_energy
                .map(|energy| (checkpoint.kinetic_energy - energy).abs())
                .unwrap_or(0.0);
            previous_energy = Some(checkpoint.kinetic_energy);
        }

        let elapsed = start.elapsed().as_secs_f64();
        let mut result = BenchmarkResult::new(self.name());
        result.execution_time = elapsed;
        result.values = history
            .checkpoints
            .iter()
            .map(|checkpoint| checkpoint.kinetic_energy)
            .collect();
        result.errors = history
            .checkpoints
            .iter()
            .map(|checkpoint| checkpoint.energy_delta)
            .collect();
        result.convergence = Self::cumulative_mean_energy(&history);

        if let Some(initial) = history.checkpoints.first() {
            result
                .metrics
                .insert("Initial Energy".to_string(), initial.kinetic_energy);
        }

        if let Some(final_checkpoint) = history.final_checkpoint() {
            result.metrics.insert(
                "Final Energy".to_string(),
                final_checkpoint.kinetic_energy,
            );
            result.metrics.insert(
                "Final Enstrophy".to_string(),
                final_checkpoint.enstrophy,
            );
            result.metrics.insert(
                "Final Dominant Shell".to_string(),
                final_checkpoint.dominant_shell as f64,
            );
            result.metrics.insert(
                "Final Forcing Epoch".to_string(),
                final_checkpoint.forcing_epoch as f64,
            );
        }

        let probe_spectrum = probe_signal_spectrum(&probe_samples, config.time_step)?;
        let probe_autocorrelation = temporal_autocorrelation(
            &probe_samples,
            config.time_step,
            probe_samples.len().saturating_sub(1),
        )?;
        let probe_dominant_bin = probe_spectrum
            .spectral_energy
            .iter()
            .enumerate()
            .max_by(|lhs, rhs| lhs.1.partial_cmp(rhs.1).unwrap_or(std::cmp::Ordering::Equal))
            .map(|(index, _)| index)
            .unwrap_or(0);

        result.metrics.insert(
            "Probe Samples".to_string(),
            probe_samples.len() as f64,
        );
        result.metrics.insert(
            "Probe Mean".to_string(),
            probe_spectrum.mean,
        );
        result.metrics.insert(
            "Probe Variance".to_string(),
            probe_spectrum.variance,
        );
        result.metrics.insert(
            "Probe Dominant Frequency".to_string(),
            probe_spectrum.frequencies_hz[probe_dominant_bin].abs(),
        );
        result.metrics.insert(
            "Probe Lag1 Autocorrelation".to_string(),
            probe_autocorrelation.values.get(1).copied().unwrap_or(1.0),
        );
        result.metrics.insert(
            "Probe Lag2 Autocorrelation".to_string(),
            probe_autocorrelation.values.get(2).copied().unwrap_or(1.0),
        );

        if !history.checkpoints.is_empty() {
            let peak_energy = history
                .checkpoints
                .iter()
                .map(|checkpoint| checkpoint.kinetic_energy)
                .fold(0.0, f64::max);
            let peak_enstrophy = history
                .checkpoints
                .iter()
                .map(|checkpoint| checkpoint.enstrophy)
                .fold(0.0, f64::max);
            let mean_energy = history
                .checkpoints
                .iter()
                .map(|checkpoint| checkpoint.kinetic_energy)
                .sum::<f64>()
                / history.checkpoints.len() as f64;
            let mean_enstrophy = history
                .checkpoints
                .iter()
                .map(|checkpoint| checkpoint.enstrophy)
                .sum::<f64>()
                / history.checkpoints.len() as f64;
            let unique_epochs = history
                .checkpoints
                .iter()
                .map(|checkpoint| checkpoint.forcing_epoch)
                .collect::<std::collections::BTreeSet<_>>()
                .len();

            result
                .metrics
                .insert("Peak Energy".to_string(), peak_energy);
            result
                .metrics
                .insert("Peak Enstrophy".to_string(), peak_enstrophy);
            result
                .metrics
                .insert("Mean Energy".to_string(), mean_energy);
            result
                .metrics
                .insert("Mean Enstrophy".to_string(), mean_enstrophy);
            result.metrics.insert(
                "Spectrum Samples".to_string(),
                history.checkpoints.len() as f64,
            );
            result.metrics.insert(
                "Forced Epochs".to_string(),
                unique_epochs as f64,
            );
        }

        result.metadata.insert(
            "dimensions".to_string(),
            format!(
                "({}, {}, {})",
                config.dimensions.0, config.dimensions.1, config.dimensions.2
            ),
        );
        result.metadata.insert(
            "domain_lengths".to_string(),
            format!(
                "({:.6}, {:.6}, {:.6})",
                config.domain_lengths().0,
                config.domain_lengths().1,
                config.domain_lengths().2
            ),
        );
        result.metadata.insert(
            "time_step".to_string(),
            format!("{:.6}", config.time_step),
        );
        result.metadata.insert(
            "final_time".to_string(),
            format!("{:.6}", config.final_time),
        );
        result.metadata.insert(
            "forcing_shell".to_string(),
            format!("{}", config.forcing_shell),
        );
        result.metadata.insert(
            "forcing_bandwidth".to_string(),
            format!("{}", config.forcing_bandwidth),
        );
        result.metadata.insert(
            "forcing_amplitude".to_string(),
            format!("{:.6}", config.forcing_amplitude),
        );
        result.metadata.insert(
            "forcing_resample_stride".to_string(),
            format!("{}", config.forcing_resample_stride),
        );
        result.metadata.insert(
            "spectrum_stride".to_string(),
            format!("{}", config.spectrum_stride),
        );
        result.metadata.insert(
            "seed".to_string(),
            format!("{}", config.seed),
        );

        Ok(ForcedTurbulenceBenchmarkReport { result, history })
    }

    fn zero_velocity_field(dimensions: (usize, usize, usize)) -> VelocityField<f64> {
        let total_points = dimensions.0 * dimensions.1 * dimensions.2;
        VelocityField {
            components: vec![Vector3::zeros(); total_points],
            dimensions,
        }
    }

    fn checkpoint(
        step: usize,
        time: f64,
        velocity: &VelocityField<f64>,
        config: &ForcedTurbulenceBenchmarkConfig,
    ) -> Result<ForcedTurbulenceCheckpoint> {
        let spectrum = kinetic_energy_spectrum(velocity)?;
        let enstrophy_spectrum = enstrophy_spectrum(velocity, config.domain_lengths())?;
        let kinetic_energy = spectrum.total_energy * config.density;
        let dominant_shell = spectrum
            .shell_energy
            .iter()
            .enumerate()
            .max_by(|lhs, rhs| lhs.1.partial_cmp(rhs.1).unwrap_or(std::cmp::Ordering::Equal))
            .map(|(index, _)| index)
            .unwrap_or(0);
        let enstrophy = enstrophy_spectrum.total_enstrophy;
        let forcing_epoch = (step / config.forcing_resample_stride) as u64;

        Ok(ForcedTurbulenceCheckpoint {
            step,
            time,
            kinetic_energy,
            energy_delta: 0.0,
            dominant_shell,
            enstrophy,
            forcing_epoch,
            spectrum,
            enstrophy_spectrum,
        })
    }

    fn cumulative_mean_energy(history: &ForcedTurbulenceBenchmarkHistory) -> Vec<f64> {
        let mut means = Vec::with_capacity(history.checkpoints.len());
        let mut running_sum = 0.0;

        for (index, checkpoint) in history.checkpoints.iter().enumerate() {
            running_sum += checkpoint.kinetic_energy;
            means.push(running_sum / (index + 1) as f64);
        }

        means
    }

    fn probe_index(dimensions: (usize, usize, usize)) -> usize {
        let (nx, ny, nz) = dimensions;
        (nz / 2) * nx * ny + (ny / 2) * nx + (nx / 2)
    }

    fn probe_sample(velocity: &VelocityField<f64>, probe_index: usize) -> f64 {
        velocity.components[probe_index].norm()
    }

    /// Validate the forced-turbulence history.
    #[must_use]
    pub fn validate_history(&self, history: &ForcedTurbulenceBenchmarkHistory) -> bool {
        let Some(final_checkpoint) = history.final_checkpoint() else {
            return false;
        };

        let saw_nonzero_energy = history
            .checkpoints
            .iter()
            .any(|checkpoint| checkpoint.kinetic_energy > 0.0);
        let saw_resampled_forcing = history
            .checkpoints
            .windows(2)
            .any(|window| window[0].forcing_epoch != window[1].forcing_epoch);

        saw_nonzero_energy
            && saw_resampled_forcing
            && final_checkpoint.kinetic_energy > 0.0
            && final_checkpoint.enstrophy.is_finite()
            && final_checkpoint.enstrophy >= 0.0
            && history
                .checkpoints
                .iter()
                .all(|checkpoint| checkpoint.kinetic_energy.is_finite())
            && history
                .checkpoints
                .iter()
                .all(|checkpoint| checkpoint.enstrophy.is_finite())
            && history
                .checkpoints
                .iter()
                .all(|checkpoint| checkpoint.spectrum.total_energy.is_finite())
            && history
                .checkpoints
                .iter()
                .all(|checkpoint| checkpoint.spectrum.total_energy >= 0.0)
            && history
                .checkpoints
                .iter()
                .all(|checkpoint| checkpoint.enstrophy_spectrum.total_enstrophy.is_finite())
            && history
                .checkpoints
                .iter()
                .all(|checkpoint| checkpoint.enstrophy_spectrum.total_enstrophy >= 0.0)
            && history
                .checkpoints
                .iter()
                .all(|checkpoint| checkpoint.dominant_shell < checkpoint.spectrum.shell_energy.len())
    }
}

impl Default for ForcedTurbulenceBenchmark3D {
    fn default() -> Self {
        Self::new(ForcedTurbulenceBenchmarkConfig::default())
    }
}

impl Benchmark<f64> for ForcedTurbulenceBenchmark3D {
    fn name(&self) -> &'static str {
        "3D Stationary Forced Turbulence"
    }

    fn description(&self) -> &'static str {
        "Apollo-backed periodic pseudospectral turbulence driven by time-resampled band-limited forcing."
    }

    fn run(&self, config: &BenchmarkConfig<f64>) -> Result<BenchmarkResult<f64>> {
        self.run_with_report(config).map(|report| report.result)
    }

    fn reference_solution(&self) -> Option<BenchmarkResult<f64>> {
        None
    }

    fn validate(&self, result: &BenchmarkResult<f64>) -> Result<bool> {
        let Some(final_energy) = result.metrics.get("Final Energy") else {
            return Ok(false);
        };
        let Some(final_enstrophy) = result.metrics.get("Final Enstrophy") else {
            return Ok(false);
        };
        let Some(peak_energy) = result.metrics.get("Peak Energy") else {
            return Ok(false);
        };
        let Some(peak_enstrophy) = result.metrics.get("Peak Enstrophy") else {
            return Ok(false);
        };
        let Some(forced_epochs) = result.metrics.get("Forced Epochs") else {
            return Ok(false);
        };
        let Some(probe_variance) = result.metrics.get("Probe Variance") else {
            return Ok(false);
        };

        let finite_history = result.values.iter().all(|value| value.is_finite())
            && result.errors.iter().all(|value| value.is_finite())
            && result.convergence.iter().all(|value| value.is_finite());

        Ok(finite_history
            && final_energy.is_finite()
            && *final_energy > 0.0
            && final_enstrophy.is_finite()
            && *final_enstrophy >= 0.0
            && peak_energy.is_finite()
            && *peak_energy >= *final_energy
            && peak_enstrophy.is_finite()
            && *peak_enstrophy >= *final_enstrophy
            && *forced_epochs >= 2.0
            && probe_variance.is_finite()
            && *probe_variance >= 0.0
            && result.values.len() >= 2)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::benchmarks::BenchmarkRunner;

    #[test]
    fn benchmark_reports_nonzero_energy_and_spectra() {
        let benchmark = ForcedTurbulenceBenchmark3D::new(ForcedTurbulenceBenchmarkConfig::new(
            (8, 8, 8),
            1.0,
            0.01,
            1.0,
            0.005,
            0.04,
            2,
            1,
            0.35,
            2,
            1,
            17,
        )
        .expect("benchmark config should be valid"));

        let runtime = BenchmarkConfig {
            resolution: 8,
            tolerance: 0.25,
            max_iterations: 16,
            reynolds_number: 120.0,
            time_step: Some(0.005),
            parallel: false,
        };

        let report = benchmark
            .run_with_report(&runtime)
            .expect("forced turbulence benchmark should run");

        assert!(!report.history.checkpoints.is_empty());
        assert!(report
            .history
            .checkpoints
            .first()
            .expect("history should contain an initial checkpoint")
            .kinetic_energy
            >= 0.0);
        assert!(report
            .history
            .checkpoints
            .last()
            .expect("history should contain a final checkpoint")
            .kinetic_energy
            > 0.0);
        assert!(report
            .history
            .checkpoints
            .last()
            .expect("history should contain a final checkpoint")
            .enstrophy
            >= 0.0);
        assert!(benchmark.validate_history(&report.history));
        assert!(report.result.metrics.contains_key("Spectrum Samples"));
        assert!(report.result.metrics.contains_key("Forced Epochs"));
        assert!(report.result.metrics.contains_key("Final Enstrophy"));
        assert!(report.result.metrics.contains_key("Peak Enstrophy"));
        assert!(report.result.metrics.contains_key("Probe Variance"));
        assert!(report.result.values.len() >= 2);
    }

    #[test]
    fn benchmark_trait_adapter_runs_through_validation_framework() {
        let benchmark = ForcedTurbulenceBenchmark3D::default();
        let runtime = BenchmarkConfig {
            resolution: 8,
            tolerance: 0.25,
            max_iterations: 16,
            reynolds_number: 120.0,
            time_step: Some(0.005),
            parallel: false,
        };

        let result = BenchmarkRunner::run_benchmark(&benchmark, &runtime)
            .expect("benchmark runner should execute the forced turbulence benchmark");

        assert!(benchmark.validate(&result).expect("validation should be computable"));
        assert!(result.metrics.contains_key("Final Energy"));
        assert!(result.metrics.contains_key("Peak Energy"));
        assert!(result.metrics.contains_key("Final Enstrophy"));
        assert!(result.metrics.contains_key("Peak Enstrophy"));
        assert!(result.metrics.contains_key("Probe Dominant Frequency"));
    }
}