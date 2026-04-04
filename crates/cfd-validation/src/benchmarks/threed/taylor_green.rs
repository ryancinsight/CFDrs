//! 3D Taylor-Green vortex benchmark with Apollo-backed periodic DNS.
//!
//! This benchmark bridges the analytical Taylor-Green solution in
//! `cfd-validation` with the Apollo-backed pseudospectral DNS stepper in
//! `cfd-3d`. It records an energy history and a spectrum history so the
//! validation layer can exercise both the time integrator and the FFT-backed
//! diagnostics on the same run.

use super::super::{Benchmark, BenchmarkConfig, BenchmarkResult};
use crate::analytical::{AnalyticalSolution, TaylorGreenVortex};
use cfd_3d::spectral::{
    kinetic_energy_spectrum, KineticEnergySpectrum, PeriodicPseudospectralDns3D,
    PeriodicPseudospectralDnsConfig,
};
use cfd_core::error::{Error, Result};
use cfd_core::physics::fluid_dynamics::VelocityField;
use serde::{Deserialize, Serialize};
use std::time::Instant;

/// Configuration for the Taylor-Green DNS benchmark.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct TaylorGreenBenchmarkConfig {
    /// Grid dimensions `(nx, ny, nz)`.
    pub dimensions: (usize, usize, usize),
    /// Characteristic length scale `L` used by the analytical solution.
    pub length_scale: f64,
    /// Characteristic velocity scale `U`.
    pub velocity_scale: f64,
    /// Kinematic viscosity `nu`.
    pub kinematic_viscosity: f64,
    /// Fluid density used to scale the energy budget.
    pub density: f64,
    /// Explicit DNS time step.
    pub time_step: f64,
    /// Final physical time for the benchmark run.
    pub final_time: f64,
    /// Sample the spectrum every `spectrum_stride` steps.
    pub spectrum_stride: usize,
    /// Relative energy error threshold used by `validate()`.
    pub validation_tolerance: f64,
}

impl TaylorGreenBenchmarkConfig {
    /// Create a validated benchmark configuration.
    pub fn new(
        dimensions: (usize, usize, usize),
        length_scale: f64,
        velocity_scale: f64,
        kinematic_viscosity: f64,
        density: f64,
        time_step: f64,
        final_time: f64,
        spectrum_stride: usize,
        validation_tolerance: f64,
    ) -> Result<Self> {
        let (nx, ny, nz) = dimensions;
        if nx == 0 || ny == 0 || nz == 0 {
            return Err(Error::InvalidConfiguration(
                "TaylorGreenBenchmarkConfig: dimensions must be greater than zero".into(),
            ));
        }
        if !length_scale.is_finite() || length_scale <= 0.0 {
            return Err(Error::InvalidConfiguration(
                "TaylorGreenBenchmarkConfig: length_scale must be finite and positive".into(),
            ));
        }
        if !velocity_scale.is_finite() || velocity_scale <= 0.0 {
            return Err(Error::InvalidConfiguration(
                "TaylorGreenBenchmarkConfig: velocity_scale must be finite and positive".into(),
            ));
        }
        if !kinematic_viscosity.is_finite() || kinematic_viscosity < 0.0 {
            return Err(Error::InvalidConfiguration(
                "TaylorGreenBenchmarkConfig: kinematic_viscosity must be finite and non-negative".into(),
            ));
        }
        if !density.is_finite() || density <= 0.0 {
            return Err(Error::InvalidConfiguration(
                "TaylorGreenBenchmarkConfig: density must be finite and positive".into(),
            ));
        }
        if !time_step.is_finite() || time_step <= 0.0 {
            return Err(Error::InvalidConfiguration(
                "TaylorGreenBenchmarkConfig: time_step must be finite and positive".into(),
            ));
        }
        if !final_time.is_finite() || final_time <= 0.0 {
            return Err(Error::InvalidConfiguration(
                "TaylorGreenBenchmarkConfig: final_time must be finite and positive".into(),
            ));
        }
        if spectrum_stride == 0 {
            return Err(Error::InvalidConfiguration(
                "TaylorGreenBenchmarkConfig: spectrum_stride must be greater than zero".into(),
            ));
        }
        if !validation_tolerance.is_finite() || validation_tolerance < 0.0 {
            return Err(Error::InvalidConfiguration(
                "TaylorGreenBenchmarkConfig: validation_tolerance must be finite and non-negative".into(),
            ));
        }

        Ok(Self {
            dimensions,
            length_scale,
            velocity_scale,
            kinematic_viscosity,
            density,
            time_step,
            final_time,
            spectrum_stride,
            validation_tolerance,
        })
    }

    /// Physical domain lengths used by the periodic DNS solver.
    #[must_use]
    pub fn domain_lengths(&self) -> (f64, f64, f64) {
        let domain_length = 2.0 * self.length_scale;
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

        if runtime.reynolds_number.is_finite() && runtime.reynolds_number > 0.0 {
            resolved.kinematic_viscosity =
                resolved.velocity_scale * resolved.length_scale / runtime.reynolds_number;
        }

        let max_steps = runtime.max_iterations.max(1);
        let configured_steps = resolved.step_count();
        let effective_steps = configured_steps.min(max_steps);
        resolved.final_time = effective_steps as f64 * resolved.time_step;

        resolved
    }
}

impl Default for TaylorGreenBenchmarkConfig {
    fn default() -> Self {
        Self {
            dimensions: (16, 16, 16),
            length_scale: 1.0,
            velocity_scale: 1.0,
            kinematic_viscosity: 0.01,
            density: 1.0,
            time_step: 0.0025,
            final_time: 0.025,
            spectrum_stride: 1,
            validation_tolerance: 0.25,
        }
    }
}

/// One Taylor-Green checkpoint with both scalar and spectral summaries.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct TaylorGreenCheckpoint {
    /// Step index for the checkpoint.
    pub step: usize,
    /// Physical time of the checkpoint.
    pub time: f64,
    /// Numerically measured kinetic energy.
    pub numerical_energy: f64,
    /// Analytic kinetic energy at the same time.
    pub analytic_energy: f64,
    /// Relative energy error against the analytic decay curve.
    pub relative_energy_error: f64,
    /// Dominant shell in the isotropic spectrum.
    pub dominant_shell: usize,
    /// Full kinetic-energy spectrum snapshot.
    pub spectrum: KineticEnergySpectrum,
}

/// Collection of Taylor-Green checkpoint data.
#[derive(Debug, Clone, Default, Serialize, Deserialize)]
pub struct TaylorGreenBenchmarkHistory {
    /// Ordered checkpoints from initial to final state.
    pub checkpoints: Vec<TaylorGreenCheckpoint>,
}

impl TaylorGreenBenchmarkHistory {
    fn push(&mut self, checkpoint: TaylorGreenCheckpoint) {
        self.checkpoints.push(checkpoint);
    }

    /// Return the last checkpoint, if any.
    #[must_use]
    pub fn final_checkpoint(&self) -> Option<&TaylorGreenCheckpoint> {
        self.checkpoints.last()
    }
}

/// Combined benchmark report containing the generic result plus rich history.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct TaylorGreenBenchmarkReport {
    /// Generic validation result expected by the benchmark framework.
    pub result: BenchmarkResult<f64>,
    /// Taylor-Green checkpoint history.
    pub history: TaylorGreenBenchmarkHistory,
}

/// Taylor-Green vortex benchmark driver.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct TaylorGreenBenchmark3D {
    /// Static benchmark configuration.
    pub config: TaylorGreenBenchmarkConfig,
}

impl TaylorGreenBenchmark3D {
    /// Create a new Taylor-Green benchmark.
    pub fn new(config: TaylorGreenBenchmarkConfig) -> Self {
        Self { config }
    }

    /// Run the benchmark and return the rich report, including checkpoint history.
    pub fn run_with_report(
        &self,
        runtime: &BenchmarkConfig<f64>,
    ) -> Result<TaylorGreenBenchmarkReport> {
        let config = self.config.with_runtime_overrides(runtime);
        let solution = TaylorGreenVortex::create(
            config.length_scale,
            config.velocity_scale,
            config.kinematic_viscosity,
            config.density,
            true,
        );

        let dns = PeriodicPseudospectralDns3D::new(PeriodicPseudospectralDnsConfig::new(
            config.dimensions,
            config.domain_lengths(),
            config.kinematic_viscosity,
            config.time_step,
        )?)?;

        let mut velocity = Self::initial_velocity_field(&config, &solution);
        let mut history = TaylorGreenBenchmarkHistory::default();
        let start = Instant::now();
        let step_count = config.step_count();
        let mut current_time = 0.0_f64;

        for step in 0..=step_count {
            if step == 0 || step % config.spectrum_stride == 0 || step == step_count {
                let checkpoint = Self::checkpoint(step, current_time, &velocity, &solution)?;
                history.push(checkpoint);
            }

            if step < step_count {
                velocity = dns.advance_euler(&velocity)?;
                current_time += config.time_step;
            }
        }

        let energy_scale = history
            .checkpoints
            .first()
            .map(|checkpoint| {
                if checkpoint.numerical_energy.abs() > f64::EPSILON {
                    checkpoint.analytic_energy / checkpoint.numerical_energy
                } else {
                    1.0
                }
            })
            .unwrap_or(1.0);

        for checkpoint in &mut history.checkpoints {
            checkpoint.relative_energy_error = if checkpoint.analytic_energy.abs() > f64::EPSILON {
                ((checkpoint.numerical_energy * energy_scale) - checkpoint.analytic_energy).abs()
                    / checkpoint.analytic_energy.abs()
            } else {
                0.0
            };
        }

        let elapsed = start.elapsed().as_secs_f64();
        let mut result = BenchmarkResult::new(self.name());
        result.execution_time = elapsed;
        result.values = history
            .checkpoints
            .iter()
            .map(|checkpoint| checkpoint.numerical_energy * energy_scale)
            .collect();
        result.errors = history
            .checkpoints
            .iter()
            .map(|checkpoint| checkpoint.relative_energy_error)
            .collect();
        result.convergence = Self::energy_deltas(&history, energy_scale);

        if let Some(initial) = history.checkpoints.first() {
            result
                .metrics
                .insert("Initial Numerical Energy".to_string(), initial.numerical_energy * energy_scale);
            result
                .metrics
                .insert("Initial Analytic Energy".to_string(), initial.analytic_energy);
        }

        if let Some(final_checkpoint) = history.final_checkpoint() {
            result.metrics.insert(
                "Final Numerical Energy".to_string(),
                final_checkpoint.numerical_energy * energy_scale,
            );
            result.metrics.insert(
                "Final Analytic Energy".to_string(),
                final_checkpoint.analytic_energy,
            );
            result.metrics.insert(
                "Final Relative Energy Error".to_string(),
                final_checkpoint.relative_energy_error,
            );
            result.metrics.insert(
                "Final Dominant Shell".to_string(),
                final_checkpoint.dominant_shell as f64,
            );
        }

        result
            .metrics
            .insert("Energy Scale Factor".to_string(), energy_scale);

        if !history.checkpoints.is_empty() {
            let max_relative_error = history
                .checkpoints
                .iter()
                .map(|checkpoint| checkpoint.relative_energy_error)
                .fold(0.0, f64::max);
            let mean_relative_error = history
                .checkpoints
                .iter()
                .map(|checkpoint| checkpoint.relative_energy_error)
                .sum::<f64>()
                / history.checkpoints.len() as f64;
            result.metrics.insert(
                "Max Relative Energy Error".to_string(),
                max_relative_error,
            );
            result.metrics.insert(
                "Mean Relative Energy Error".to_string(),
                mean_relative_error,
            );
            result.metrics.insert(
                "Spectrum Samples".to_string(),
                history.checkpoints.len() as f64,
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
            "step_count".to_string(),
            format!("{}", step_count),
        );
        result.metadata.insert(
            "spectrum_stride".to_string(),
            format!("{}", config.spectrum_stride),
        );
        result.metadata.insert(
            "validation_tolerance".to_string(),
            format!("{:.6}", config.validation_tolerance),
        );
        result.metadata.insert(
            "energy_scale_factor".to_string(),
            format!("{:.6}", energy_scale),
        );

        Ok(TaylorGreenBenchmarkReport { result, history })
    }

    fn initial_velocity_field(
        config: &TaylorGreenBenchmarkConfig,
        solution: &TaylorGreenVortex<f64>,
    ) -> VelocityField<f64> {
        Self::velocity_field_at_time(config, solution, 0.0)
    }

    fn velocity_field_at_time(
        config: &TaylorGreenBenchmarkConfig,
        solution: &TaylorGreenVortex<f64>,
        time: f64,
    ) -> VelocityField<f64> {
        let (lx, ly, lz) = config.domain_lengths();
        let (nx, ny, nz) = config.dimensions;
        let mut components = Vec::with_capacity(nx * ny * nz);

        for k in 0..nz {
            let z = lz * k as f64 / nz as f64;
            for j in 0..ny {
                let y = ly * j as f64 / ny as f64;
                for i in 0..nx {
                    let x = lx * i as f64 / nx as f64;
                    components.push(solution.evaluate(x, y, z, time));
                }
            }
        }

        VelocityField {
            components,
            dimensions: config.dimensions,
        }
    }

    fn checkpoint(
        step: usize,
        time: f64,
        velocity: &VelocityField<f64>,
        solution: &TaylorGreenVortex<f64>,
    ) -> Result<TaylorGreenCheckpoint> {
        let spectrum = kinetic_energy_spectrum(velocity)?;
        let numerical_energy = spectrum.total_energy * solution.density;
        let analytic_energy = solution.kinetic_energy(time);
        let relative_energy_error = if analytic_energy.abs() > f64::EPSILON {
            (numerical_energy - analytic_energy).abs() / analytic_energy.abs()
        } else {
            0.0
        };
        let dominant_shell = spectrum
            .shell_energy
            .iter()
            .enumerate()
            .max_by(|lhs, rhs| lhs.1.partial_cmp(rhs.1).unwrap_or(std::cmp::Ordering::Equal))
            .map(|(index, _)| index)
            .unwrap_or(0);

        Ok(TaylorGreenCheckpoint {
            step,
            time,
            numerical_energy,
            analytic_energy,
            relative_energy_error,
            dominant_shell,
            spectrum,
        })
    }

    fn energy_deltas(history: &TaylorGreenBenchmarkHistory, energy_scale: f64) -> Vec<f64> {
        let mut deltas = Vec::with_capacity(history.checkpoints.len());
        if let Some(first) = history.checkpoints.first() {
            deltas.push(0.0);
            let mut previous = first.numerical_energy * energy_scale;
            for checkpoint in history.checkpoints.iter().skip(1) {
                let current = checkpoint.numerical_energy * energy_scale;
                deltas.push((current - previous).abs());
                previous = current;
            }
        }
        deltas
    }

    /// Validate the benchmark-specific energy history.
    #[must_use]
    pub fn validate_history(&self, history: &TaylorGreenBenchmarkHistory) -> bool {
        let Some(final_checkpoint) = history.final_checkpoint() else {
            return false;
        };

        let max_relative_error = history
            .checkpoints
            .iter()
            .map(|checkpoint| checkpoint.relative_energy_error)
            .fold(0.0, f64::max);

        max_relative_error <= self.config.validation_tolerance
            && final_checkpoint.numerical_energy.is_finite()
            && final_checkpoint.analytic_energy.is_finite()
            && history
                .checkpoints
                .iter()
                .all(|checkpoint| checkpoint.numerical_energy.is_finite())
            && history
                .checkpoints
                .iter()
                .all(|checkpoint| checkpoint.analytic_energy.is_finite())
    }
}

impl Default for TaylorGreenBenchmark3D {
    fn default() -> Self {
        Self::new(TaylorGreenBenchmarkConfig::default())
    }
}

impl Benchmark<f64> for TaylorGreenBenchmark3D {
    fn name(&self) -> &'static str {
        "3D Taylor-Green DNS"
    }

    fn description(&self) -> &'static str {
        "Apollo-backed periodic pseudospectral Taylor-Green benchmark with energy and spectrum histories."
    }

    fn run(&self, config: &BenchmarkConfig<f64>) -> Result<BenchmarkResult<f64>> {
        self.run_with_report(config).map(|report| report.result)
    }

    fn reference_solution(&self) -> Option<BenchmarkResult<f64>> {
        let mut reference = BenchmarkResult::new(self.name());
        let solution = TaylorGreenVortex::create(
            self.config.length_scale,
            self.config.velocity_scale,
            self.config.kinematic_viscosity,
            self.config.density,
            true,
        );
        reference
            .metrics
            .insert("Initial Analytic Energy".to_string(), solution.kinetic_energy(0.0));
        Some(reference)
    }

    fn validate(&self, result: &BenchmarkResult<f64>) -> Result<bool> {
        let Some(final_relative_error) = result.metrics.get("Final Relative Energy Error") else {
            return Ok(false);
        };
        let Some(initial_energy) = result.metrics.get("Initial Numerical Energy") else {
            return Ok(false);
        };
        let Some(final_energy) = result.metrics.get("Final Numerical Energy") else {
            return Ok(false);
        };

        let finite_history = result.values.iter().all(|value| value.is_finite())
            && result.errors.iter().all(|value| value.is_finite())
            && result.convergence.iter().all(|value| value.is_finite());

        let decays_overall = *final_energy <= *initial_energy * (1.0 + self.config.validation_tolerance);

        Ok(finite_history
            && final_relative_error.is_finite()
            && *final_relative_error <= self.config.validation_tolerance
            && decays_overall)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::benchmarks::BenchmarkRunner;

    #[test]
    fn benchmark_report_records_energy_decay_and_spectra() {
        let benchmark = TaylorGreenBenchmark3D::new(TaylorGreenBenchmarkConfig::new(
            (8, 8, 8),
            1.0,
            1.0,
            0.01,
            1.0,
            0.0025,
            0.025,
            1,
            0.25,
        )
        .expect("benchmark config should be valid"));

        let runtime = BenchmarkConfig {
            resolution: 8,
            tolerance: 0.25,
            max_iterations: 32,
            reynolds_number: 100.0,
            time_step: Some(0.0025),
            parallel: false,
        };

        let report = benchmark
            .run_with_report(&runtime)
            .expect("Taylor-Green benchmark should run");

        assert!(!report.history.checkpoints.is_empty());
        assert_eq!(report.result.values.len(), report.history.checkpoints.len());
        assert_eq!(report.result.errors.len(), report.history.checkpoints.len());
        assert!(report
            .history
            .checkpoints
            .first()
            .expect("history should contain an initial checkpoint")
            .numerical_energy
            > report
                .history
                .checkpoints
                .last()
                .expect("history should contain a final checkpoint")
                .numerical_energy);
        assert!(benchmark.validate_history(&report.history));
        assert!(report
            .history
            .checkpoints
            .iter()
            .all(|checkpoint| checkpoint.dominant_shell <= checkpoint.spectrum.shell_energy.len()));
    }

    #[test]
    fn benchmark_trait_adapter_runs_through_validation_framework() {
        let benchmark = TaylorGreenBenchmark3D::default();
        let runtime = BenchmarkConfig {
            resolution: 8,
            tolerance: 0.25,
            max_iterations: 32,
            reynolds_number: 100.0,
            time_step: Some(0.0025),
            parallel: false,
        };

        let result = BenchmarkRunner::run_benchmark(&benchmark, &runtime)
            .expect("benchmark runner should execute the Taylor-Green benchmark");

        assert!(benchmark.validate(&result).expect("validation should be computable"));
        assert!(result.metrics.contains_key("Final Relative Energy Error"));
        assert!(result.metrics.contains_key("Spectrum Samples"));
        assert!(result.values.len() >= 2);
    }

    #[test]
    fn benchmark_initial_condition_matches_analytical_taylor_green() {
        let benchmark = TaylorGreenBenchmark3D::default();
        let solution = TaylorGreenVortex::create(
            benchmark.config.length_scale,
            benchmark.config.velocity_scale,
            benchmark.config.kinematic_viscosity,
            benchmark.config.density,
            true,
        );
        let config = benchmark.config.clone();
        let field = TaylorGreenBenchmark3D::initial_velocity_field(&config, &solution);

        assert_eq!(field.dimensions, config.dimensions);
        assert_eq!(field.components.len(), config.dimensions.0 * config.dimensions.1 * config.dimensions.2);
        assert!(field.components.iter().any(|sample| sample.norm() > 0.0));
    }
}
