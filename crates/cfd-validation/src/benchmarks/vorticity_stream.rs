//! 2D vorticity-streamfunction cavity benchmark.
//!
//! This benchmark exercises the 2D vorticity-stream solver on a lid-driven
//! cavity, reporting convergence, centerline profiles, enstrophy, and basic
//! divergence diagnostics for fast incompressible microfluidic analysis.

use super::cavity::LidDrivenCavity;
use super::{Benchmark, BenchmarkConfig, BenchmarkResult};
use cfd_2d::grid::array2d::Array2D;
use cfd_2d::grid::StructuredGrid2D;
use cfd_2d::physics::vorticity_stream::VorticityStreamConfig;
use cfd_2d::physics::VorticityStreamSolver;
use cfd_core::compute::solver::SolverConfig;
use cfd_core::error::{Error, Result};
use nalgebra::Vector2;
use serde::{Deserialize, Serialize};
use std::time::Instant;

/// Centerline profile pair: `(y, u)` vertical and `(x, v)` horizontal samples.
type CenterlineProfiles = (Vec<(f64, f64)>, Vec<(f64, f64)>);

/// Configuration for the vorticity-stream cavity benchmark.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)
]
pub struct VorticityStreamCavityConfig {
    /// Number of grid points `(nx, ny)`.
    pub grid_points: (usize, usize),
    /// Cavity size.
    pub size: f64,
    /// Lid velocity.
    pub lid_velocity: f64,
    /// Reynolds number.
    pub reynolds: f64,
    /// Time step.
    pub time_step: f64,
    /// Maximum number of explicit steps.
    pub max_steps: usize,
    /// Sample every `sample_stride` steps.
    pub sample_stride: usize,
    /// Stream-function convergence tolerance.
    pub stream_tolerance: f64,
    /// Maximum SOR iterations for the Poisson solve.
    pub poisson_iterations: usize,
}

impl VorticityStreamCavityConfig {
    /// Create a validated benchmark configuration.
    pub fn new(
        grid_points: (usize, usize),
        size: f64,
        lid_velocity: f64,
        reynolds: f64,
        time_step: f64,
        max_steps: usize,
        sample_stride: usize,
        stream_tolerance: f64,
        poisson_iterations: usize,
    ) -> Result<Self> {
        let (nx, ny) = grid_points;
        if nx < 3 || ny < 3 {
            return Err(Error::InvalidConfiguration(
                "VorticityStreamCavityConfig: grid_points must be at least 3x3".into(),
            ));
        }
        if !size.is_finite() || size <= 0.0 {
            return Err(Error::InvalidConfiguration(
                "VorticityStreamCavityConfig: size must be finite and positive".into(),
            ));
        }
        if !lid_velocity.is_finite() || lid_velocity <= 0.0 {
            return Err(Error::InvalidConfiguration(
                "VorticityStreamCavityConfig: lid_velocity must be finite and positive".into(),
            ));
        }
        if !reynolds.is_finite() || reynolds <= 0.0 {
            return Err(Error::InvalidConfiguration(
                "VorticityStreamCavityConfig: reynolds must be finite and positive".into(),
            ));
        }
        if !time_step.is_finite() || time_step <= 0.0 {
            return Err(Error::InvalidConfiguration(
                "VorticityStreamCavityConfig: time_step must be finite and positive".into(),
            ));
        }
        if max_steps == 0 {
            return Err(Error::InvalidConfiguration(
                "VorticityStreamCavityConfig: max_steps must be greater than zero".into(),
            ));
        }
        if sample_stride == 0 {
            return Err(Error::InvalidConfiguration(
                "VorticityStreamCavityConfig: sample_stride must be greater than zero".into(),
            ));
        }
        if !stream_tolerance.is_finite() || stream_tolerance <= 0.0 {
            return Err(Error::InvalidConfiguration(
                "VorticityStreamCavityConfig: stream_tolerance must be finite and positive".into(),
            ));
        }
        if poisson_iterations == 0 {
            return Err(Error::InvalidConfiguration(
                "VorticityStreamCavityConfig: poisson_iterations must be greater than zero".into(),
            ));
        }

        Ok(Self {
            grid_points,
            size,
            lid_velocity,
            reynolds,
            time_step,
            max_steps,
            sample_stride,
            stream_tolerance,
            poisson_iterations,
        })
    }

    /// Physical domain bounds.
    #[must_use]
    pub fn domain_bounds(&self) -> (f64, f64, f64, f64) {
        (0.0, self.size, 0.0, self.size)
    }

    /// Total step count used for the transient solve.
    #[must_use]
    pub fn step_count(&self) -> usize {
        self.max_steps.max(1)
    }

    fn with_runtime_overrides(&self, runtime: &BenchmarkConfig<f64>) -> Self {
        let mut resolved = self.clone();

        if runtime.resolution >= 3 {
            resolved.grid_points = (runtime.resolution, runtime.resolution);
        }

        if let Some(time_step) = runtime.time_step {
            if time_step.is_finite() && time_step > 0.0 {
                resolved.time_step = time_step;
            }
        }

        if runtime.tolerance.is_finite() && runtime.tolerance > 0.0 {
            resolved.stream_tolerance = runtime.tolerance;
        }

        if runtime.reynolds_number.is_finite() && runtime.reynolds_number > 0.0 {
            resolved.reynolds = runtime.reynolds_number;
        }

        let capped_steps = runtime.max_iterations.max(1);
        resolved.max_steps = resolved.max_steps.min(capped_steps);
        resolved.poisson_iterations = resolved.poisson_iterations.min(capped_steps);

        resolved
    }

    fn has_ghia_reference(&self) -> bool {
        (self.reynolds - 100.0).abs() <= 1e-9
    }
}

impl Default for VorticityStreamCavityConfig {
    fn default() -> Self {
        Self {
            grid_points: (33, 33),
            size: 1.0,
            lid_velocity: 1.0,
            reynolds: 100.0,
            time_step: 0.001,
            max_steps: 120,
            sample_stride: 6,
            stream_tolerance: 1e-6,
            poisson_iterations: 120,
        }
    }
}

/// One checkpoint of the vorticity-stream cavity solve.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct VorticityStreamCavityCheckpoint {
    /// Step index.
    pub step: usize,
    /// Physical time.
    pub time: f64,
    /// Maximum change in stream function since the previous step.
    pub residual: f64,
    /// Stream function at the cavity center.
    pub stream_center: f64,
    /// Vorticity at the cavity center.
    pub vorticity_center: f64,
    /// Maximum speed across the grid.
    pub max_speed: f64,
    /// Maximum discrete divergence across the grid.
    pub max_divergence: f64,
    /// Total enstrophy of the current flow field.
    pub enstrophy: f64,
}

/// Checkpoint history for the vorticity-stream cavity benchmark.
#[derive(Debug, Clone, Default, Serialize, Deserialize)]
pub struct VorticityStreamCavityHistory {
    /// Ordered checkpoints from initial to final state.
    pub checkpoints: Vec<VorticityStreamCavityCheckpoint>,
}

impl VorticityStreamCavityHistory {
    fn push(&mut self, checkpoint: VorticityStreamCavityCheckpoint) {
        self.checkpoints.push(checkpoint);
    }

    /// Return the last checkpoint, if any.
    #[must_use]
    pub fn final_checkpoint(&self) -> Option<&VorticityStreamCavityCheckpoint> {
        self.checkpoints.last()
    }
}

/// Rich report for the vorticity-stream cavity benchmark.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct VorticityStreamCavityReport {
    /// Generic validation result expected by the benchmark framework.
    pub result: BenchmarkResult<f64>,
    /// Checkpoint history.
    pub history: VorticityStreamCavityHistory,
}

/// 2D lid-driven cavity benchmark solved with the vorticity-stream formulation.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct VorticityStreamCavityBenchmark {
    /// Static benchmark configuration.
    pub config: VorticityStreamCavityConfig,
}

impl VorticityStreamCavityBenchmark {
    /// Create a new benchmark.
    pub fn new(config: VorticityStreamCavityConfig) -> Self {
        Self { config }
    }

    /// Run the benchmark and return the rich report.
    pub fn run_with_report(
        &self,
        runtime: &BenchmarkConfig<f64>,
    ) -> Result<VorticityStreamCavityReport> {
        let config = self.config.with_runtime_overrides(runtime);
        let (nx, ny) = config.grid_points;
        let grid = StructuredGrid2D::new(nx, ny, 0.0, config.size, 0.0, config.size)?;

        let solver_config = VorticityStreamConfig {
            base: SolverConfig::builder()
                .max_iterations(config.poisson_iterations)
                .tolerance(config.stream_tolerance)
                .time_step(config.time_step)
                .parallel(false)
                .build(),
            time_step: config.time_step,
            stream_tolerance: config.stream_tolerance,
            vorticity_tolerance: config.stream_tolerance,
            sor_omega: 1.85,
        };
        let mut solver = VorticityStreamSolver::new(solver_config, &grid, config.reynolds);
        solver.initialize_lid_driven_cavity(config.lid_velocity)?;

        let reference_profiles = self.reference_profiles(&config);
        let mut residual_history = Vec::with_capacity(config.step_count());
        let mut history = VorticityStreamCavityHistory::default();
        let start = Instant::now();

        history.push(Self::checkpoint(0, 0.0, &solver, &grid)?);

        for step in 1..=config.step_count() {
            let previous_stream = solver.stream_function().clone();
            solver.step()?;
            let residual = Self::stream_residual(solver.stream_function(), &previous_stream);
            residual_history.push(residual);

            if step % config.sample_stride == 0
                || step == config.step_count()
                || residual <= config.stream_tolerance
            {
                history.push(Self::checkpoint(
                    step,
                    step as f64 * config.time_step,
                    &solver,
                    &grid,
                )?);
            }

            if residual <= config.stream_tolerance && step >= 2 {
                break;
            }
        }

        let elapsed = start.elapsed().as_secs_f64();
        let final_checkpoint = history
            .final_checkpoint()
            .ok_or_else(|| Error::InvalidConfiguration("VorticityStreamCavityBenchmark: history must contain at least one checkpoint".into()))?;

        let mut result = BenchmarkResult::new(self.name());
        result.execution_time = elapsed;
        result.convergence.clone_from(&residual_history);
        result.metadata.insert(
            "grid_points".to_string(),
            format!("({}, {})", config.grid_points.0, config.grid_points.1),
        );
        result.metadata.insert(
            "size".to_string(),
            format!("{:.6}", config.size),
        );
        result.metadata.insert(
            "lid_velocity".to_string(),
            format!("{:.6}", config.lid_velocity),
        );
        result.metadata.insert(
            "reynolds".to_string(),
            format!("{:.6}", config.reynolds),
        );
        result.metadata.insert(
            "time_step".to_string(),
            format!("{:.6}", config.time_step),
        );
        result.metadata.insert(
            "max_steps".to_string(),
            format!("{}", config.max_steps),
        );
        result.metadata.insert(
            "sample_stride".to_string(),
            format!("{}", config.sample_stride),
        );
        result.metadata.insert(
            "poisson_iterations".to_string(),
            format!("{}", config.poisson_iterations),
        );

        let (vertical_positions, vertical_reference_values, horizontal_positions, horizontal_reference_values) =
            reference_profiles.as_ref().map_or_else(
                || {
                    (
                        Self::uniform_positions(ny, config.size),
                        Vec::new(),
                        Self::uniform_positions(nx, config.size),
                        Vec::new(),
                    )
                },
                |reference| {
                    (
                        reference.0.iter().map(|(position, _)| *position).collect(),
                        reference.0.iter().map(|(_, value)| *value).collect(),
                        reference.1.iter().map(|(position, _)| *position).collect(),
                        reference.1.iter().map(|(_, value)| *value).collect(),
                    )
                },
            );

        let final_vertical_u = Self::sample_vertical_centerline(
            solver.velocity_field(),
            &grid,
            config.size * 0.5,
            &vertical_positions,
        );
        let final_horizontal_v = Self::sample_horizontal_centerline(
            solver.velocity_field(),
            &grid,
            config.size * 0.5,
            &horizontal_positions,
        );

        let mut reference_errors = Vec::new();
        if !vertical_reference_values.is_empty() && !horizontal_reference_values.is_empty() {
            reference_errors.extend(
                final_vertical_u
                    .iter()
                    .zip(vertical_reference_values.iter())
                    .map(|(sample, expected)| (sample - expected).abs()),
            );
            reference_errors.extend(
                final_horizontal_v
                    .iter()
                    .zip(horizontal_reference_values.iter())
                    .map(|(sample, expected)| (sample - expected).abs()),
            );

            result.metrics.insert(
                "Centerline U RMSE".to_string(),
                Self::rmse(&final_vertical_u, &vertical_reference_values),
            );
            result.metrics.insert(
                "Centerline V RMSE".to_string(),
                Self::rmse(&final_horizontal_v, &horizontal_reference_values),
            );
        } else {
            reference_errors.extend(residual_history.iter().copied());
        }

        result.values.clone_from(&final_vertical_u);
        result.values.extend(final_horizontal_v.iter().copied());
        result.errors = reference_errors;

        result.metrics.insert(
            "Final Residual".to_string(),
            final_checkpoint.residual,
        );
        result.metrics.insert(
            "Final Stream Center".to_string(),
            final_checkpoint.stream_center,
        );
        result.metrics.insert(
            "Final Vorticity Center".to_string(),
            final_checkpoint.vorticity_center,
        );
        result.metrics.insert(
            "Final Max Speed".to_string(),
            final_checkpoint.max_speed,
        );
        result.metrics.insert(
            "Final Divergence Max".to_string(),
            final_checkpoint.max_divergence,
        );
        result.metrics.insert(
            "Final Enstrophy".to_string(),
            final_checkpoint.enstrophy,
        );
        result.metrics.insert(
            "Checkpoint Count".to_string(),
            history.checkpoints.len() as f64,
        );
        result.metrics.insert(
            "Residual Samples".to_string(),
            residual_history.len() as f64,
        );

        Ok(VorticityStreamCavityReport { result, history })
    }

    fn reference_profiles(
        &self,
        config: &VorticityStreamCavityConfig,
    ) -> Option<CenterlineProfiles> {
        if !config.has_ghia_reference() {
            return None;
        }

        let cavity = LidDrivenCavity::new(config.size, config.lid_velocity, config.reynolds);
        let vertical = cavity.ghia_u_centerline(config.reynolds);
        let horizontal = cavity.ghia_v_centerline(config.reynolds);

        if vertical.is_empty() || horizontal.is_empty() {
            return None;
        }

        Some((
            vertical
                .into_iter()
                .map(|(position, value)| (position * config.size, value * config.lid_velocity))
                .collect(),
            horizontal
                .into_iter()
                .map(|(position, value)| (position * config.size, value * config.lid_velocity))
                .collect(),
        ))
    }

    fn checkpoint(
        step: usize,
        time: f64,
        solver: &VorticityStreamSolver<f64>,
        grid: &StructuredGrid2D<f64>,
    ) -> Result<VorticityStreamCavityCheckpoint> {
        Ok(VorticityStreamCavityCheckpoint {
            step,
            time,
            residual: 0.0,
            stream_center: solver.stream_at_center(),
            vorticity_center: solver.vorticity_at_center(),
            max_speed: Self::max_speed(solver.velocity_field()),
            max_divergence: Self::max_divergence(solver.velocity_field(), grid),
            enstrophy: Self::enstrophy(solver.vorticity(), grid),
        })
    }

    fn stream_residual(current: &Array2D<f64>, previous: &Array2D<f64>) -> f64 {
        current
            .iter()
            .zip(previous.iter())
            .map(|(current_value, previous_value)| (current_value - previous_value).abs())
            .fold(0.0, f64::max)
    }

    fn max_speed(field: &Array2D<Vector2<f64>>) -> f64 {
        field
            .iter()
            .map(|velocity| velocity.norm())
            .fold(0.0, f64::max)
    }

    fn max_divergence(field: &Array2D<Vector2<f64>>, grid: &StructuredGrid2D<f64>) -> f64 {
        let two_dx = 2.0 * grid.dx;
        let two_dy = 2.0 * grid.dy;
        let mut max_divergence: f64 = 0.0;

        for i in 1..grid.nx - 1 {
            for j in 1..grid.ny - 1 {
                let du_dx = (field[(i + 1, j)].x - field[(i - 1, j)].x) / two_dx;
                let dv_dy = (field[(i, j + 1)].y - field[(i, j - 1)].y) / two_dy;
                max_divergence = max_divergence.max((du_dx + dv_dy).abs());
            }
        }

        max_divergence
    }

    fn enstrophy(field: &Array2D<f64>, grid: &StructuredGrid2D<f64>) -> f64 {
        let cell_area = grid.dx * grid.dy;
        0.5 * field.iter().map(|vorticity| vorticity * vorticity).sum::<f64>() * cell_area
    }

    fn uniform_positions(count: usize, size: f64) -> Vec<f64> {
        if count <= 1 {
            return vec![0.0];
        }

        let spacing = size / (count.saturating_sub(1) as f64);
        (0..count).map(|index| index as f64 * spacing).collect()
    }

    fn sample_vertical_centerline(
        field: &Array2D<Vector2<f64>>,
        grid: &StructuredGrid2D<f64>,
        x: f64,
        y_positions: &[f64],
    ) -> Vec<f64> {
        y_positions
            .iter()
            .map(|y| Self::sample_velocity(field, grid, x, *y).x)
            .collect()
    }

    fn sample_horizontal_centerline(
        field: &Array2D<Vector2<f64>>,
        grid: &StructuredGrid2D<f64>,
        y: f64,
        x_positions: &[f64],
    ) -> Vec<f64> {
        x_positions
            .iter()
            .map(|x| Self::sample_velocity(field, grid, *x, y).y)
            .collect()
    }

    fn sample_velocity(
        field: &Array2D<Vector2<f64>>,
        grid: &StructuredGrid2D<f64>,
        x: f64,
        y: f64,
    ) -> Vector2<f64> {
        let x_min = grid.bounds.0;
        let x_max = grid.bounds.1;
        let y_min = grid.bounds.2;
        let y_max = grid.bounds.3;
        let clamped_x = x.clamp(x_min, x_max);
        let clamped_y = y.clamp(y_min, y_max);

        let x_position = if grid.dx > 0.0 {
            (clamped_x - x_min) / grid.dx
        } else {
            0.0
        };
        let y_position = if grid.dy > 0.0 {
            (clamped_y - y_min) / grid.dy
        } else {
            0.0
        };

        let i0 = x_position.floor() as usize;
        let j0 = y_position.floor() as usize;
        let i0 = i0.min(grid.nx - 2);
        let j0 = j0.min(grid.ny - 2);
        let i1 = i0 + 1;
        let j1 = j0 + 1;
        let tx = (x_position - i0 as f64).clamp(0.0, 1.0);
        let ty = (y_position - j0 as f64).clamp(0.0, 1.0);

        let v00 = field[(i0, j0)];
        let v10 = field[(i1, j0)];
        let v01 = field[(i0, j1)];
        let v11 = field[(i1, j1)];

        let interpolate = |a: f64, b: f64, c: f64, d: f64| {
            let lower = a * (1.0 - tx) + b * tx;
            let upper = c * (1.0 - tx) + d * tx;
            lower * (1.0 - ty) + upper * ty
        };

        Vector2::new(
            interpolate(v00.x, v10.x, v01.x, v11.x),
            interpolate(v00.y, v10.y, v01.y, v11.y),
        )
    }

    fn rmse(samples: &[f64], reference: &[f64]) -> f64 {
        if samples.is_empty() || reference.is_empty() {
            return 0.0;
        }

        let count = samples.len().min(reference.len());
        let sum_sq = samples
            .iter()
            .zip(reference.iter())
            .take(count)
            .map(|(sample, expected)| (sample - expected).powi(2))
            .sum::<f64>();
        (sum_sq / count as f64).sqrt()
    }

    fn reference_solution_values(&self) -> Option<Vec<f64>> {
        let config = &self.config;
        let reference_profiles = self.reference_profiles(config)?;
        let mut values = reference_profiles
            .0
            .iter()
            .map(|(_, value)| *value)
            .collect::<Vec<_>>();
        values.extend(reference_profiles.1.iter().map(|(_, value)| *value));
        Some(values)
    }

    /// Validate a checkpoint history.
    #[must_use]
    pub fn validate_history(&self, history: &VorticityStreamCavityHistory) -> bool {
        let Some(final_checkpoint) = history.final_checkpoint() else {
            return false;
        };

        history
            .checkpoints
            .iter()
            .all(|checkpoint| {
                checkpoint.residual.is_finite()
                    && checkpoint.stream_center.is_finite()
                    && checkpoint.vorticity_center.is_finite()
                    && checkpoint.max_speed.is_finite()
                    && checkpoint.max_divergence.is_finite()
                    && checkpoint.enstrophy.is_finite()
            })
            && final_checkpoint.max_speed > 0.0
            && final_checkpoint.enstrophy > 0.0
    }
}

impl Default for VorticityStreamCavityBenchmark {
    fn default() -> Self {
        Self::new(VorticityStreamCavityConfig::default())
    }
}

impl Benchmark<f64> for VorticityStreamCavityBenchmark {
    fn name(&self) -> &'static str {
        "2D Vorticity-Stream Cavity"
    }

    fn description(&self) -> &'static str {
        "Lid-driven cavity solved with the vorticity-streamfunction formulation for fast incompressible flow analysis."
    }

    fn run(&self, config: &BenchmarkConfig<f64>) -> Result<BenchmarkResult<f64>> {
        self.run_with_report(config).map(|report| report.result)
    }

    fn reference_solution(&self) -> Option<BenchmarkResult<f64>> {
        let mut result = BenchmarkResult::new(self.name());
        result.values = self.reference_solution_values()?;
        result.metadata.insert(
            "reference".to_string(),
            "Ghia et al. (1982) lid-driven cavity centerline data".to_string(),
        );
        result.metadata.insert(
            "layout".to_string(),
            "vertical centerline U followed by horizontal centerline V".to_string(),
        );
        Some(result)
    }

    fn validate(&self, result: &BenchmarkResult<f64>) -> Result<bool> {
        let Some(final_residual) = result.metrics.get("Final Residual") else {
            return Ok(false);
        };
        let Some(final_enstrophy) = result.metrics.get("Final Enstrophy") else {
            return Ok(false);
        };
        let Some(final_max_speed) = result.metrics.get("Final Max Speed") else {
            return Ok(false);
        };
        let Some(final_divergence) = result.metrics.get("Final Divergence Max") else {
            return Ok(false);
        };

        let finite_history = result.values.iter().all(|value| value.is_finite())
            && result.errors.iter().all(|value| value.is_finite())
            && result.convergence.iter().all(|value| value.is_finite());

        let reference_metrics_ok = result
            .metrics
            .get("Centerline U RMSE")
            .is_none_or(|value| value.is_finite() && *value >= 0.0)
            && result
                .metrics
                .get("Centerline V RMSE")
                .is_none_or(|value| value.is_finite() && *value >= 0.0);

        Ok(finite_history
            && final_residual.is_finite()
            && *final_residual >= 0.0
            && final_enstrophy.is_finite()
            && *final_enstrophy > 0.0
            && final_max_speed.is_finite()
            && *final_max_speed > 0.0
            && final_divergence.is_finite()
            && *final_divergence >= 0.0
            && result.values.len() >= 2
            && !result.convergence.is_empty()
            && reference_metrics_ok)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::benchmarks::BenchmarkRunner;

    #[test]
    fn benchmark_reports_centerline_profiles_and_diagnostics() {
        let benchmark = VorticityStreamCavityBenchmark::default();
        let runtime = BenchmarkConfig {
            resolution: 17,
            tolerance: 0.01,
            max_iterations: 16,
            reynolds_number: 100.0,
            time_step: Some(0.001),
            parallel: false,
        };

        let report = benchmark
            .run_with_report(&runtime)
            .expect("vorticity-stream cavity benchmark should run");

        assert!(!report.history.checkpoints.is_empty());
        assert!(report.history.final_checkpoint().is_some());
        assert!(benchmark.validate_history(&report.history));
        assert!(report.result.metrics.contains_key("Final Residual"));
        assert!(report.result.metrics.contains_key("Final Enstrophy"));
        assert!(report.result.metrics.contains_key("Final Max Speed"));
        assert!(report.result.metrics.contains_key("Final Divergence Max"));
        assert!(report.result.metrics.contains_key("Centerline U RMSE"));
        assert!(report.result.metrics.contains_key("Centerline V RMSE"));
    }

    #[test]
    fn benchmark_trait_adapter_runs_through_validation_framework() {
        let benchmark = VorticityStreamCavityBenchmark::default();
        let runtime = BenchmarkConfig {
            resolution: 17,
            tolerance: 0.01,
            max_iterations: 16,
            reynolds_number: 100.0,
            time_step: Some(0.001),
            parallel: false,
        };

        let result = BenchmarkRunner::run_benchmark(&benchmark, &runtime)
            .expect("benchmark runner should execute the vorticity-stream cavity benchmark");

        assert!(benchmark.validate(&result).expect("validation should be computable"));
        assert!(result.metrics.contains_key("Final Residual"));
        assert!(result.metrics.contains_key("Final Enstrophy"));
        assert!(result.metrics.contains_key("Final Max Speed"));
    }
}