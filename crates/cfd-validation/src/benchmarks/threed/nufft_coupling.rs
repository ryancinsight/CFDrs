//! 3D NUFFT coupling benchmark for irregular probes and Lagrangian marker spreading.

use super::super::{Benchmark, BenchmarkConfig, BenchmarkResult};
use cfd_3d::{LagrangianPoint, NufftMarkerCoupler3D};
use cfd_core::error::{Error, Result};
use cfd_core::physics::fluid_dynamics::VelocityField;
use nalgebra::Vector3;
use serde::{Deserialize, Serialize};
use std::collections::HashMap;
use std::time::Instant;

/// Configuration for the NUFFT coupling benchmark.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct NufftCouplingBenchmarkConfig {
    /// Grid dimensions `(nx, ny, nz)`.
    pub dimensions: (usize, usize, usize),
    /// Characteristic physical length used for the periodic domain.
    pub length_scale: f64,
    /// Number of irregular probe locations.
    pub probe_count: usize,
    /// Number of Lagrangian markers.
    pub marker_count: usize,
    /// Weight assigned to each marker.
    pub marker_weight: f64,
    /// Scale factor for the synthetic marker force field.
    pub marker_force_scale: f64,
    /// Validation tolerance for probe interpolation and force balance.
    pub validation_tolerance: f64,
}

impl NufftCouplingBenchmarkConfig {
    /// Create a validated NUFFT coupling benchmark configuration.
    pub fn new(
        dimensions: (usize, usize, usize),
        length_scale: f64,
        probe_count: usize,
        marker_count: usize,
        marker_weight: f64,
        marker_force_scale: f64,
        validation_tolerance: f64,
    ) -> Result<Self> {
        let (nx, ny, nz) = dimensions;
        if nx == 0 || ny == 0 || nz == 0 {
            return Err(Error::InvalidConfiguration(
                "NufftCouplingBenchmarkConfig: dimensions must be greater than zero".into(),
            ));
        }
        if !length_scale.is_finite() || length_scale <= 0.0 {
            return Err(Error::InvalidConfiguration(
                "NufftCouplingBenchmarkConfig: length_scale must be finite and positive".into(),
            ));
        }
        if probe_count == 0 {
            return Err(Error::InvalidConfiguration(
                "NufftCouplingBenchmarkConfig: probe_count must be greater than zero".into(),
            ));
        }
        if marker_count == 0 {
            return Err(Error::InvalidConfiguration(
                "NufftCouplingBenchmarkConfig: marker_count must be greater than zero".into(),
            ));
        }
        if !marker_weight.is_finite() || marker_weight <= 0.0 {
            return Err(Error::InvalidConfiguration(
                "NufftCouplingBenchmarkConfig: marker_weight must be finite and positive".into(),
            ));
        }
        if !marker_force_scale.is_finite() || marker_force_scale < 0.0 {
            return Err(Error::InvalidConfiguration(
                "NufftCouplingBenchmarkConfig: marker_force_scale must be finite and non-negative".into(),
            ));
        }
        if !validation_tolerance.is_finite() || validation_tolerance < 0.0 {
            return Err(Error::InvalidConfiguration(
                "NufftCouplingBenchmarkConfig: validation_tolerance must be finite and non-negative".into(),
            ));
        }

        Ok(Self {
            dimensions,
            length_scale,
            probe_count,
            marker_count,
            marker_weight,
            marker_force_scale,
            validation_tolerance,
        })
    }

    /// Return the periodic cubic domain lengths.
    #[must_use]
    pub fn domain_lengths(&self) -> (f64, f64, f64) {
        (self.length_scale, self.length_scale, self.length_scale)
    }

    fn with_runtime_overrides(&self, runtime: &BenchmarkConfig<f64>) -> Self {
        let mut resolved = self.clone();

        if runtime.resolution > 0 {
            resolved.dimensions = (runtime.resolution, runtime.resolution, runtime.resolution);
        }

        if runtime.tolerance.is_finite() && runtime.tolerance >= 0.0 {
            resolved.validation_tolerance = runtime.tolerance;
        }

        let budget = runtime.max_iterations.max(1);
        resolved.probe_count = resolved.probe_count.min(budget);
        resolved.marker_count = resolved.marker_count.min(budget);

        resolved
    }
}

impl Default for NufftCouplingBenchmarkConfig {
    fn default() -> Self {
        Self {
            dimensions: (16, 16, 16),
            length_scale: 1.0,
            probe_count: 12,
            marker_count: 4,
            marker_weight: 1.0,
            marker_force_scale: 0.5,
            validation_tolerance: 1e-9,
        }
    }
}

/// One irregular-probe checkpoint.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct NufftCouplingCheckpoint {
    /// Probe index.
    pub probe_index: usize,
    /// Probe position.
    pub position: [f64; 3],
    /// Analytic velocity at the probe.
    pub analytic_velocity: [f64; 3],
    /// NUFFT-sampled velocity at the probe.
    pub sampled_velocity: [f64; 3],
    /// Absolute probe interpolation error.
    pub absolute_error: f64,
}

/// History collected during the NUFFT coupling benchmark.
#[derive(Debug, Clone, Default, Serialize, Deserialize)]
pub struct NufftCouplingBenchmarkHistory {
    /// Irregular probe checkpoints.
    pub probe_checkpoints: Vec<NufftCouplingCheckpoint>,
    /// Force-balance residual after marker spreading.
    pub force_balance_residual: [f64; 3],
}

/// Combined benchmark report containing the generic result plus rich history.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct NufftCouplingBenchmarkReport {
    /// Generic validation result expected by the benchmark framework.
    pub result: BenchmarkResult<f64>,
    /// NUFFT coupling history.
    pub history: NufftCouplingBenchmarkHistory,
}

/// NUFFT coupling benchmark driver.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct NufftCouplingBenchmark3D {
    /// Static benchmark configuration.
    pub config: NufftCouplingBenchmarkConfig,
}

impl NufftCouplingBenchmark3D {
    /// Create a new NUFFT coupling benchmark.
    #[must_use]
    pub fn new(config: NufftCouplingBenchmarkConfig) -> Self {
        Self { config }
    }

    /// Run the benchmark and return the rich report.
    pub fn run_with_report(
        &self,
        runtime: &BenchmarkConfig<f64>,
    ) -> Result<NufftCouplingBenchmarkReport> {
        let config = self.config.with_runtime_overrides(runtime);
        let coupler = NufftMarkerCoupler3D::from_dimensions(
            config.dimensions,
            config.domain_lengths(),
        )?;

        let field = Self::build_velocity_field(&coupler, config.length_scale);
        let probes = Self::probe_positions(config.probe_count, config.domain_lengths());
        let markers = Self::marker_points(
            config.marker_count,
            config.domain_lengths(),
            config.marker_weight,
            config.marker_force_scale,
        );

        let start = Instant::now();
        let sampled = coupler.sample_velocity_at_positions(&field, &probes)?;
        let spread = coupler.spread_marker_forces(&markers)?;
        let elapsed = start.elapsed().as_secs_f64();

        let mut history = NufftCouplingBenchmarkHistory::default();
        let mut errors = Vec::with_capacity(probes.len());
        let mut convergence = Vec::with_capacity(probes.len());
        let mut cumulative_error = 0.0_f64;

        for (index, (position, sample)) in probes.iter().zip(sampled.iter()).enumerate() {
            let analytic = Self::analytic_velocity(position, config.length_scale);
            let absolute_error = (sample - analytic).norm();
            cumulative_error = cumulative_error.max(absolute_error);

            history.probe_checkpoints.push(NufftCouplingCheckpoint {
                probe_index: index,
                position: [position.x, position.y, position.z],
                analytic_velocity: [analytic.x, analytic.y, analytic.z],
                sampled_velocity: [sample.x, sample.y, sample.z],
                absolute_error,
            });
            errors.push(absolute_error);
            convergence.push(cumulative_error);
        }

        let total_weighted_force = Self::total_weighted_force(&markers);
        let spread_sum = Self::velocity_field_sum(&spread);
        let residual = spread_sum - total_weighted_force;
        history.force_balance_residual = [residual.x, residual.y, residual.z];

        let max_probe_error = errors.iter().copied().fold(0.0_f64, f64::max);
        let mean_probe_error = if errors.is_empty() {
            0.0
        } else {
            errors.iter().sum::<f64>() / errors.len() as f64
        };
        let force_balance_linf = residual.x.abs().max(residual.y.abs()).max(residual.z.abs());
        let force_balance_l1 = residual.x.abs() + residual.y.abs() + residual.z.abs();

        let mut result = BenchmarkResult::new(self.name());
        result.values.clone_from(&errors);
        result.errors = errors;
        result.convergence = convergence;
        result.execution_time = elapsed;
        result.metrics = HashMap::from([
            ("probe_count".to_string(), config.probe_count as f64),
            ("marker_count".to_string(), config.marker_count as f64),
            ("max_probe_error".to_string(), max_probe_error),
            ("mean_probe_error".to_string(), mean_probe_error),
            ("force_balance_linf".to_string(), force_balance_linf),
            ("force_balance_l1".to_string(), force_balance_l1),
        ]);
        result.metadata = HashMap::from([
            ("dimensions".to_string(), format!("{:?}", config.dimensions)),
            ("length_scale".to_string(), config.length_scale.to_string()),
            ("marker_weight".to_string(), config.marker_weight.to_string()),
            ("marker_force_scale".to_string(), config.marker_force_scale.to_string()),
        ]);

        Ok(NufftCouplingBenchmarkReport { result, history })
    }

    fn build_velocity_field(
        coupler: &NufftMarkerCoupler3D,
        length_scale: f64,
    ) -> VelocityField<f64> {
        let grid = coupler.grid();
        let mut components = Vec::with_capacity(grid.nx * grid.ny * grid.nz);

        for k in 0..grid.nz {
            let z = k as f64 * grid.dz;
            for j in 0..grid.ny {
                let y = j as f64 * grid.dy;
                for i in 0..grid.nx {
                    let x = i as f64 * grid.dx;
                    components.push(Self::analytic_velocity(&Vector3::new(x, y, z), length_scale));
                }
            }
        }

        VelocityField {
            components,
            dimensions: (grid.nx, grid.ny, grid.nz),
        }
    }

    fn analytic_velocity(position: &Vector3<f64>, length_scale: f64) -> Vector3<f64> {
        let phase_x = std::f64::consts::TAU * position.x / length_scale;
        let phase_y = std::f64::consts::TAU * position.y / length_scale;
        let phase_z = std::f64::consts::TAU * position.z / length_scale;

        Vector3::new(
            phase_x.sin() * phase_y.cos() * phase_z.cos(),
            phase_x.cos() * phase_y.sin() * phase_z.cos(),
            phase_x.cos() * phase_y.cos() * phase_z.sin(),
        )
    }

    fn probe_positions(count: usize, lengths: (f64, f64, f64)) -> Vec<Vector3<f64>> {
        let (lx, ly, lz) = lengths;
        (0..count)
            .map(|index| {
                let offset = index as f64 + 0.5;
                let x = (offset * 0.381_966_011_250_105_1).fract() * lx;
                let y = (offset * 0.527_864_045_000_420_6 + 0.125).fract() * ly;
                let z = (offset * 0.618_033_988_749_894_8 + 0.25).fract() * lz;
                Vector3::new(x, y, z)
            })
            .collect()
    }

    fn marker_points(
        count: usize,
        lengths: (f64, f64, f64),
        marker_weight: f64,
        marker_force_scale: f64,
    ) -> Vec<LagrangianPoint<f64>> {
        Self::probe_positions(count, lengths)
            .into_iter()
            .enumerate()
            .map(|(index, position)| {
                let phase_x = std::f64::consts::TAU * position.x / lengths.0;
                let phase_y = std::f64::consts::TAU * position.y / lengths.1;
                let phase_z = std::f64::consts::TAU * position.z / lengths.2;
                let force = Vector3::new(
                    marker_force_scale * (phase_x.cos() + 0.25 * phase_y.sin()),
                    marker_force_scale * (phase_y.cos() - 0.25 * phase_z.sin()),
                    marker_force_scale * (phase_z.cos() + 0.25 * phase_x.sin()),
                );

                let mut point = LagrangianPoint::new(position, marker_weight + index as f64 * 0.0);
                point.force = force;
                point
            })
            .collect()
    }

    fn total_weighted_force(markers: &[LagrangianPoint<f64>]) -> Vector3<f64> {
        markers.iter().fold(Vector3::zeros(), |acc, marker| {
            acc + marker.force * marker.weight
        })
    }

    fn velocity_field_sum(field: &VelocityField<f64>) -> Vector3<f64> {
        field
            .components
            .iter()
            .fold(Vector3::zeros(), |acc, value| acc + value)
    }
}

impl Benchmark<f64> for NufftCouplingBenchmark3D {
    fn name(&self) -> &'static str {
        "nufft_coupling_3d"
    }

    fn description(&self) -> &'static str {
        "NUFFT-backed irregular probe sampling and Lagrangian marker spreading"
    }

    fn run(&self, config: &BenchmarkConfig<f64>) -> Result<BenchmarkResult<f64>> {
        self.run_with_report(config).map(|report| report.result)
    }

    fn reference_solution(&self) -> Option<BenchmarkResult<f64>> {
        None
    }

    fn validate(&self, result: &BenchmarkResult<f64>) -> Result<bool> {
        let max_probe_error = result
            .metrics
            .get("max_probe_error")
            .copied()
            .unwrap_or(f64::INFINITY);
        let force_balance_linf = result
            .metrics
            .get("force_balance_linf")
            .copied()
            .unwrap_or(f64::INFINITY);

        Ok(max_probe_error <= self.config.validation_tolerance
            && force_balance_linf <= self.config.validation_tolerance)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn benchmark_reports_probe_sampling_and_force_balance() {
        let benchmark = NufftCouplingBenchmark3D::new(NufftCouplingBenchmarkConfig::default());
        let runtime = BenchmarkConfig {
            resolution: 8,
            tolerance: 1e-8,
            max_iterations: 16,
            reynolds_number: 100.0,
            time_step: None,
            parallel: false,
        };

        let report = benchmark
            .run_with_report(&runtime)
            .expect("benchmark should run");

        assert_eq!(report.history.probe_checkpoints.len(), 12);
        assert!(report.result.metrics["max_probe_error"] <= 1e-8);
        assert!(report.result.metrics["force_balance_linf"] <= 1e-8);
        assert!(report
            .history
            .probe_checkpoints
            .iter()
            .all(|checkpoint| checkpoint.absolute_error.is_finite()));
    }

    #[test]
    fn benchmark_trait_adapter_runs_through_validation_framework() {
        let benchmark = NufftCouplingBenchmark3D::new(NufftCouplingBenchmarkConfig::default());
        let runtime = BenchmarkConfig::default();
        let result = benchmark.run(&runtime).expect("benchmark should run");

        assert!(benchmark.validate(&result).expect("validation should succeed"));
    }
}