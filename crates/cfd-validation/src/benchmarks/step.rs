//! Backward-facing-step benchmark adapter.
//!
//! Geometry, SIMPLE execution, and signed wall-shear extraction are owned by
//! `cfd-2d`. This module maps the validation framework's common configuration
//! and result types onto that provider contract.

use super::{Benchmark, BenchmarkConfig, BenchmarkResult};
use cfd_core::error::{Error, Result};
use cfd_core::CfdScalar;
use eunomia::{FloatElement, RealField};

/// Backward-facing-step benchmark geometry.
pub struct BackwardFacingStep<T: RealField + Copy> {
    /// Height of the upstream solid step.
    pub step_height: T,
    /// Height of the downstream channel.
    pub channel_height: T,
    /// Fluid length before the step face.
    pub upstream_length: T,
    /// Fluid length after the step face.
    pub downstream_length: T,
    /// Total streamwise channel length.
    pub channel_length: T,
    /// Uniform inlet velocity.
    pub inlet_velocity: T,
}

impl<T: CfdScalar + Copy + FloatElement> BackwardFacingStep<T> {
    /// Create a backward-facing-step benchmark with explicit upstream and downstream lengths.
    pub fn new(
        step_height: T,
        channel_height: T,
        upstream_length: T,
        downstream_length: T,
        inlet_velocity: T,
    ) -> Self {
        Self {
            step_height,
            channel_height,
            upstream_length,
            downstream_length,
            channel_length: upstream_length + downstream_length,
            inlet_velocity,
        }
    }
}

impl<T: CfdScalar + Copy + FloatElement> Benchmark<T> for BackwardFacingStep<T> {
    fn name(&self) -> &'static str {
        "Backward Facing Step"
    }

    fn description(&self) -> &'static str {
        "2D laminar flow over a backward-facing step"
    }

    fn run(&self, config: &BenchmarkConfig<T>) -> Result<BenchmarkResult<T>> {
        let start = std::time::Instant::now();
        let geometry = cfd_2d::solvers::ns_fvm::BackwardFacingStepGeometry::new(
            self.step_height,
            self.channel_height,
            self.upstream_length,
            self.downstream_length,
            self.inlet_velocity,
        )?;
        let mut simple = cfd_2d::solvers::ns_fvm::SIMPLEConfig::default();
        simple.max_iterations = config.max_iterations;
        simple.tolerance = config.tolerance;
        let provider_config = cfd_2d::solvers::ns_fvm::BackwardFacingStepConfig {
            resolution: config.resolution,
            reynolds_number: config.reynolds_number,
            simple,
        };
        let solved = cfd_2d::solvers::ns_fvm::BackwardFacingStepSolver::new(geometry)
            .solve(&provider_config)?;

        let mut metadata = std::collections::HashMap::new();
        metadata.insert(
            "boundary_contract".to_string(),
            "masked upstream step, no-slip walls, parabolic inlet, fixed-pressure outlet"
                .to_string(),
        );
        metadata.insert(
            "wall_shear_samples".to_string(),
            solved.wall_shear.len().to_string(),
        );

        Ok(BenchmarkResult {
            name: self.name().to_string(),
            values: vec![solved.reattachment_length],
            errors: vec![],
            convergence: vec![solved.solve.residual],
            execution_time: start.elapsed().as_secs_f64(),
            metrics: std::collections::HashMap::new(),
            metadata,
        })
    }

    fn reference_solution(&self) -> Option<BenchmarkResult<T>> {
        // The default benchmark configuration uses Re_h = 100. For the
        // expansion-ratio-two geometry, the published Re_h = 100 reference is
        // x_r / h ≈ 3.0; runtime results never substitute this value for a
        // solved wall-shear crossing. See Ziegner (2014), Table 5.1, p. 52:
        // https://www5.in.tum.de/pub/ziegner_14.pdf
        let reference_reattachment = <T as FloatElement>::from_f64(3.0) * self.step_height;
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
        let Some(&computed_reattachment) = result.values.first() else {
            return Ok(false);
        };
        let Some(&last_residual) = result.convergence.last() else {
            return Ok(false);
        };
        if !computed_reattachment.is_finite() || !last_residual.is_finite() {
            return Ok(false);
        }

        let reference = self.reference_solution().ok_or_else(|| {
            Error::InvalidInput("backward-facing-step reference is unavailable".to_string())
        })?;
        let reference_reattachment = reference.values[0];
        let tolerance = <T as FloatElement>::from_f64(0.30);
        let relative_error =
            (computed_reattachment - reference_reattachment).abs() / reference_reattachment;
        let within_reference = relative_error <= tolerance;
        let physically_reasonable = computed_reattachment > <T as FloatElement>::from_f64(0.0)
            && computed_reattachment < <T as FloatElement>::from_f64(20.0) * self.step_height;
        let converged = last_residual <= <T as FloatElement>::from_f64(1e-4);
        Ok(within_reference && physically_reasonable && converged)
    }
}
