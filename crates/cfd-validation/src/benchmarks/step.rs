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

    /// Return the published Armaly reference for a supported benchmark case.
    ///
    /// The reference is keyed by the configured Reynolds number and the
    /// expansion geometry. Armaly et al. show that the laminar reattachment
    /// curve is nonlinear and depends on more than Reynolds number, so this
    /// adapter exposes only the committed anchor cases instead of inventing an
    /// interpolation between them. The source is Armaly et al., *Experimental
    /// and theoretical investigation of backward-facing step flow*, Journal of
    /// Fluid Mechanics 127 (1983), Fig. 4 and the corresponding two-dimensional
    /// comparison data: <https://courses.washington.edu/me431/handouts/armaly-jfm-83.pdf>.
    #[must_use]
    pub fn reference_solution_for(
        &self,
        config: &BenchmarkConfig<T>,
    ) -> Option<BenchmarkResult<T>> {
        self.reference_solution_for_reynolds(config.reynolds_number)
    }

    fn reference_solution_for_reynolds(&self, reynolds_number: T) -> Option<BenchmarkResult<T>> {
        let two = <T as FloatElement>::from_f64(2.0);
        if self.channel_height != two * self.step_height {
            return None;
        }

        let (reference_ratio, experimental_spread, source_case) =
            if reynolds_number == <T as FloatElement>::from_f64(100.0) {
                (
                    <T as FloatElement>::from_f64(2.84),
                    <T as FloatElement>::from_f64(0.13),
                    "Armaly 2D Re=100",
                )
            } else if reynolds_number == <T as FloatElement>::from_f64(389.0) {
                (
                    <T as FloatElement>::from_f64(7.83),
                    <T as FloatElement>::from_f64(0.42),
                    "Armaly 2D Re=389",
                )
            } else {
                return None;
            };

        let reference_reattachment = reference_ratio * self.step_height;
        let reference_error = experimental_spread * self.step_height;
        let mut metrics = std::collections::HashMap::new();
        metrics.insert("reynolds_number".to_string(), reynolds_number);
        metrics.insert("step_height".to_string(), self.step_height);
        metrics.insert("channel_height".to_string(), self.channel_height);
        let mut metadata = std::collections::HashMap::new();
        metadata.insert("reference_case".to_string(), source_case.to_string());
        metadata.insert(
            "reference_source".to_string(),
            "Armaly et al. (1983), JFM 127, Fig. 4".to_string(),
        );

        Some(BenchmarkResult {
            name: "Backward Facing Step (Armaly Reference)".to_string(),
            values: vec![reference_reattachment],
            errors: vec![reference_error],
            convergence: vec![],
            execution_time: 0.0,
            metrics,
            metadata,
        })
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
        simple.pressure_sweep_cap = Some(33);
        simple.pressure_sor_relaxation = Some(<T as FloatElement>::from_f64(1.7));
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
            metrics: {
                let mut metrics = std::collections::HashMap::new();
                metrics.insert("reynolds_number".to_string(), config.reynolds_number);
                metrics.insert("step_height".to_string(), self.step_height);
                metrics.insert("channel_height".to_string(), self.channel_height);
                metrics
            },
            metadata,
        })
    }

    fn reference_solution(&self) -> Option<BenchmarkResult<T>> {
        self.reference_solution_for(&BenchmarkConfig::default())
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

        let Some(&reynolds_number) = result.metrics.get("reynolds_number") else {
            return Err(Error::InvalidInput(
                "backward-facing-step result is missing reynolds_number metadata".to_string(),
            ));
        };
        let reference = self
            .reference_solution_for_reynolds(reynolds_number)
            .ok_or_else(|| {
                Error::InvalidInput("backward-facing-step reference is unavailable".to_string())
            })?;
        let Some(&reference_reattachment) = reference.values.first() else {
            return Err(Error::InvalidInput(
                "backward-facing-step reference has no reattachment value".to_string(),
            ));
        };
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
