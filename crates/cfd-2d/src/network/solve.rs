use cfd_core::error::Result as CfdResult;
use nalgebra::RealField;
use num_traits::{Float, FromPrimitive, ToPrimitive};

use super::channel::solve_channel_entry;
use super::types::{Channel2dResult, Network2DSolver, Network2dResult, ProjectedNetwork2dResult};

impl<T> Network2DSolver<T>
where
    T: RealField + Copy + Float + FromPrimitive + ToPrimitive + std::fmt::Debug,
{
    /// Solve every channel domain in parallel and return per-channel results.
    pub fn solve_all(&mut self, tolerance: f64) -> CfdResult<Network2dResult<T>> {
        use rayon::prelude::*;

        let reference_trace = self.reference_trace.clone();
        let separation_tracking_enabled = self.separation_tracking_enabled;

        let per_channel: Vec<CfdResult<Channel2dResult<T>>> = self
            .channels
            .par_iter_mut()
            .zip(reference_trace.channel_traces.par_iter())
            .map(|(entry, channel_trace)| {
                solve_channel_entry(
                    entry,
                    tolerance,
                    channel_trace.flow_rate_m3_s.to_f64().unwrap_or(0.0),
                    channel_trace,
                    separation_tracking_enabled,
                )
            })
            .collect();

        let mut results = Vec::with_capacity(per_channel.len());
        let mut total_hi = T::zero();
        let mut total_outlet_error_pct = T::zero();
        let mut max_outlet_error_pct = T::zero();
        let mut converged_count = 0usize;

        for result in per_channel {
            let channel = result?;
            if channel.solve_result.converged {
                converged_count += 1;
            }
            if channel.field_outlet_flow_error_pct > max_outlet_error_pct {
                max_outlet_error_pct = channel.field_outlet_flow_error_pct;
            }
            total_outlet_error_pct += channel.field_outlet_flow_error_pct;
            total_hi += channel.hemolysis_index;
            results.push(channel);
        }

        let mean_outlet_error_pct = if results.is_empty() {
            T::zero()
        } else {
            total_outlet_error_pct
                / T::from_usize(results.len()).expect("analytical constant conversion")
        };

        Ok(Network2dResult {
            channels: results,
            total_hemolysis_index: total_hi,
            converged_count,
            max_field_outlet_flow_error_pct: max_outlet_error_pct,
            mean_field_outlet_flow_error_pct: mean_outlet_error_pct,
            reference_trace,
        })
    }

    /// Solve the projected schematics-driven network and return the compatibility
    /// result plus projection metadata.
    pub fn solve_projected(
        &mut self,
        tolerance: f64,
    ) -> CfdResult<ProjectedNetwork2dResult<T>> {
        let result = self.solve_all(tolerance)?;
        let projection = self.projection_summary_ref().clone();
        Ok(ProjectedNetwork2dResult { result, projection })
    }
}
