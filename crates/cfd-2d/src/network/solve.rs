use cfd_core::error::Result as CfdResult;
use cfd_core::physics::hemolysis::HemolysisModel;
use nalgebra::RealField;
use num_traits::{Float, FromPrimitive, ToPrimitive};

use super::postprocess::{extract_field_outlet_flow_rate, extract_field_wall_shear};
use super::types::{Channel2dResult, Network2DSolver, Network2dResult};

/// Minimum cross-section area (m²) used to guard against division by zero
/// when computing velocity from flow rate.
const MIN_CROSS_SECTION_AREA_M2: f64 = 1e-18;

impl<T> Network2DSolver<T>
where
    T: RealField + Copy + Float + FromPrimitive + ToPrimitive + std::fmt::Debug,
{
    /// Solve every channel domain in parallel and return per-channel results.
    pub fn solve_all(&mut self, tolerance: f64) -> CfdResult<Network2dResult<T>> {
        use rayon::prelude::*;

        let tol_t = T::from_f64(tolerance).unwrap_or_else(T::one);

        let per_channel: Vec<CfdResult<Channel2dResult<T>>> = self
            .channels
            .par_iter_mut()
            .map(|entry| {
                let area = entry.cross_section_area_m2;
                let u_mean = entry.flow_rate_m3_s / area.max(MIN_CROSS_SECTION_AREA_M2);
                let u_inlet =
                    T::from_f64(u_mean).unwrap_or_else(T::one);

                entry.solver.config.tolerance = tol_t;

                let solve_result = entry.solver.solve(u_inlet).map_err(|e| {
                    cfd_core::error::Error::InvalidInput(format!(
                        "Network2DSolver: channel '{}' failed: {e}",
                        entry.id
                    ))
                })?;

                tracing::debug!(
                    "  [2D] Channel {:<20} | u_in: {:.3e} | inner iter: {} | inner res: {:.3e} | Conv: {}",
                    entry.id,
                    u_inlet.to_f64().unwrap_or(0.0),
                    solve_result.iterations,
                    solve_result.residual.to_f64().unwrap_or(0.0),
                    solve_result.converged
                );
                let shear_rate = entry.cross_section.wall_shear_rate(u_mean);
                let shear_pa = entry.viscosity_pa_s * shear_rate;
                let (field_max, field_mean) = extract_field_wall_shear(&entry.solver);

                let field_outlet_flow = extract_field_outlet_flow_rate(
                    &entry.solver,
                    T::from_f64(entry.cross_section_area_m2).unwrap_or_else(T::zero),
                );
                let outlet_flow_error = field_outlet_flow - entry.reference_trace.flow_rate_m3_s;
                let outlet_flow_error_pct = if Float::abs(entry.reference_trace.flow_rate_m3_s)
                    > T::from_f64(MIN_CROSS_SECTION_AREA_M2).unwrap_or_else(T::zero)
                {
                    Float::abs(outlet_flow_error / entry.reference_trace.flow_rate_m3_s)
                        * T::from_f64(100.0).unwrap_or_else(T::one)
                } else {
                    T::zero()
                };

                let t_s = area * entry.length_m / entry.flow_rate_m3_s.max(1e-30);
                let hi = HemolysisModel::giersiepen_millifluidic()
                    .damage_index(shear_pa, t_s)
                    .unwrap_or(0.0);
                let hi_t = T::from_f64(hi).unwrap_or_else(T::zero);

                Ok(Channel2dResult {
                    channel_id: entry.id.clone(),
                    therapy_zone: entry.therapy_zone,
                    is_venturi_throat: entry.is_venturi_throat,
                    solve_result,
                    wall_shear_pa: T::from_f64(shear_pa).unwrap_or_else(T::zero),
                    field_wall_shear_max_pa: field_max,
                    field_wall_shear_mean_pa: field_mean,
                    field_outlet_flow_m3_s: field_outlet_flow,
                    field_outlet_flow_error_m3_s: outlet_flow_error,
                    field_outlet_flow_error_pct: outlet_flow_error_pct,
                    transit_time_s: T::from_f64(t_s).unwrap_or_else(T::zero),
                    hemolysis_index: hi_t,
                    reference_trace: entry.reference_trace.clone(),
                })
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
            total_outlet_error_pct / T::from_usize(results.len()).unwrap_or_else(T::one)
        };

        Ok(Network2dResult {
            channels: results,
            total_hemolysis_index: total_hi,
            converged_count,
            max_field_outlet_flow_error_pct: max_outlet_error_pct,
            mean_field_outlet_flow_error_pct: mean_outlet_error_pct,
            reference_trace: self.reference_trace.clone(),
        })
    }
}
