use cfd_core::error::Result as CfdResult;
use cfd_core::physics::hemolysis::HemolysisModel;
use cfd_schematics::domain::model::CrossSectionSpec;
use nalgebra::RealField;
use num_traits::{Float, FromPrimitive, ToPrimitive};

use crate::solvers::cell_tracking::physics::{CellTrackerConfig, VelocityFieldInterpolator};
use crate::solvers::cell_tracking::population::{CellPopulation, TrackedCell};
use crate::solvers::cell_tracking::tracker::CellTracker;
use crate::solvers::ns_fvm::SolveResult;

use super::postprocess::{
    extract_field_inlet_outlet_pressure, extract_field_outlet_flow_rate, extract_field_wall_shear,
};
use super::types::{Channel2dEntry, Channel2dResult};
use super::ChannelReferenceTrace;

const MIN_CROSS_SECTION_AREA_M2: f64 = 1e-18;
const MIN_FLOW_RATE_M3_S: f64 = 1e-30;

pub(crate) fn solve_channel_entry<T>(
    entry: &mut Channel2dEntry<T>,
    tolerance: f64,
    flow_rate_m3_s: f64,
    reference_trace: &ChannelReferenceTrace<T>,
    separation_tracking_enabled: bool,
) -> CfdResult<Channel2dResult<T>>
where
    T: RealField + Copy + Float + FromPrimitive + ToPrimitive + std::fmt::Debug,
{
    let tol_t = T::from_f64(tolerance).expect("analytical constant conversion");
    entry.flow_rate_m3_s = flow_rate_m3_s;
    entry.reference_trace = reference_trace.clone();

    let area = entry.cross_section_area_m2;
    let area_safe = area.max(MIN_CROSS_SECTION_AREA_M2);
    let u_mean = flow_rate_m3_s / area_safe;
    let u_inlet = T::from_f64(u_mean).expect("analytical constant conversion");

    entry.solver.config.tolerance = tol_t;

    let shear_rate = entry.cross_section.wall_shear_rate(u_mean);
    let shear_pa = entry.viscosity_pa_s * shear_rate;
    let solve_attempt = entry.solver.solve(u_inlet);

    let (
        field_wall_shear_max_pa,
        field_wall_shear_mean_pa,
        field_inlet_pressure_pa,
        field_outlet_pressure_pa,
        field_pressure_drop_pa,
        field_effective_resistance_pa_s_per_m3,
        field_outlet_flow_m3_s,
        field_outlet_flow_error_m3_s,
        field_outlet_flow_error_pct,
        field_separation_efficiency_pct,
        solve_result,
    ) = match solve_attempt {
        Ok(solve_result) => {
            if solve_result.converged {
                tracing::debug!(
                    "  [2D] Channel {:<20} | u_in: {:.3e} | inner iter: {} | inner res: {:.3e} | Conv: {}",
                    entry.id,
                    u_inlet.to_f64().unwrap_or(0.0),
                    solve_result.iterations,
                    solve_result.residual.to_f64().unwrap_or(0.0),
                    solve_result.converged
                );
            } else {
                tracing::warn!(
                    channel = %entry.id,
                    residual = solve_result.residual.to_f64().unwrap_or(0.0),
                    "2D solve did not converge; preserving extracted field metrics"
                );
            }

            let (field_wall_shear_max_pa, field_wall_shear_mean_pa) =
                extract_field_wall_shear(&entry.solver);
            let (field_inlet_pressure_pa, field_outlet_pressure_pa) =
                extract_field_inlet_outlet_pressure(&entry.solver);
            let field_pressure_drop_pa = field_inlet_pressure_pa - field_outlet_pressure_pa;
            let field_outlet_flow_m3_s = extract_field_outlet_flow_rate(
                &entry.solver,
                T::from_f64(area).expect("analytical constant conversion"),
            );
            let field_outlet_flow_error_m3_s =
                field_outlet_flow_m3_s - reference_trace.flow_rate_m3_s;
            let field_outlet_flow_error_pct = if Float::abs(reference_trace.flow_rate_m3_s)
                > T::from_f64(MIN_FLOW_RATE_M3_S).expect("analytical constant conversion")
            {
                Float::abs(field_outlet_flow_error_m3_s / reference_trace.flow_rate_m3_s)
                    * T::from_f64(100.0).expect("analytical constant conversion")
            } else {
                T::zero()
            };

            let flow_abs =
                Float::abs(T::from_f64(flow_rate_m3_s).expect("analytical constant conversion"));
            let mut field_effective_resistance_pa_s_per_m3 = reference_trace.resistance_pa_s_per_m3;
            if flow_abs > T::from_f64(MIN_FLOW_RATE_M3_S).expect("analytical constant conversion") {
                field_effective_resistance_pa_s_per_m3 =
                    Float::abs(field_pressure_drop_pa) / flow_abs;
            }

            let mut field_separation_efficiency_pct = None;
            if separation_tracking_enabled
                && matches!(entry.cross_section, CrossSectionSpec::Rectangular { .. })
            {
                let cnf = CellTrackerConfig {
                    viscosity: entry.viscosity_pa_s,
                    fluid_density: 1060.0,
                    hydraulic_diameter_m: entry.cross_section.hydraulic_diameter(),
                    u_max: u_mean.abs(),
                    ..Default::default()
                };
                let tracker = CellTracker::new(&entry.solver, cnf);

                let n_per_pop = 50usize;
                let mut cells = Vec::with_capacity(n_per_pop * 2);
                let (x_min, _, y_min, y_max) = entry.solver.bounds();
                let h_bounds = y_max - y_min;

                for i in 0..n_per_pop {
                    let y = y_min
                        + h_bounds * 0.05
                        + h_bounds * 0.9 * (i as f64 / (n_per_pop - 1) as f64);
                    cells.push(TrackedCell {
                        population: CellPopulation::CTC,
                        x: x_min + 1e-6,
                        y,
                        vx: u_mean.abs(),
                        vy: 0.0,
                        id: i,
                    });
                    cells.push(TrackedCell {
                        population: CellPopulation::RBC,
                        x: x_min + 1e-6,
                        y,
                        vx: u_mean.abs(),
                        vy: 0.0,
                        id: n_per_pop + i,
                    });
                }

                let trajectories = tracker.trace_cells(&cells, 1e-5, 5_000_000);
                let mut sum_y_ctc = 0.0;
                let mut count_ctc = 0usize;
                let mut sum_y_rbc = 0.0;
                let mut count_rbc = 0usize;

                for traj in &trajectories {
                    if let Some(pos) = traj.positions.last() {
                        let y_tilde = (pos[1] - y_min) / h_bounds.max(1e-12);
                        match traj.population {
                            CellPopulation::CTC => {
                                sum_y_ctc += y_tilde;
                                count_ctc += 1;
                            }
                            CellPopulation::RBC | CellPopulation::WBC => {
                                sum_y_rbc += y_tilde;
                                count_rbc += 1;
                            }
                        }
                    }
                }

                let y_tilde_ctc = if count_ctc > 0 {
                    sum_y_ctc / f64::from(count_ctc as u32)
                } else {
                    0.5
                };
                let y_tilde_rbc = if count_rbc > 0 {
                    sum_y_rbc / f64::from(count_rbc as u32)
                } else {
                    0.5
                };
                let diff_pct = (y_tilde_ctc - y_tilde_rbc).abs() * 100.0;
                field_separation_efficiency_pct =
                    Some(T::from_f64(diff_pct).unwrap_or_else(T::zero));
            }

            (
                field_wall_shear_max_pa,
                field_wall_shear_mean_pa,
                field_inlet_pressure_pa,
                field_outlet_pressure_pa,
                field_pressure_drop_pa,
                field_effective_resistance_pa_s_per_m3,
                field_outlet_flow_m3_s,
                field_outlet_flow_error_m3_s,
                field_outlet_flow_error_pct,
                field_separation_efficiency_pct,
                solve_result,
            )
        }
        Err(error) => {
            tracing::warn!(
                channel = %entry.id,
                error = %error,
                "2D solve failed; using reference trace fallback metrics"
            );
            (
                T::from_f64(shear_pa).expect("analytical constant conversion"),
                T::from_f64(shear_pa).expect("analytical constant conversion"),
                reference_trace.pressure_drop_pa,
                T::zero(),
                reference_trace.pressure_drop_pa,
                reference_trace.resistance_pa_s_per_m3,
                reference_trace.flow_rate_m3_s,
                T::zero(),
                T::zero(),
                None,
                SolveResult::new(0, T::zero(), false),
            )
        }
    };

    let t_s = area * entry.length_m / flow_rate_m3_s.abs().max(MIN_FLOW_RATE_M3_S);
    let hi = HemolysisModel::giersiepen_millifluidic()
        .damage_index(shear_pa, t_s)
        .unwrap_or(0.0);
    let hi_t = T::from_f64(hi).expect("analytical constant conversion");

    Ok(Channel2dResult {
        channel_id: entry.id.clone(),
        therapy_zone: entry.therapy_zone,
        is_venturi_throat: entry.is_venturi_throat,
        solve_result,
        projection: entry.projection.clone(),
        wall_shear_pa: T::from_f64(shear_pa).expect("analytical constant conversion"),
        field_wall_shear_max_pa,
        field_wall_shear_mean_pa,
        field_inlet_pressure_pa,
        field_outlet_pressure_pa,
        field_pressure_drop_pa,
        field_effective_resistance_pa_s_per_m3,
        field_outlet_flow_m3_s,
        field_outlet_flow_error_m3_s,
        field_outlet_flow_error_pct,
        transit_time_s: T::from_f64(t_s).expect("analytical constant conversion"),
        field_separation_efficiency_pct,
        hemolysis_index: hi_t,
        reference_trace: reference_trace.clone(),
    })
}
