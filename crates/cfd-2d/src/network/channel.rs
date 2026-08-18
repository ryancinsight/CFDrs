use aequitas::systems::si::quantities::{DynamicViscosity, Length, MassDensity, Time, Velocity};
use cfd_core::error::{Error, Result as CfdResult};
use cfd_core::physics::hemolysis::HemolysisModel;
use cfd_core::CfdScalar;
use cfd_schematics::domain::model::CrossSectionSpec;
use eunomia::{FloatElement, NumericElement};

use crate::scalar;
use crate::solvers::cell_tracking::physics::{CellTrackerConfig, VelocityFieldInterpolator};
use crate::solvers::cell_tracking::population::{CellPopulation, TrackedCell};
use crate::solvers::cell_tracking::tracker::CellTracker;

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
    T: CfdScalar + Copy + FloatElement + eunomia::RealField + std::fmt::Debug,
{
    let zero: T = scalar::zero();
    let tol_t = <T as FloatElement>::from_f64(tolerance);
    entry.flow_rate_m3_s = flow_rate_m3_s;
    entry.reference_trace = reference_trace.clone();

    let area = entry.cross_section_area_m2;
    let area_safe = area.max(MIN_CROSS_SECTION_AREA_M2);
    let u_mean = flow_rate_m3_s / area_safe;
    let flow_direction = if flow_rate_m3_s.is_sign_negative() {
        -1.0
    } else {
        1.0
    };
    // The field solver owns a left velocity-inlet/right pressure-outlet
    // geometry. Solve the magnitude in that orientation and restore the
    // directed sign only on the extracted flow metric.
    let u_inlet = <T as FloatElement>::from_f64(u_mean.abs());

    entry.solver.config.tolerance = tol_t;

    let shear_rate = entry.cross_section.wall_shear_rate(u_mean);
    let shear_pa = entry.viscosity_pa_s * shear_rate.abs();
    // Network channels use the normalized parabolic profile that defines the
    // hydraulic reference coupling. Geometry-specific solvers retain the
    // uniform profile through `solve` unless they select this mode explicitly.
    let solve_result = entry.solver.solve_parabolic_inlet(u_inlet)?;

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
    ) = {
        if solve_result.converged {
            tracing::debug!(
                "  [2D] Channel {:<20} | u_in: {:.3e} | inner iter: {} | inner res: {:.3e} | Conv: {}",
                entry.id,
                <T as NumericElement>::to_f64(u_inlet),
                solve_result.iterations,
                <T as NumericElement>::to_f64(solve_result.residual),
                solve_result.converged
            );
        } else {
            tracing::warn!(
                channel = %entry.id,
                residual = <T as NumericElement>::to_f64(solve_result.residual),
                "2D solve did not converge; preserving extracted field metrics"
            );
        }

        let (field_wall_shear_max_pa, field_wall_shear_mean_pa) =
            extract_field_wall_shear(&entry.solver);
        let (field_inlet_pressure_pa, field_outlet_pressure_pa) =
            extract_field_inlet_outlet_pressure(&entry.solver);
        let field_pressure_drop_pa = field_inlet_pressure_pa - field_outlet_pressure_pa;
        let field_outlet_flow_m3_s = <T as FloatElement>::from_f64(flow_direction)
            * extract_field_outlet_flow_rate(&entry.solver, <T as FloatElement>::from_f64(area));
        let field_outlet_flow_error_m3_s = field_outlet_flow_m3_s - reference_trace.flow_rate_m3_s;
        let min_flow: T = <T as FloatElement>::from_f64(MIN_FLOW_RATE_M3_S);
        let hundred: T = <T as FloatElement>::from_f64(100.0);
        let field_outlet_flow_error_pct =
            if <T as NumericElement>::abs(reference_trace.flow_rate_m3_s) > min_flow {
                <T as NumericElement>::abs(
                    field_outlet_flow_error_m3_s / reference_trace.flow_rate_m3_s,
                ) * hundred
            } else {
                zero
            };

        let flow_abs = <T as NumericElement>::abs(<T as FloatElement>::from_f64(flow_rate_m3_s));
        let mut field_effective_resistance_pa_s_per_m3 = reference_trace.resistance_pa_s_per_m3;
        if flow_abs > min_flow {
            field_effective_resistance_pa_s_per_m3 =
                <T as NumericElement>::abs(field_pressure_drop_pa) / flow_abs;
        }

        let mut field_separation_efficiency_pct = None;
        if separation_tracking_enabled
            && matches!(entry.cross_section, CrossSectionSpec::Rectangular { .. })
        {
            let cnf = CellTrackerConfig {
                viscosity: DynamicViscosity::from_base(entry.viscosity_pa_s),
                fluid_density: MassDensity::from_base(1060.0),
                hydraulic_diameter: entry.cross_section.hydraulic_diameter(),
                u_max: Velocity::from_base(u_mean.abs()),
                ..Default::default()
            };
            let tracker = CellTracker::new(&entry.solver, cnf);

            let n_per_pop = 50usize;
            let mut cells = Vec::with_capacity(n_per_pop * 2);
            let (x_min_q, _, y_min_q, y_max_q) = entry.solver.bounds();
            let x_min = x_min_q.into_base();
            let y_min = y_min_q.into_base();
            let y_max = y_max_q.into_base();
            let h_bounds = y_max - y_min;

            for i in 0..n_per_pop {
                let y =
                    y_min + h_bounds * 0.05 + h_bounds * 0.9 * (i as f64 / (n_per_pop - 1) as f64);
                cells.push(TrackedCell {
                    population: CellPopulation::CTC,
                    x: Length::from_base(x_min + 1e-6),
                    y: Length::from_base(y),
                    vx: Velocity::from_base(u_mean.abs()),
                    vy: Velocity::from_base(0.0),
                    id: i,
                });
                cells.push(TrackedCell {
                    population: CellPopulation::RBC,
                    x: Length::from_base(x_min + 1e-6),
                    y: Length::from_base(y),
                    vx: Velocity::from_base(u_mean.abs()),
                    vy: Velocity::from_base(0.0),
                    id: n_per_pop + i,
                });
            }

            let trajectories = tracker.trace_cells(&cells, Time::from_base(1e-5), 5_000_000);
            let mut sum_y_ctc = 0.0;
            let mut count_ctc = 0usize;
            let mut sum_y_rbc = 0.0;
            let mut count_rbc = 0usize;

            for traj in &trajectories {
                if let Some(pos) = traj.positions.last() {
                    let y_tilde = (pos.y.into_base() - y_min) / h_bounds.max(1e-12);
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
            field_separation_efficiency_pct = Some(<T as FloatElement>::from_f64(diff_pct));
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
    };

    let t_s = area * entry.length_m / flow_rate_m3_s.abs().max(MIN_FLOW_RATE_M3_S);
    let hi = HemolysisModel::giersiepen_millifluidic()
        .damage_index(shear_pa, t_s)
        .map_err(|error| {
            Error::InvalidConfiguration(format!(
                "2D channel hemolysis inputs invalid (shear stress {shear_pa}, exposure time {t_s}): {error}"
            ))
        })?;
    let hi_t = <T as FloatElement>::from_f64(hi);

    Ok(Channel2dResult {
        channel_id: entry.id.clone(),
        therapy_zone: entry.therapy_zone,
        is_venturi_throat: entry.is_venturi_throat,
        solve_result,
        projection: entry.projection.clone(),
        wall_shear_pa: <T as FloatElement>::from_f64(shear_pa),
        field_wall_shear_max_pa,
        field_wall_shear_mean_pa,
        field_inlet_pressure_pa,
        field_outlet_pressure_pa,
        field_pressure_drop_pa,
        field_effective_resistance_pa_s_per_m3,
        field_outlet_flow_m3_s,
        field_outlet_flow_error_m3_s,
        field_outlet_flow_error_pct,
        transit_time_s: <T as FloatElement>::from_f64(t_s),
        field_separation_efficiency_pct,
        hemolysis_index: hi_t,
        reference_trace: reference_trace.clone(),
    })
}
