//! Shared Milestone 12 multi-fidelity venturi validation.
//!
//! Delegates the actual 1D / 2D / 3D pressure-drop computations to
//! [`cfd_validation::numerical::venturi_cross_fidelity`] and converts the
//! results into the report-facing [`ValidationRow`] type.

use std::path::Path;

use cfd_validation::numerical::venturi_cross_fidelity::{
    validate_venturi, VenturiCrossFidelityResult, VenturiValidationInput,
};

use crate::constraints::BLOOD_DENSITY_KG_M3;
use crate::reporting::{Milestone12ReportDesign, ValidationRow};

fn total_loss_coefficient(dp_pa: f64, inlet_velocity_m_s: f64) -> f64 {
    let dynamic_pressure_pa = 0.5 * BLOOD_DENSITY_KG_M3 * inlet_velocity_m_s * inlet_velocity_m_s;
    if dynamic_pressure_pa > 1.0e-18 {
        dp_pa / dynamic_pressure_pa
    } else {
        0.0
    }
}

fn sanitize_report_scalar(value: f64) -> f64 {
    if value.is_finite() {
        value
    } else {
        0.0
    }
}

fn validation_row_from_result(
    track: &str,
    id: &str,
    topology: &str,
    score: f64,
    sigma_1d: f64,
    result: &VenturiCrossFidelityResult,
) -> ValidationRow {
    let inlet_velocity_m_s = result.input.inlet_velocity_1d();

    ValidationRow {
        track: track.to_string(),
        id: id.to_string(),
        topology: topology.to_string(),
        k_loss_1d: sanitize_report_scalar(total_loss_coefficient(
            result.dp_1d_pa,
            inlet_velocity_m_s,
        )),
        dp_1d_bernoulli_pa: sanitize_report_scalar(result.dp_1d_pa),
        k_loss_2d: sanitize_report_scalar(total_loss_coefficient(
            result.dp_2d_pa,
            inlet_velocity_m_s,
        )),
        dp_2d_fvm_pa: sanitize_report_scalar(result.dp_2d_pa),
        k_loss_3d: sanitize_report_scalar(total_loss_coefficient(
            result.dp_3d_pa,
            inlet_velocity_m_s,
        )),
        dp_3d_fem_pa: sanitize_report_scalar(result.dp_3d_pa),
        agreement_1d_2d_pct: sanitize_report_scalar(result.diff_1d_2d_pct),
        agreement_2d_3d_pct: sanitize_report_scalar(result.diff_2d_3d_pct),
        mass_error_3d_pct: sanitize_report_scalar(result.mass_error_3d_pct),
        sigma_1d: sanitize_report_scalar(sigma_1d),
        sigma_2d: sanitize_report_scalar(result.sigma_2d),
        score: sanitize_report_scalar(score),
        two_d_converged: result.two_d_converged,
        high_re_laminar_mismatch: result.high_re_laminar_mismatch,
    }
}

fn validate_venturi_candidate(
    track: &str,
    design: &Milestone12ReportDesign,
    out_dir: &Path,
) -> Result<ValidationRow, Box<dyn std::error::Error>> {
    let candidate = &design.candidate;
    let metrics = &design.metrics;
    let topology = candidate
        .topology_spec()
        .map_err(|error| error.to_string())?;
    let venturi = topology
        .venturi_placements
        .first()
        .ok_or_else(|| format!("validation requires venturi topology: {}", candidate.id))?;

    let parallel_paths = (metrics.active_venturi_throat_count
        / metrics.serial_venturi_stages_per_path.max(1))
    .max(1) as f64;
    let q = candidate.operating_point.flow_rate_m3_s * metrics.therapy_channel_fraction
        / parallel_paths;

    let input = VenturiValidationInput {
        label: candidate.id.clone(),
        inlet_diameter_m: venturi.throat_geometry.inlet_width_m,
        throat_diameter_m: venturi.throat_geometry.throat_width_m,
        throat_length_m: venturi.throat_geometry.throat_length_m,
        flow_rate_m3_s: q,
        inlet_gauge_pa: candidate.operating_point.inlet_gauge_pa,
    };

    let result = validate_venturi(&input);

    // Write per-case artifact.
    std::fs::write(
        out_dir.join(format!("{}_cross_fidelity.json", candidate.id)),
        serde_json::to_string_pretty(&result)?,
    )?;

    let row = validation_row_from_result(
        track,
        &candidate.id,
        &topology.short_code(),
        design.score,
        metrics.cavitation_number,
        &result,
    );
    std::fs::write(
        out_dir.join(format!("{}_validation.json", candidate.id)),
        serde_json::to_string_pretty(&row)?,
    )?;
    Ok(row)
}

/// Run the canonical Milestone 12 venturi validation on the selected Option 2
/// and GA-ranked designs, writing per-design artifacts and the aggregate row set.
///
/// # Errors
/// Returns an error if the underlying 2D/3D solvers fail or any artifact cannot
/// be written.
pub fn run_milestone12_validation(
    out_dir: &Path,
    option2: &Milestone12ReportDesign,
    ga: &Milestone12ReportDesign,
) -> Result<Vec<ValidationRow>, Box<dyn std::error::Error>> {
    std::fs::create_dir_all(out_dir)?;
    let rows = vec![
        validate_venturi_candidate("Option 2 Combined Selective Venturi", option2, out_dir)?,
        validate_venturi_candidate("GA HydroSDT Venturi", ga, out_dir)?,
    ];
    std::fs::write(
        out_dir.join("milestone12_validation_rows.json"),
        serde_json::to_string_pretty(&rows)?,
    )?;
    Ok(rows)
}

#[cfg(test)]
mod tests {
    use super::validation_row_from_result;
    use cfd_validation::numerical::venturi_cross_fidelity::{
        Fidelity2DResult, Fidelity3DResult, FidelityBreakdown1D, VenturiCrossFidelityResult,
        VenturiValidationInput,
    };

    #[test]
    fn validation_row_carries_dimensionless_loss_coefficients() {
        let input = VenturiValidationInput {
            label: "demo".to_string(),
            inlet_diameter_m: 2.0e-3,
            throat_diameter_m: 1.0e-3,
            throat_length_m: 1.0e-3,
            flow_rate_m3_s: 1.2e-6,
            inlet_gauge_pa: 30_000.0,
        };
        let result = VenturiCrossFidelityResult {
            label: "demo".to_string(),
            input: input.clone(),
            dp_1d_pa: 1200.0,
            breakdown_1d: FidelityBreakdown1D {
                dp_contraction_pa: 600.0,
                dp_friction_pa: 450.0,
                dp_expansion_loss_pa: 300.0,
                dp_recovery_pa: 150.0,
                dp_total_pa: 1200.0,
                discharge_coefficient: 0.96,
                friction_factor: 0.032,
                re_throat: 1800.0,
            },
            dp_2d_pa: 1000.0,
            result_2d: Fidelity2DResult {
                u_inlet: 0.5,
                u_throat: 2.0,
                dp_throat_pa: 900.0,
                dp_recovery_pa: 100.0,
                sigma: 0.6,
                converged: true,
            },
            dp_3d_pa: 900.0,
            result_3d: Fidelity3DResult {
                u_inlet: 0.5,
                u_throat: 2.0,
                dp_throat_pa: 850.0,
                dp_recovery_pa: 50.0,
                mass_error: 0.01,
                resolution: (48, 24),
            },
            diff_1d_2d_pct: 16.6666666667,
            diff_2d_3d_pct: 10.0,
            mass_error_3d_pct: 1.0,
            two_d_converged: true,
            high_re_laminar_mismatch: false,
            re_throat: 1800.0,
            sigma_2d: 0.6,
        };

        let row = validation_row_from_result("track", "id", "topo", 0.7, 0.5, &result);
        let dynamic_pressure_pa =
            0.5 * crate::constraints::BLOOD_DENSITY_KG_M3 * input.inlet_velocity_1d().powi(2);

        assert!((row.k_loss_1d - 1200.0 / dynamic_pressure_pa).abs() < 1.0e-12);
        assert!((row.k_loss_2d - 1000.0 / dynamic_pressure_pa).abs() < 1.0e-12);
        assert!((row.k_loss_3d - 900.0 / dynamic_pressure_pa).abs() < 1.0e-12);
        assert_eq!(row.sigma_1d, 0.5);
        assert_eq!(row.sigma_2d, 0.6);
    }
}
