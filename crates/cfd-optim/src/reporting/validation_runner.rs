//! Shared Milestone 12 multi-fidelity venturi validation.
//!
//! Delegates the actual 1D / 2D / 3D pressure-drop computations to
//! [`cfd_validation::numerical::venturi_cross_fidelity`] and converts the
//! results into the report-facing [`ValidationRow`] type.

use std::path::Path;

use cfd_validation::numerical::venturi_cross_fidelity::{
    validate_venturi, VenturiValidationInput,
};

use crate::reporting::{Milestone12ReportDesign, ValidationRow};

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

    let row = ValidationRow {
        track: track.to_string(),
        id: candidate.id.clone(),
        topology: topology.short_code(),
        dp_1d_bernoulli_pa: result.dp_1d_pa,
        dp_2d_fvm_pa: result.dp_2d_pa,
        dp_3d_fem_pa: result.dp_3d_pa,
        agreement_1d_2d_pct: result.diff_1d_2d_pct,
        agreement_2d_3d_pct: result.diff_2d_3d_pct,
        mass_error_3d_pct: result.mass_error_3d_pct,
        sigma_1d: metrics.cavitation_number,
        sigma_2d: result.sigma_2d,
        score: design.score,
        two_d_converged: result.two_d_converged,
        high_re_laminar_mismatch: result.high_re_laminar_mismatch,
    };
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
