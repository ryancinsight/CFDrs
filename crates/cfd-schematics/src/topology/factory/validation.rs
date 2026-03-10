//! Validation logic for [`BlueprintTopologySpec`] and venturi throat geometry.

use crate::topology::model::{BlueprintTopologySpec, ChannelRouteSpec, ThroatGeometrySpec};

/// Validate an entire topology spec before building.
///
/// Called by [`BlueprintTopologyFactory::build`] and
/// [`NetworkBlueprint::validate`].
pub fn validate_spec(spec: &BlueprintTopologySpec) -> Result<(), String> {
    if spec.design_name.is_empty() {
        return Err("BlueprintTopologySpec.design_name is empty".into());
    }
    let (bw, bh) = spec.box_dims_mm;
    if bw <= 0.0 || bh <= 0.0 {
        return Err(format!(
            "Box dimensions must be positive: ({bw}, {bh})"
        ));
    }
    if spec.inlet_width_m <= 0.0 {
        return Err(format!(
            "inlet_width_m must be positive: {}",
            spec.inlet_width_m
        ));
    }
    if spec.outlet_width_m <= 0.0 {
        return Err(format!(
            "outlet_width_m must be positive: {}",
            spec.outlet_width_m
        ));
    }
    if spec.trunk_length_m <= 0.0 {
        return Err(format!(
            "trunk_length_m must be positive: {}",
            spec.trunk_length_m
        ));
    }
    if spec.outlet_tail_length_m <= 0.0 {
        return Err(format!(
            "outlet_tail_length_m must be positive: {}",
            spec.outlet_tail_length_m
        ));
    }

    // Validate series channels
    for ch in &spec.series_channels {
        validate_route(&format!("series:{}", ch.channel_id), &ch.route)?;
    }

    // Validate parallel channels
    for ch in &spec.parallel_channels {
        validate_route(&format!("parallel:{}", ch.channel_id), &ch.route)?;
    }

    // Validate split stages
    for stage in &spec.split_stages {
        if stage.branches.is_empty() {
            return Err(format!(
                "Split stage '{}' has no branches",
                stage.stage_id
            ));
        }
        let expected = stage.split_kind.branch_count();
        if stage.branches.len() != expected {
            return Err(format!(
                "Split stage '{}' expects {} branches for {:?}, found {}",
                stage.stage_id,
                expected,
                stage.split_kind,
                stage.branches.len()
            ));
        }
        for branch in &stage.branches {
            validate_route(
                &format!("{}:{}", stage.stage_id, branch.label),
                &branch.route,
            )?;
        }
    }

    // Validate venturi placements
    for vp in &spec.venturi_placements {
        if vp.serial_throat_count == 0 {
            return Err(format!(
                "Venturi placement '{}' has zero serial_throat_count",
                vp.placement_id
            ));
        }
        validate_throat_geometry(&vp.throat_geometry, &vp.placement_id)?;
    }

    Ok(())
}

/// Validate a single channel route specification.
pub fn validate_route(label: &str, route: &ChannelRouteSpec) -> Result<(), String> {
    if route.length_m <= 0.0 {
        return Err(format!("{label}: length_m must be positive: {}", route.length_m));
    }
    if route.width_m <= 0.0 {
        return Err(format!("{label}: width_m must be positive: {}", route.width_m));
    }
    if route.height_m <= 0.0 {
        return Err(format!("{label}: height_m must be positive: {}", route.height_m));
    }
    if let Some(ref serp) = route.serpentine {
        if serp.segments == 0 {
            return Err(format!("{label}: serpentine segments must be > 0"));
        }
        if serp.bend_radius_m <= 0.0 {
            return Err(format!("{label}: serpentine bend_radius_m must be positive"));
        }
        if serp.segment_length_m <= 0.0 {
            return Err(format!("{label}: serpentine segment_length_m must be positive"));
        }
    }
    Ok(())
}

/// Validate venturi throat geometry constraints.
pub fn validate_throat_geometry(
    geometry: &ThroatGeometrySpec,
    placement_id: &str,
) -> Result<(), String> {
    if geometry.throat_width_m <= 0.0 {
        return Err(format!(
            "Venturi '{}': throat_width_m must be positive",
            placement_id
        ));
    }
    if geometry.throat_height_m <= 0.0 {
        return Err(format!(
            "Venturi '{}': throat_height_m must be positive",
            placement_id
        ));
    }
    if geometry.throat_length_m <= 0.0 {
        return Err(format!(
            "Venturi '{}': throat_length_m must be positive",
            placement_id
        ));
    }
    if geometry.inlet_width_m <= 0.0 {
        return Err(format!(
            "Venturi '{}': inlet_width_m must be positive",
            placement_id
        ));
    }
    if geometry.outlet_width_m <= 0.0 {
        return Err(format!(
            "Venturi '{}': outlet_width_m must be positive",
            placement_id
        ));
    }
    if geometry.convergent_half_angle_deg <= 0.0 || geometry.convergent_half_angle_deg >= 90.0 {
        return Err(format!(
            "Venturi '{}': convergent_half_angle_deg must be in (0, 90): {}",
            placement_id, geometry.convergent_half_angle_deg
        ));
    }
    if geometry.divergent_half_angle_deg <= 0.0 || geometry.divergent_half_angle_deg >= 90.0 {
        return Err(format!(
            "Venturi '{}': divergent_half_angle_deg must be in (0, 90): {}",
            placement_id, geometry.divergent_half_angle_deg
        ));
    }
    Ok(())
}
