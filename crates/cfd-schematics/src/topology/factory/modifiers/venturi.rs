//! Post-processing pass: insert venturi throat geometry into existing channels.

use crate::domain::model::NetworkBlueprint;
use crate::geometry::metadata::VenturiGeometryMetadata;
use crate::topology::model::BlueprintTopologySpec;

/// Applies venturi placements to matching channels in the built blueprint.
///
/// For each `VenturiPlacementSpec`, locates the target channel by ID and
/// attaches `VenturiGeometryMetadata` to it.
pub fn apply_venturi_placements(
    blueprint: &mut NetworkBlueprint,
    spec: &BlueprintTopologySpec,
) -> Result<(), String> {
    for vp in &spec.venturi_placements {
        let target = blueprint
            .channels
            .iter_mut()
            .find(|ch| ch.id.as_str() == vp.target_channel_id);

        match target {
            Some(channel) => {
                let meta = VenturiGeometryMetadata {
                    throat_width_m: vp.throat_geometry.throat_width_m,
                    throat_height_m: vp.throat_geometry.throat_height_m,
                    throat_length_m: vp.throat_geometry.throat_length_m,
                    inlet_width_m: vp.throat_geometry.inlet_width_m,
                    outlet_width_m: vp.throat_geometry.outlet_width_m,
                    convergent_half_angle_deg: vp.throat_geometry.convergent_half_angle_deg,
                    divergent_half_angle_deg: vp.throat_geometry.divergent_half_angle_deg,
                    throat_position: 0.5,
                };
                channel.venturi_geometry = Some(meta);
            }
            None => {
                return Err(format!(
                    "Venturi placement '{}' targets channel '{}' which does not exist in the blueprint",
                    vp.placement_id, vp.target_channel_id
                ));
            }
        }
    }
    Ok(())
}
