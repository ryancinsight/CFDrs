//! Post-processing pass: insert venturi throat geometry into existing channels.

use crate::domain::model::NetworkBlueprint;
use crate::topology::model::{BlueprintTopologySpec, VenturiConfig};

/// Applies venturi placements to matching channels in the built blueprint.
///
/// For each `VenturiPlacementSpec`, locates the target channel by ID and
/// attaches `VenturiGeometryMetadata` to it.
pub fn apply_venturi_placements(
    blueprint: &mut NetworkBlueprint,
    spec: &BlueprintTopologySpec,
) -> Result<(), String> {
    for vp in &spec.venturi_placements {
        blueprint.add_venturi(&VenturiConfig {
            target_channel_ids: vec![vp.target_channel_id.clone()],
            serial_throat_count: vp.serial_throat_count,
            throat_geometry: vp.throat_geometry.clone(),
            placement_mode: vp.placement_mode,
        })?;
    }
    Ok(())
}
