use super::NetworkBlueprint;
use crate::domain::therapy_metadata::TherapyZone;
use crate::geometry::metadata::{ChannelVenturiSpec, MetadataContainer, VenturiGeometryMetadata};
use crate::topology::{
    TopologyOptimizationStage, TreatmentActuationMode, VenturiConfig, VenturiPlacementSpec,
};

impl NetworkBlueprint {
    pub fn add_venturi(&mut self, config: &VenturiConfig) -> Result<(), String> {
        let target_channel_ids = if config.target_channel_ids.is_empty() {
            self.treatment_channel_ids()
        } else {
            config.target_channel_ids.clone()
        };
        if target_channel_ids.is_empty() {
            return Err("venturi augmentation requires at least one target channel".to_string());
        }

        let mut placements = Vec::with_capacity(target_channel_ids.len());
        for (index, channel_id) in target_channel_ids.into_iter().enumerate() {
            let channel = self
                .channels
                .iter_mut()
                .find(|channel| channel.id.as_str() == channel_id)
                .ok_or_else(|| {
                    format!(
                        "venturi augmentation target '{}' does not exist in blueprint '{}'",
                        channel_id, self.name
                    )
                })?;

            let resolved_inlet_width_m = if config.throat_geometry.inlet_width_m > 0.0 {
                config.throat_geometry.inlet_width_m
            } else {
                channel.effective_width_m()
            };
            let resolved_outlet_width_m = if config.throat_geometry.outlet_width_m > 0.0 {
                config.throat_geometry.outlet_width_m
            } else {
                channel.effective_width_m()
            };
            let venturi_geometry = VenturiGeometryMetadata {
                throat_width_m: config.throat_geometry.throat_width_m,
                throat_height_m: config.throat_geometry.throat_height_m,
                throat_length_m: config.throat_geometry.throat_length_m,
                inlet_width_m: resolved_inlet_width_m,
                outlet_width_m: resolved_outlet_width_m,
                convergent_half_angle_deg: config.throat_geometry.convergent_half_angle_deg,
                divergent_half_angle_deg: config.throat_geometry.divergent_half_angle_deg,
                throat_position: 0.5,
            };
            channel.venturi_geometry = Some(venturi_geometry.clone());
            if channel.metadata.is_none() {
                channel.metadata = Some(MetadataContainer::new());
            }
            let metadata = channel
                .metadata
                .as_mut()
                .expect("channel metadata container must exist");
            metadata.insert(venturi_geometry);
            metadata.insert(ChannelVenturiSpec {
                n_throats: config.serial_throat_count,
                is_ctc_stream: channel
                    .therapy_zone
                    .is_some_and(|zone| zone == TherapyZone::CancerTarget),
                throat_width_m: config.throat_geometry.throat_width_m,
                height_m: config.throat_geometry.throat_height_m,
                inter_throat_spacing_m: config.throat_geometry.throat_length_m,
            });

            placements.push(VenturiPlacementSpec {
                placement_id: format!("venturi_{index}"),
                target_channel_id: channel_id,
                serial_throat_count: config.serial_throat_count,
                throat_geometry: config.throat_geometry.clone(),
                placement_mode: config.placement_mode,
            });
        }

        if let Some(topology) = &mut self.topology {
            topology.venturi_placements = placements;
            topology.treatment_mode = TreatmentActuationMode::VenturiCavitation;
        }
        if let Some(lineage) = &mut self.lineage {
            lineage.current_stage =
                TopologyOptimizationStage::AsymmetricSplitVenturiCavitationSelectivity;
        }
        Ok(())
    }
}
