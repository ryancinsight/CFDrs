//! Fluid volume summary for channel systems.

use serde::{Deserialize, Serialize};

use crate::domain::model::NetworkBlueprint;

/// Summary of the fluid volume within a channel system.
#[derive(Debug, Clone, Default, Serialize, Deserialize)]
pub struct FluidVolumeSummary {
    /// Total centerline length in the schematic [mm].
    pub total_channel_length_mm: f64,
    /// Total fluid volume in the system [mm^3].
    pub total_fluid_volume_mm3: f64,
    /// Total fluid volume in the system [uL].
    pub total_fluid_volume_ul: f64,
    /// Number of channels contributing to the volume.
    pub channel_count: usize,
    /// Human-readable legend label for report and plot overlays.
    pub display_label: String,
}

/// Per-channel fluid-volume summary derived from the authoritative blueprint geometry.
#[derive(Debug, Clone, Default, Serialize, Deserialize)]
pub struct ChannelFluidVolumeSummary {
    /// Blueprint channel identifier.
    pub channel_id: String,
    /// Upstream blueprint node identifier.
    pub from_node_id: String,
    /// Downstream blueprint node identifier.
    pub to_node_id: String,
    /// Total channel centerline length in the schematic [mm].
    pub centerline_length_mm: f64,
    /// True blueprint cross-sectional area [mm^2].
    pub cross_section_area_mm2: f64,
    /// Total fluid volume for this channel [mm^3].
    pub fluid_volume_mm3: f64,
    /// Total fluid volume for this channel [uL].
    pub fluid_volume_ul: f64,
}

impl NetworkBlueprint {
    #[must_use]
    pub fn channel_fluid_volume_summaries(&self) -> Vec<ChannelFluidVolumeSummary> {
        self.channels
            .iter()
            .map(|channel| {
                let centerline_length_mm = if channel.path.len() >= 2 {
                    channel
                        .path
                        .windows(2)
                        .map(|segment| {
                            let dx = segment[1].0 - segment[0].0;
                            let dy = segment[1].1 - segment[0].1;
                            dx.hypot(dy)
                        })
                        .sum()
                } else {
                    channel.length_m * 1000.0
                };
                let cross_section_area_mm2 = channel.cross_section.area() * 1.0e6;
                let fluid_volume_mm3 = centerline_length_mm * cross_section_area_mm2;

                ChannelFluidVolumeSummary {
                    channel_id: channel.id.as_str().to_string(),
                    from_node_id: channel.from.to_string(),
                    to_node_id: channel.to.to_string(),
                    centerline_length_mm,
                    cross_section_area_mm2,
                    fluid_volume_mm3,
                    fluid_volume_ul: fluid_volume_mm3,
                }
            })
            .collect()
    }

    #[must_use]
    pub fn fluid_volume_summary(&self) -> FluidVolumeSummary {
        let channel_summaries = self.channel_fluid_volume_summaries();
        let total_channel_length_mm = channel_summaries
            .iter()
            .map(|channel| channel.centerline_length_mm)
            .sum();
        let total_fluid_volume_mm3 = channel_summaries
            .iter()
            .map(|channel| channel.fluid_volume_mm3)
            .sum();

        let total_fluid_volume_ul = total_fluid_volume_mm3;

        FluidVolumeSummary {
            total_channel_length_mm,
            total_fluid_volume_mm3,
            total_fluid_volume_ul,
            channel_count: channel_summaries.len(),
            display_label: format!(
                "Volume: {:.3} uL over {:.2} mm",
                total_fluid_volume_ul, total_channel_length_mm
            ),
        }
    }
}
