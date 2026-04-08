//! Private build-pipeline helpers for BlueprintTopologyFactory.
use crate::config::{ChannelTypeConfig, GeometryConfig};
use crate::domain::model::{EdgeId, NetworkBlueprint};
use crate::domain::therapy_metadata::TherapyZone;
use crate::geometry::generator::{
    create_parallel_geometry_from_spec, create_primitive_selective_tree_geometry_from_spec,
    create_series_geometry_from_spec, GeometryGeneratorBuilder,
};
use crate::geometry::metadata::ChannelVisualRole;
use crate::topology::model::{
    BlueprintTopologySpec, TopologyLineageMetadata, VenturiPlacementSpec,
};
use super::BlueprintTopologyFactory;

impl BlueprintTopologyFactory {
    pub(super) fn build_series_path(
        spec: &BlueprintTopologySpec,
        lineage: TopologyLineageMetadata,
    ) -> crate::domain::model::NetworkBlueprint {
        let mut blueprint = create_series_geometry_from_spec(spec);
        blueprint.lineage = Some(lineage);
        blueprint.render_hints = Some(Self::render_hints_for_spec(spec));
        blueprint
    }

    pub(super) fn build_parallel_path(
        spec: &BlueprintTopologySpec,
        lineage: TopologyLineageMetadata,
    ) -> crate::domain::model::NetworkBlueprint {
        let mut blueprint = create_parallel_geometry_from_spec(spec);
        blueprint.lineage = Some(lineage);
        blueprint.render_hints = Some(Self::render_hints_for_spec(spec));
        blueprint
    }

    /// Construct a split-tree blueprint via the canonical GeometryGeneratorBuilder.
    pub(super) fn build_split_tree(
        spec: &BlueprintTopologySpec,
        lineage: TopologyLineageMetadata,
    ) -> Result<NetworkBlueprint, String> {
        if let Some(mut blueprint) = create_primitive_selective_tree_geometry_from_spec(spec)? {
            blueprint.lineage = Some(lineage);
            return Ok(blueprint);
        }

        let splits = Self::spec_to_split_types(spec)?;

        let representative_width_m = Self::representative_channel_width_m(spec);
        let representative_height_m = Self::representative_channel_height_m(spec);

        let geo_config = GeometryConfig {
            channel_width: representative_width_m * 1e3,   // m → mm
            channel_height: representative_height_m * 1e3, // m → mm
            ..GeometryConfig::default()
        };

        let mut blueprint = GeometryGeneratorBuilder::new(
            spec.box_dims_mm,
            &splits,
            &geo_config,
            &ChannelTypeConfig::AllStraight,
        )
        .with_topology_spec(spec.clone())
        .with_lineage(lineage)
        .with_render_hints(Self::render_hints_for_spec(spec))
        .build();

        blueprint.name.clone_from(&spec.design_name);

        // Reconcile auto-generated channel IDs with the spec's semantic names
        Self::reconcile_channel_ids(&mut blueprint, spec);

        Ok(blueprint)
    }

    /// Reconcile auto-generated channel IDs (`ch_0`, `ch_1`, ...) with the
    /// topology spec's semantic channel names.
    ///
    /// For series-path specs, channels are renamed 1:1 in positional order.
    /// For split-tree specs, only the trunk and outlet channels can be
    /// reconciled positionally; branch channels keep their auto-generated IDs.
    pub(super) fn reconcile_channel_ids(blueprint: &mut NetworkBlueprint, spec: &BlueprintTopologySpec) {
        if !spec.series_channels.is_empty() {
            // Series path: channels map 1:1 in order
            let named_ids: Vec<String> = spec
                .series_channels
                .iter()
                .map(|ch| ch.channel_id.clone())
                .collect();
            for (i, name) in named_ids.into_iter().enumerate() {
                if let Some(channel) = blueprint.channels.get_mut(i) {
                    channel.id = EdgeId(name);
                }
            }
        } else if !spec.parallel_channels.is_empty() {
            let named_ids: Vec<String> = spec
                .parallel_channels
                .iter()
                .map(|ch| ch.channel_id.clone())
                .collect();
            for (i, name) in named_ids.into_iter().enumerate() {
                if let Some(channel) = blueprint.channels.get_mut(i) {
                    channel.id = EdgeId(name);
                }
            }
        }
    }

    pub(super) fn has_materialized_venturi_geometry(
        blueprint: &NetworkBlueprint,
        spec: &BlueprintTopologySpec,
    ) -> bool {
        !spec.venturi_placements.is_empty()
            && spec.venturi_placements.iter().all(|placement| {
                blueprint.channels.iter().any(|channel| {
                    (channel.id.as_str() == placement.target_channel_id
                        || channel.id.as_str().starts_with(&placement.target_channel_id))
                        && channel.venturi_geometry.is_some()
                })
            })
    }

    pub(super) fn resolve_materialized_venturi_targets(
        blueprint: &NetworkBlueprint,
        spec: &BlueprintTopologySpec,
    ) -> BlueprintTopologySpec {
        let expanded = Self::selective_materialized_venturi_placements(blueprint, spec);
        if expanded.is_empty() {
            return spec.clone();
        }

        let mut resolved = spec.clone();
        resolved.venturi_placements = expanded;
        resolved
    }

    pub(super) fn selective_materialized_venturi_placements(
        blueprint: &NetworkBlueprint,
        spec: &BlueprintTopologySpec,
    ) -> Vec<VenturiPlacementSpec> {
        if spec.split_stages.is_empty() || spec.venturi_placements.is_empty() {
            return Vec::new();
        }

        let preferred = Self::leading_merge_side_treatment_channels(blueprint, true);
        let fallback = Self::leading_merge_side_treatment_channels(blueprint, false);
        let actual_ids = if fallback.len() > preferred.len() || preferred.is_empty() {
            fallback
        } else {
            preferred
        };

        if actual_ids.is_empty() {
            return Vec::new();
        }

        let template = spec
            .venturi_placements
            .iter()
            .max_by_key(|placement| placement.serial_throat_count)
            .or_else(|| spec.venturi_placements.first())
            .expect("non-empty venturi placements required for selective expansion");

        actual_ids
            .into_iter()
            .enumerate()
            .map(|(idx, materialized_id)| VenturiPlacementSpec {
                placement_id: format!("venturi_{idx}"),
                target_channel_id: materialized_id,
                serial_throat_count: template.serial_throat_count,
                throat_geometry: template.throat_geometry.clone(),
                placement_mode: template.placement_mode,
            })
            .collect()
    }

    pub(super) fn leading_merge_side_treatment_channels(
        blueprint: &NetworkBlueprint,
        require_venturi: bool,
    ) -> Vec<String> {
        let mut channels = blueprint
            .channels
            .iter()
            .filter(|channel| {
                channel.therapy_zone == Some(TherapyZone::CancerTarget)
                    && channel.visual_role != Some(ChannelVisualRole::Trunk)
                    && (!require_venturi || channel.venturi_geometry.is_some())
            })
            .filter_map(|channel| {
                let max_x = channel
                    .path
                    .iter()
                    .map(|(x, _)| *x)
                    .fold(f64::NEG_INFINITY, f64::max);
                (max_x > blueprint.box_dims.0 * 0.5 + 1.0e-6).then(|| {
                    let centroid_x =
                        channel.path.iter().map(|(x, _)| *x).sum::<f64>() / channel.path.len() as f64;
                    let centroid_y =
                        channel.path.iter().map(|(_, y)| *y).sum::<f64>() / channel.path.len() as f64;
                    (channel.id.0.clone(), centroid_x, centroid_y)
                })
            })
            .collect::<Vec<_>>();

        channels.sort_by(|left, right| {
            left.1
                .total_cmp(&right.1)
                .then(left.2.total_cmp(&right.2))
                .then(left.0.cmp(&right.0))
        });
        let Some((_, lead_x, _)) = channels.first() else {
            return Vec::new();
        };
        let lead_x = *lead_x;
        channels
            .into_iter()
            .filter(|(_, centroid_x, _)| (*centroid_x - lead_x).abs() <= 1.0e-6)
            .map(|(channel_id, _, _)| channel_id)
            .collect()
    }

}
