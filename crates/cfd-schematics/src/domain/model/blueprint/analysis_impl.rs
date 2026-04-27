use super::{ChannelOverlapAnalysis, NetworkBlueprint};
use crate::topology::BlueprintTopologyFactory;

impl NetworkBlueprint {
    #[must_use]
    pub fn inlet_count(&self) -> usize {
        self.nodes
            .iter()
            .filter(|node| matches!(node.kind, super::super::NodeKind::Inlet))
            .count()
    }

    #[must_use]
    pub fn outlet_count(&self) -> usize {
        self.nodes
            .iter()
            .filter(|node| matches!(node.kind, super::super::NodeKind::Outlet))
            .count()
    }

    #[must_use]
    pub fn junction_count(&self) -> usize {
        self.nodes
            .iter()
            .filter(|node| matches!(node.kind, super::super::NodeKind::Junction))
            .count()
    }

    #[must_use]
    pub fn pipe_count(&self) -> usize {
        self.channels
            .iter()
            .filter(|channel| matches!(channel.kind, super::super::EdgeKind::Pipe))
            .count()
    }

    #[must_use]
    pub fn total_length_m(&self) -> f64 {
        self.channels.iter().map(|channel| channel.length_m).sum()
    }

    #[must_use]
    pub fn length_in_zone(&self, zone: crate::domain::therapy_metadata::TherapyZone) -> f64 {
        use crate::domain::therapy_metadata::TherapyZoneMetadata;
        self.channels
            .iter()
            .filter(|channel| {
                channel
                    .metadata
                    .as_ref()
                    .and_then(|metadata| metadata.get::<TherapyZoneMetadata>())
                    .is_some_and(|therapy_zone| therapy_zone.zone == zone)
            })
            .map(|channel| channel.length_m)
            .sum()
    }

    #[must_use]
    pub fn venturi_channels(&self) -> Vec<&super::super::ChannelSpec> {
        self.channels
            .iter()
            .filter(|channel| {
                channel.venturi_geometry.is_some()
                    || channel.metadata.as_ref().is_some_and(
                        crate::geometry::metadata::MetadataContainer::contains::<
                            crate::geometry::metadata::VenturiGeometryMetadata,
                        >,
                    )
            })
            .collect()
    }

    pub fn resolve_channel_overlaps(&mut self) -> crate::geometry::IntersectionResult {
        crate::geometry::insert_intersection_nodes(self)
    }

    #[must_use]
    pub fn unresolved_channel_overlap_count(&self) -> usize {
        crate::geometry::unresolved_intersection_count(self)
    }

    #[must_use]
    pub fn has_unresolved_channel_overlaps(&self) -> bool {
        self.unresolved_channel_overlap_count() > 0
    }

    #[must_use]
    pub fn channel_overlap_analysis(&self) -> ChannelOverlapAnalysis {
        let channels_with_paths: Vec<(usize, f64)> = self
            .channels
            .iter()
            .enumerate()
            .filter(|(_, channel)| channel.path.len() >= 2)
            .map(|(index, channel)| {
                let width_mm = match channel.cross_section {
                    super::super::CrossSectionSpec::Rectangular { width_m, .. } => width_m * 1e3,
                    super::super::CrossSectionSpec::Circular { diameter_m } => diameter_m * 1e3,
                };
                (index, width_mm)
            })
            .collect();

        let mut worst = ChannelOverlapAnalysis::default();
        for (left_index, (channel_a_index, width_a)) in channels_with_paths.iter().enumerate() {
            let path_a = &self.channels[*channel_a_index].path;
            for (channel_b_index, width_b) in &channels_with_paths[left_index + 1..] {
                let path_b = &self.channels[*channel_b_index].path;
                let min_dist = min_path_distance(path_a, path_b);
                let half_sum = (width_a + width_b) * 0.5;
                if min_dist < half_sum {
                    let narrower = width_a.min(*width_b).max(1e-9);
                    let overlap_frac = ((half_sum - min_dist) / narrower).clamp(0.0, 1.0);
                    if overlap_frac > worst.max_overlap_fraction {
                        let wider = width_a.max(*width_b).max(1e-9);
                        worst.max_overlap_fraction = overlap_frac;
                        worst.width_ratio_at_worst = wider / narrower;
                        worst.overlap_pair_count += 1;
                    }
                }
            }
        }
        worst
    }

    #[must_use]
    pub fn max_channel_overlap_fraction(&self) -> f64 {
        self.channel_overlap_analysis().max_overlap_fraction
    }

    pub fn validate(&self) -> Result<(), String> {
        if self.nodes.is_empty() {
            return Err("NetworkBlueprint has no nodes".to_string());
        }
        if self.channels.is_empty() {
            return Err("NetworkBlueprint has no channels".to_string());
        }
        let mut node_ids = std::collections::HashSet::with_capacity(self.nodes.len());
        node_ids.extend(self.nodes.iter().map(|node| node.id.as_str()));
        for channel in &self.channels {
            if !node_ids.contains(channel.from.as_str()) {
                return Err(format!(
                    "Channel '{}' references unknown from-node '{}'",
                    channel.id.as_str(),
                    channel.from.as_str()
                ));
            }
            if !node_ids.contains(channel.to.as_str()) {
                return Err(format!(
                    "Channel '{}' references unknown to-node '{}'",
                    channel.id.as_str(),
                    channel.to.as_str()
                ));
            }
        }
        let overlap_count = self.unresolved_channel_overlap_count();
        if overlap_count > 0 {
            return Err(format!(
                "NetworkBlueprint '{}' contains {overlap_count} unresolved interior channel crossing(s)",
                self.name
            ));
        }
        if let Some(topology) = &self.topology {
            BlueprintTopologyFactory::validate_spec(topology)?;
            if topology.is_selective_routing() && !self.is_geometry_authored() {
                return Err(format!(
                    "NetworkBlueprint '{}' carries selective split-tree topology '{}' but was not authored through create_geometry()",
                    self.name,
                    topology.stage_sequence_label()
                ));
            }
        }
        Ok(())
    }

    #[must_use]
    pub fn describe(&self) -> String {
        let mut inlets = 0usize;
        let mut outlets = 0usize;
        let mut junctions = 0usize;
        for node in &self.nodes {
            match node.kind {
                super::super::NodeKind::Inlet => inlets += 1,
                super::super::NodeKind::Outlet => outlets += 1,
                super::super::NodeKind::Junction => junctions += 1,
                super::super::NodeKind::Reservoir => {}
            }
        }
        let mut pipes = 0usize;
        let mut total_len = 0.0f64;
        for channel in &self.channels {
            if matches!(channel.kind, super::super::EdgeKind::Pipe) {
                pipes += 1;
            }
            total_len += channel.length_m;
        }
        format!(
            "NetworkBlueprint '{}'\n\
             ─────────────────────────────────────\n\
             Nodes   : {} total  ({} inlets, {} outlets, {} junctions)\n\
             Channels: {}  ({} pipes)\n\
             Length  : {:.6} m total\n",
            self.name,
            self.nodes.len(),
            inlets,
            outlets,
            junctions,
            self.channels.len(),
            pipes,
            total_len,
        )
    }
}

fn min_path_distance(path_a: &[(f64, f64)], path_b: &[(f64, f64)]) -> f64 {
    let mut min_distance = f64::INFINITY;
    for &(ax, ay) in path_a {
        for &(bx, by) in path_b {
            let distance = ((ax - bx).powi(2) + (ay - by).powi(2)).sqrt();
            if distance < min_distance {
                min_distance = distance;
            }
        }
    }
    min_distance
}
