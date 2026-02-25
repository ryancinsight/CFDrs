//! Network topology analysis for `NetworkBlueprint`.
//!
//! Classifies a blueprint into a topology class and provides graph-traversal
//! helpers used by the mesh pipeline to lay out channel segments.

use std::collections::HashMap;

use cfd_schematics::{ChannelSpec, NetworkBlueprint, NodeKind};
use cfd_schematics::domain::therapy_metadata::{TherapyZone, TherapyZoneMetadata};

/// Topology classification of a `NetworkBlueprint`.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum TopologyClass {
    /// A single linear chain of `n_segments` channels (inlet → … → outlet).
    LinearChain {
        /// Number of channel segments in the chain.
        n_segments: usize,
    },
    /// Symmetric bifurcation: inlet → parent → (d1, d2) → parent_out → outlet.
    Bifurcation,
    /// Symmetric trifurcation: inlet → parent → (d1, d2, d3) → outlets.
    Trifurcation,
    /// Venturi chain: linear chain where at least one channel is tagged
    /// `TherapyZone::CancerTarget`.
    VenturiChain,
    /// Any topology not covered by the above.
    Complex,
}

/// Graph-analysis overlay on a `NetworkBlueprint`.
pub struct NetworkTopology<'bp> {
    bp: &'bp NetworkBlueprint,
    /// `node_id → degree` (total number of adjacent channels).
    degrees: HashMap<String, usize>,
}

impl<'bp> NetworkTopology<'bp> {
    /// Build the topology analysis for `bp`.
    pub fn new(bp: &'bp NetworkBlueprint) -> Self {
        let mut degrees: HashMap<String, usize> = HashMap::new();
        for node in &bp.nodes {
            degrees.insert(node.id.to_string(), 0);
        }
        for ch in &bp.channels {
            *degrees.entry(ch.from.to_string()).or_insert(0) += 1;
            *degrees.entry(ch.to.to_string()).or_insert(0) += 1;
        }
        Self { bp, degrees }
    }

    /// Total degree (number of adjacent channels) of a node.
    pub fn degree(&self, node_id: &str) -> usize {
        self.degrees.get(node_id).copied().unwrap_or(0)
    }

    /// Channels that originate **from** `node_id`.
    pub fn outgoing_channels(&self, node_id: &str) -> Vec<&ChannelSpec> {
        self.bp
            .channels
            .iter()
            .filter(|c| c.from.as_str() == node_id)
            .collect()
    }

    /// Channels that terminate **at** `node_id`.
    pub fn incoming_channels(&self, node_id: &str) -> Vec<&ChannelSpec> {
        self.bp
            .channels
            .iter()
            .filter(|c| c.to.as_str() == node_id)
            .collect()
    }

    /// ID of the unique inlet node, or `None` if the blueprint has != 1 inlet.
    pub fn inlet_node_id(&self) -> Option<&str> {
        let inlets: Vec<_> = self
            .bp
            .nodes
            .iter()
            .filter(|n| matches!(n.kind, NodeKind::Inlet))
            .collect();
        if inlets.len() == 1 {
            Some(inlets[0].id.as_str())
        } else {
            None
        }
    }

    /// IDs of all outlet nodes.
    pub fn outlet_node_ids(&self) -> Vec<&str> {
        self.bp
            .nodes
            .iter()
            .filter(|n| matches!(n.kind, NodeKind::Outlet))
            .map(|n| n.id.as_str())
            .collect()
    }

    /// Classify the network topology.
    pub fn classify(&self) -> TopologyClass {
        let n_in = self.bp.inlet_count();
        let n_out = self.bp.outlet_count();
        let n_junc = self.bp.junction_count();
        let n_ch = self.bp.channels.len();

        // VenturiChain: linear chain where one channel has CancerTarget metadata.
        if n_in == 1 && n_out == 1 && n_junc == 2 && n_ch == 3 {
            let has_cancer = self.bp.channels.iter().any(|c| {
                c.metadata
                    .as_ref()
                    .and_then(|m| m.get::<TherapyZoneMetadata>())
                    .map(|m| m.zone == TherapyZone::CancerTarget)
                    .unwrap_or(false)
            });
            if has_cancer && self.linear_path_channels().is_some() {
                return TopologyClass::VenturiChain;
            }
        }

        // Bifurcation: n_in=1, n_out=1, n_junc=2, n_ch=4
        if n_in == 1 && n_out == 1 && n_junc == 2 && n_ch == 4 {
            let junction_ids: Vec<_> = self
                .bp
                .nodes
                .iter()
                .filter(|n| matches!(n.kind, NodeKind::Junction))
                .map(|n| n.id.as_str())
                .collect();
            let both_degree3 = junction_ids
                .iter()
                .all(|id| self.degree(id) == 3);
            if both_degree3 {
                return TopologyClass::Bifurcation;
            }
        }

        // Trifurcation: n_in=1, n_out=3, n_junc=1, n_ch=4
        if n_in == 1 && n_out == 3 && n_junc == 1 && n_ch == 4 {
            return TopologyClass::Trifurcation;
        }

        // LinearChain: all junctions have degree 2, linear path exists
        if n_in == 1 && n_out == 1 {
            let all_junc_degree2 = self
                .bp
                .nodes
                .iter()
                .filter(|n| matches!(n.kind, NodeKind::Junction))
                .all(|n| self.degree(n.id.as_str()) == 2);
            if all_junc_degree2 {
                if let Some(path) = self.linear_path_channels() {
                    return TopologyClass::LinearChain {
                        n_segments: path.len(),
                    };
                }
            }
        }

        TopologyClass::Complex
    }

    /// Return channels in traversal order for a linear topology (inlet → outlet),
    /// or `None` if the topology is not linear.
    pub fn linear_path_channels(&self) -> Option<Vec<&ChannelSpec>> {
        let inlet_id = self.inlet_node_id()?;
        let outlet_ids = self.outlet_node_ids();
        if outlet_ids.len() != 1 {
            return None;
        }
        let outlet_id = outlet_ids[0];

        let mut path = Vec::new();
        let mut current = inlet_id;

        loop {
            let outgoing = self.outgoing_channels(current);
            if outgoing.len() != 1 {
                break;
            }
            let ch = outgoing[0];
            path.push(ch);
            current = ch.to.as_str();
            if current == outlet_id {
                return Some(path);
            }
        }
        None
    }

    /// Return `(parent_in, daughter_1, daughter_2, parent_out)` channels for a
    /// `Bifurcation` topology, or `None` otherwise.
    pub fn bifurcation_channels(
        &self,
    ) -> Option<(&ChannelSpec, &ChannelSpec, &ChannelSpec, &ChannelSpec)> {
        let inlet_id = self.inlet_node_id()?;

        // parent_in: inlet → diverging_junction
        let parent_in = self.outgoing_channels(inlet_id).into_iter().next()?;
        let div_junc = parent_in.to.as_str();

        // Two daughters from the diverging junction
        let mut daughters: Vec<&ChannelSpec> = self.outgoing_channels(div_junc);
        if daughters.len() != 2 {
            return None;
        }
        daughters.sort_by(|a, b| a.id.as_str().cmp(b.id.as_str()));
        let d1 = daughters[0];
        let d2 = daughters[1];

        // They should converge to the same junction
        if d1.to.as_str() != d2.to.as_str() {
            return None;
        }
        let conv_junc = d1.to.as_str();

        // parent_out: converging_junction → outlet
        let parent_out = self.outgoing_channels(conv_junc).into_iter().next()?;

        Some((parent_in, d1, d2, parent_out))
    }

    /// Return `(parent, daughter_1, daughter_2, daughter_3)` channels for a
    /// `Trifurcation` topology, or `None` otherwise.
    pub fn trifurcation_channels(
        &self,
    ) -> Option<(&ChannelSpec, &ChannelSpec, &ChannelSpec, &ChannelSpec)> {
        let inlet_id = self.inlet_node_id()?;
        let parent = self.outgoing_channels(inlet_id).into_iter().next()?;
        let junc = parent.to.as_str();
        let mut daughters: Vec<&ChannelSpec> = self.outgoing_channels(junc);
        if daughters.len() != 3 {
            return None;
        }
        daughters.sort_by(|a, b| a.id.as_str().cmp(b.id.as_str()));
        Some((parent, daughters[0], daughters[1], daughters[2]))
    }
}

#[cfg(test)]
mod tests {
    use cfd_schematics::interface::presets::{
        serpentine_chain, symmetric_bifurcation, symmetric_trifurcation, venturi_chain,
    };

    use super::*;

    #[test]
    fn venturi_classifies_as_venturi_chain() {
        let bp = venturi_chain("v", 0.030, 0.004, 0.002);
        let topo = NetworkTopology::new(&bp);
        assert_eq!(topo.classify(), TopologyClass::VenturiChain);
    }

    #[test]
    fn bifurcation_classifies_correctly() {
        let bp = symmetric_bifurcation("b", 0.010, 0.010, 0.004, 0.003);
        let topo = NetworkTopology::new(&bp);
        assert_eq!(topo.classify(), TopologyClass::Bifurcation);
    }

    #[test]
    fn trifurcation_classifies_correctly() {
        let bp = symmetric_trifurcation("t", 0.010, 0.008, 0.004, 0.003);
        let topo = NetworkTopology::new(&bp);
        assert_eq!(topo.classify(), TopologyClass::Trifurcation);
    }

    #[test]
    fn serpentine_3_classifies_as_linear_chain() {
        let bp = serpentine_chain("s", 3, 0.010, 0.004);
        let topo = NetworkTopology::new(&bp);
        assert_eq!(topo.classify(), TopologyClass::LinearChain { n_segments: 3 });
    }

    #[test]
    fn degree_checks() {
        let bp = symmetric_bifurcation("b", 0.010, 0.010, 0.004, 0.003);
        let topo = NetworkTopology::new(&bp);
        assert_eq!(topo.degree("inlet"), 1);
        assert_eq!(topo.degree("diverging_junction"), 3);
        assert_eq!(topo.degree("outlet"), 1);
    }

    #[test]
    fn linear_path_order() {
        let bp = serpentine_chain("s", 3, 0.010, 0.004);
        let topo = NetworkTopology::new(&bp);
        let path = topo.linear_path_channels().expect("should have linear path");
        assert_eq!(path.len(), 3);
        assert_eq!(path[0].id.as_str(), "segment_1");
        assert_eq!(path[2].id.as_str(), "segment_3");
    }

    #[test]
    fn outlet_node_ids_count() {
        let bp = symmetric_trifurcation("t", 0.010, 0.008, 0.004, 0.003);
        let topo = NetworkTopology::new(&bp);
        assert_eq!(topo.outlet_node_ids().len(), 3);
    }

    #[test]
    fn bifurcation_channels_returns_four() {
        let bp = symmetric_bifurcation("b", 0.010, 0.010, 0.004, 0.003);
        let topo = NetworkTopology::new(&bp);
        let chans = topo.bifurcation_channels();
        assert!(chans.is_some());
    }
}
