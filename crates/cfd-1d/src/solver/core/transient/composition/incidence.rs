//! Node-to-edge incidence of a network in compressed-row form.

use crate::domain::network::Network;
use cfd_core::CfdScalar;
use cfd_core::physics::fluid::FluidTrait;
use petgraph::visit::EdgeRef;
use std::ops::Index;

/// One edge as seen from a node it touches.
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
pub(super) struct IncidentEdge {
    pub(super) edge_index: usize,
    pub(super) source: usize,
    pub(super) target: usize,
}

/// Edges incident on each node, stored contiguously.
///
/// Node `n`'s edges are `edges[offsets[n]..offsets[n + 1]]`, in graph edge
/// order; a self-loop appears once. The mixing sweeps visit every node once
/// per fixed-point iteration, so one contiguous array replaces a heap
/// allocation per node on that path.
#[derive(Debug)]
pub(super) struct NodeIncidence {
    offsets: Vec<usize>,
    edges: Vec<IncidentEdge>,
}

impl NodeIncidence {
    /// Build the incidence of every node of `network` by a counting sort
    /// over its edges.
    pub(super) fn from_network<T: CfdScalar + Copy, F: FluidTrait<T> + Clone>(
        network: &Network<T, F>,
    ) -> Self {
        let node_count = network.node_count();
        let mut offsets = vec![0usize; node_count + 1];
        for edge_ref in network.graph.edge_references() {
            let source = edge_ref.source().index();
            let target = edge_ref.target().index();
            offsets[source + 1] += 1;
            if source != target {
                offsets[target + 1] += 1;
            }
        }
        for node in 0..node_count {
            offsets[node + 1] += offsets[node];
        }

        let mut cursor = offsets[..node_count].to_vec();
        let mut edges = vec![IncidentEdge::default(); offsets[node_count]];
        for edge_ref in network.graph.edge_references() {
            let edge = IncidentEdge {
                edge_index: edge_ref.id().index(),
                source: edge_ref.source().index(),
                target: edge_ref.target().index(),
            };
            edges[cursor[edge.source]] = edge;
            cursor[edge.source] += 1;
            if edge.source != edge.target {
                edges[cursor[edge.target]] = edge;
                cursor[edge.target] += 1;
            }
        }

        Self { offsets, edges }
    }
}

impl Index<usize> for NodeIncidence {
    type Output = [IncidentEdge];

    fn index(&self, node: usize) -> &[IncidentEdge] {
        &self.edges[self.offsets[node]..self.offsets[node + 1]]
    }
}

#[cfg(test)]
mod tests {
    use super::{IncidentEdge, NodeIncidence};
    use crate::domain::network::{Network, NetworkBuilder};
    use cfd_core::physics::fluid::database::water_20c;

    fn incident(edge_index: usize, source: usize, target: usize) -> IncidentEdge {
        IncidentEdge {
            edge_index,
            source,
            target,
        }
    }

    #[test]
    fn rows_list_incident_edges_in_graph_edge_order() {
        let mut builder = NetworkBuilder::<f64>::new();
        let inlet = builder.add_inlet("inlet".to_string());
        let junction = builder.add_junction("junction".to_string());
        let outlet = builder.add_outlet("outlet".to_string());
        builder.connect_with_pipe(inlet, junction, "a".to_string());
        builder.connect_with_pipe(junction, outlet, "b".to_string());
        builder.connect_with_pipe(inlet, junction, "c".to_string());
        builder.connect_with_pipe(inlet, outlet, "d".to_string());
        let graph = builder.build().expect("three-node network must validate");
        let network = Network::new(graph, water_20c::<f64>().expect("water properties exist"));

        let incidence = NodeIncidence::from_network(&network);

        assert_eq!(
            &incidence[0],
            &[incident(0, 0, 1), incident(2, 0, 1), incident(3, 0, 2)]
        );
        assert_eq!(
            &incidence[1],
            &[incident(0, 0, 1), incident(1, 1, 2), incident(2, 0, 1)]
        );
        assert_eq!(&incidence[2], &[incident(1, 1, 2), incident(3, 0, 2)]);
    }
}
