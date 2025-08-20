//! Network topology and graph structures for 1D CFD simulations.
//!
//! This module provides the core network representation for microfluidic and millifluidic
//! systems using graph-based data structures optimized for CFD analysis.

mod node;
mod edge;
mod builder;
mod analysis;
mod validation;

pub use node::{Node, NodeType};
pub use edge::{ChannelProperties, Edge};
pub use builder::NetworkBuilder;
pub use analysis::NetworkAnalysis;
pub use validation::NetworkValidator;

use cfd_core::{Error, Result, Fluid};
use nalgebra::{RealField, DVector};
use num_traits::FromPrimitive;
use petgraph::{Graph, Directed};
use petgraph::graph::{EdgeIndex, NodeIndex};
use serde::{Deserialize, Serialize};
use std::collections::HashMap;

/// Type alias for boundary conditions in 1D networks
pub type BoundaryCondition<T> = cfd_core::boundary::BoundaryCondition<T>;

/// Network graph type using petgraph
pub type NetworkGraph<N, E> = Graph<N, E, Directed>;

/// Main network structure for 1D CFD analysis
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct Network<T: RealField> {
    /// Graph structure
    graph: NetworkGraph<Node<T>, ChannelProperties<T>>,
    /// Node index mapping
    node_indices: HashMap<String, NodeIndex>,
    /// Edge index mapping
    edge_indices: HashMap<String, EdgeIndex>,
    /// Boundary conditions
    boundary_conditions: HashMap<NodeIndex, BoundaryCondition<T>>,
    /// Fluid properties
    fluid: Fluid<T>,
    /// Solution state - pressures at nodes
    pressures: DVector<T>,
    /// Solution state - flow rates at edges
    flow_rates: DVector<T>,
}

impl<T: RealField + FromPrimitive> Network<T> {
    /// Create a new empty network
    pub fn new(fluid: Fluid<T>) -> Self {
        Self {
            graph: Graph::new(),
            node_indices: HashMap::new(),
            edge_indices: HashMap::new(),
            boundary_conditions: HashMap::new(),
            fluid,
            pressures: DVector::zeros(0),
            flow_rates: DVector::zeros(0),
        }
    }

    /// Add a node to the network
    pub fn add_node(&mut self, node: Node<T>) -> NodeIndex {
        let id = node.id.clone();
        let idx = self.graph.add_node(node);
        self.node_indices.insert(id, idx);
        idx
    }

    /// Add an edge (channel) to the network
    pub fn add_edge(&mut self, from: &str, to: &str, properties: ChannelProperties<T>) -> Result<EdgeIndex> {
        let from_idx = self.node_indices.get(from)
            .ok_or_else(|| Error::InvalidConfiguration(format!("Node {} not found", from)))?;
        let to_idx = self.node_indices.get(to)
            .ok_or_else(|| Error::InvalidConfiguration(format!("Node {} not found", to)))?;
        
        let edge_id = format!("{}_{}", from, to);
        let idx = self.graph.add_edge(*from_idx, *to_idx, properties);
        self.edge_indices.insert(edge_id, idx);
        Ok(idx)
    }

    /// Set boundary condition at a node
    pub fn set_boundary_condition(&mut self, node: &str, condition: BoundaryCondition<T>) -> Result<()> {
        let idx = self.node_indices.get(node)
            .ok_or_else(|| Error::InvalidConfiguration(format!("Node {} not found", node)))?;
        self.boundary_conditions.insert(*idx, condition);
        Ok(())
    }

    /// Get the number of nodes
    pub fn node_count(&self) -> usize {
        self.graph.node_count()
    }

    /// Get the number of edges
    pub fn edge_count(&self) -> usize {
        self.graph.edge_count()
    }

    /// Get characteristic length scale
    pub fn characteristic_length(&self) -> T {
        if self.edge_count() == 0 {
            T::one()
        } else {
            let total_length: T = self.graph.edge_weights()
                .map(|e| e.length)
                .fold(T::zero(), |acc, l| acc + l);
            total_length / T::from_usize(self.edge_count()).unwrap_or_else(T::one)
        }
    }

    /// Get fluid properties
    pub fn fluid(&self) -> &Fluid<T> {
        &self.fluid
    }

    /// Get pressures vector
    pub fn pressures(&self) -> &DVector<T> {
        &self.pressures
    }

    /// Get flow rates vector
    pub fn flow_rates(&self) -> &DVector<T> {
        &self.flow_rates
    }

    /// Update solution from solver
    pub fn update_from_solution(&mut self, solution: DVector<T>) -> Result<()> {
        if solution.len() != self.node_count() {
            return Err(Error::InvalidConfiguration(
                "Solution dimension mismatch".to_string()
            ));
        }
        self.pressures = solution;
        self.calculate_flow_rates()?;
        Ok(())
    }

    /// Calculate flow rates from pressures
    fn calculate_flow_rates(&mut self) -> Result<()> {
        let mut flows = Vec::with_capacity(self.edge_count());
        
        for edge in self.graph.edge_references() {
            let from_idx = edge.source();
            let to_idx = edge.target();
            let properties = edge.weight();
            
            let p1 = self.pressures[from_idx.index()];
            let p2 = self.pressures[to_idx.index()];
            
            // Q = ΔP / R (Hagen-Poiseuille for laminar flow)
            let flow = (p1 - p2) / properties.resistance;
            flows.push(flow);
        }
        
        self.flow_rates = DVector::from_vec(flows);
        Ok(())
    }

    /// Get boundary conditions iterator
    pub fn boundary_conditions(&self) -> impl Iterator<Item = (usize, &BoundaryCondition<T>)> {
        self.boundary_conditions.iter()
            .map(|(idx, bc)| (idx.index(), bc))
    }

    /// Get edges for parallel processing
    pub fn edges_parallel(&self) -> impl ParallelIterator<Item = EdgeData<T>> + '_ {
        use rayon::prelude::*;
        
        self.graph.edge_references()
            .par_bridge()
            .map(|edge| EdgeData {
                nodes: (edge.source().index(), edge.target().index()),
                conductance: T::one() / edge.weight().resistance,
            })
    }
}

/// Edge data for parallel processing
pub struct EdgeData<T: RealField> {
    pub nodes: (usize, usize),
    pub conductance: T,
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_network_creation() {
        let fluid = Fluid::water().expect("Failed to create water fluid");
        let mut network = Network::<f64>::new(fluid);
        
        let node1 = Node::new("inlet".to_string(), NodeType::Inlet);
        let node2 = Node::new("outlet".to_string(), NodeType::Outlet);
        
        network.add_node(node1);
        network.add_node(node2);
        
        assert_eq!(network.node_count(), 2);
    }
}