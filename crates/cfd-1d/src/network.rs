//! Network topology for 1D CFD simulations

use cfd_core::{Error, Result, Fluid};
use nalgebra::{RealField, DVector};
use num_traits::FromPrimitive;
use petgraph::{Graph, Directed};
use petgraph::graph::{EdgeIndex, NodeIndex};
use petgraph::visit::EdgeRef;
use serde::{Deserialize, Serialize};
use std::collections::HashMap;

pub type BoundaryCondition<T> = cfd_core::boundary::BoundaryCondition<T>;
pub type NetworkGraph<N, E> = Graph<N, E, Directed>;

/// Node types in the network
#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize)]
pub enum NodeType {
    Inlet,
    Outlet,
    Junction,
    Reservoir,
}

/// Edge types in the network
#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize)]
pub enum EdgeType {
    Pipe,
    Valve,
    Pump,
}

/// Node in the network
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct Node<T: RealField> {
    pub id: String,
    pub node_type: NodeType,
    pub position: (T, T),
}

impl<T: RealField> Node<T> {
    pub fn new(id: String, node_type: NodeType) -> Self {
        Self {
            id,
            node_type,
            position: (T::zero(), T::zero()),
        }
    }
}

/// Edge in the network (channel/pipe)
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct Edge<T: RealField> {
    pub id: String,
    pub edge_type: EdgeType,
    pub properties: ChannelProperties<T>,
}

/// Channel properties
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct ChannelProperties<T: RealField> {
    pub resistance: T,
    pub length: T,
    pub area: T,
}

impl<T: RealField> ChannelProperties<T> {
    pub fn new(resistance: T, length: T, area: T) -> Self {
        Self { resistance, length, area }
    }
}

/// Node properties (for compatibility)
pub type NodeProperties<T> = Node<T>;

/// Edge properties (for compatibility)
pub type EdgeProperties<T> = ChannelProperties<T>;

/// Network metadata
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct NetworkMetadata {
    pub name: String,
    pub description: String,
}

/// Network builder
pub struct NetworkBuilder<T: RealField> {
    network: Network<T>,
}

impl<T: RealField + FromPrimitive + Copy> NetworkBuilder<T> {
    pub fn new(fluid: Fluid<T>) -> Self {
        Self {
            network: Network::new(fluid),
        }
    }

    pub fn add_node(mut self, node: Node<T>) -> Self {
        self.network.add_node(node);
        self
    }

    pub fn add_edge(mut self, from: &str, to: &str, properties: ChannelProperties<T>) -> Result<Self> {
        self.network.add_edge(from, to, properties)?;
        Ok(self)
    }

    pub fn build(self) -> Network<T> {
        self.network
    }
}

/// Main network structure
#[derive(Debug, Clone)]
pub struct Network<T: RealField> {
    graph: NetworkGraph<Node<T>, ChannelProperties<T>>,
    node_indices: HashMap<String, NodeIndex>,
    edge_indices: HashMap<String, EdgeIndex>,
    boundary_conditions: HashMap<NodeIndex, BoundaryCondition<T>>,
    fluid: Fluid<T>,
    pressures: DVector<T>,
    flow_rates: DVector<T>,
}

impl<T: RealField + FromPrimitive + Copy> Network<T> {
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

    pub fn add_node(&mut self, node: Node<T>) -> NodeIndex {
        let id = node.id.clone();
        let idx = self.graph.add_node(node);
        self.node_indices.insert(id, idx);
        idx
    }

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

    pub fn set_boundary_condition(&mut self, node: &str, condition: BoundaryCondition<T>) -> Result<()> {
        let idx = self.node_indices.get(node)
            .ok_or_else(|| Error::InvalidConfiguration(format!("Node {} not found", node)))?;
        self.boundary_conditions.insert(*idx, condition);
        Ok(())
    }

    pub fn node_count(&self) -> usize {
        self.graph.node_count()
    }

    pub fn edge_count(&self) -> usize {
        self.graph.edge_count()
    }

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

    pub fn fluid(&self) -> &Fluid<T> {
        &self.fluid
    }

    pub fn pressures(&self) -> &DVector<T> {
        &self.pressures
    }

    pub fn flow_rates(&self) -> &DVector<T> {
        &self.flow_rates
    }

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

    fn calculate_flow_rates(&mut self) -> Result<()> {
        let mut flows = Vec::with_capacity(self.edge_count());
        
        for edge in self.graph.edge_references() {
            let from_idx = edge.source();
            let to_idx = edge.target();
            let properties = edge.weight();
            
            let p1 = self.pressures[from_idx.index()];
            let p2 = self.pressures[to_idx.index()];
            
            let flow = (p1 - p2) / properties.resistance;
            flows.push(flow);
        }
        
        self.flow_rates = DVector::from_vec(flows);
        Ok(())
    }

    pub fn boundary_conditions(&self) -> impl Iterator<Item = (usize, &BoundaryCondition<T>)> {
        self.boundary_conditions.iter()
            .map(|(idx, bc)| (idx.index(), bc))
    }

    pub fn edges_parallel(&self) -> impl rayon::iter::ParallelIterator<Item = EdgeData<T>> + '_ {
        use rayon::prelude::*;
        
        self.graph.edge_references()
            .par_bridge()
            .map(|edge| EdgeData {
                nodes: (edge.source().index(), edge.target().index()),
                conductance: T::one() / edge.weight().resistance.clone(),
            })
    }
    
    pub fn graph(&self) -> &NetworkGraph<Node<T>, ChannelProperties<T>> {
        &self.graph
    }
    
    /// Get iterator over edges
    pub fn edges(&self) -> impl Iterator<Item = EdgeData<T>> + '_ {
        self.graph.edge_references()
            .map(|edge| EdgeData {
                nodes: (edge.source().index(), edge.target().index()),
                conductance: T::one() / edge.weight().resistance.clone(),
            })
    }
    
    /// Get iterator over nodes
    pub fn nodes(&self) -> impl Iterator<Item = &Node<T>> + '_ {
        self.graph.node_weights()
    }
    
    /// Get node index by ID
    pub fn get_node_index(&self, id: &str) -> Option<usize> {
        self.node_indices.get(id).map(|idx| idx.index())
    }
    
    /// Get edge ID by index
    pub fn get_edge_id_by_index(&self, idx: petgraph::graph::EdgeIndex) -> Option<String> {
        self.edge_indices.iter()
            .find(|(_, &edge_idx)| edge_idx == idx)
            .map(|(id, _)| id.clone())
    }
    
    /// Get edges connected to a node
    pub fn node_edges(&self, node_id: &str) -> Result<Vec<EdgeData<T>>> {
        let node_idx = self.node_indices.get(node_id)
            .ok_or_else(|| Error::InvalidConfiguration(format!("Node {} not found", node_id)))?;
        
        Ok(self.graph.edges(*node_idx)
            .map(|edge| EdgeData {
                nodes: (edge.source().index(), edge.target().index()),
                conductance: T::one() / edge.weight().resistance.clone(),
            })
            .collect())
    }
    
    pub fn node_indices(&self) -> &HashMap<String, NodeIndex> {
        &self.node_indices
    }
}

pub struct EdgeData<T: RealField> {
    pub nodes: (usize, usize),
    pub conductance: T,
}