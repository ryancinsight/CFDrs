//! Network wrapper with convenience methods

use super::{NetworkGraph, Node};
use crate::channel::ChannelGeometry;
use cfd_core::{
    boundary::BoundaryCondition,
    error::{Error, Result},
    fluid::{ConstantPropertyFluid, Fluid},
};
use nalgebra::{DVector, RealField};
use num_traits::FromPrimitive;
use petgraph::graph::{EdgeIndex, NodeIndex};
use petgraph::visit::EdgeRef;
use std::collections::HashMap;

/// Extended network with fluid properties and convenience methods
#[derive(Clone, Debug)]
pub struct Network<T: RealField + Copy> {
    /// The underlying graph
    pub graph: NetworkGraph<T>,
    /// Fluid properties for the network
    pub fluid: ConstantPropertyFluid<T>,
    /// Prescribed boundary conditions
    boundary_conditions: HashMap<NodeIndex, BoundaryCondition<T>>,
    /// Node pressures
    pub pressures: HashMap<NodeIndex, T>,
    /// Edge flow rates
    pub flow_rates: HashMap<EdgeIndex, T>,
    /// Edge properties
    pub properties: HashMap<EdgeIndex, EdgeProperties<T>>,
}

/// Properties for edges in the network
#[derive(Debug, Clone)]
pub struct EdgeProperties<T: RealField + Copy> {
    /// Edge identifier
    pub id: String,
    /// Physical component type
    pub component_type: super::ComponentType,
    /// Length of the channel
    pub length: T,
    /// Cross-sectional area
    pub area: T,
    /// Hydraulic diameter
    pub hydraulic_diameter: Option<T>,
    /// Resistance coefficient
    pub resistance: T,
    /// Channel geometry if applicable
    pub geometry: Option<ChannelGeometry<T>>,
    /// Additional properties
    pub properties: HashMap<String, T>,
}

/// Edge with properties for iteration
pub struct EdgeWithProperties<'a, T: RealField + Copy> {
    /// Edge identifier
    pub id: String,
    /// Flow rate
    pub flow_rate: T,
    /// Node indices (from, to)
    pub nodes: (NodeIndex, NodeIndex),
    /// Properties
    pub properties: &'a EdgeProperties<T>,
}

/// Edge for parallel processing
pub struct ParallelEdge<T: RealField + Copy> {
    /// Node indices (from, to) as usize
    pub nodes: (usize, usize),
    /// Conductance (1/resistance)
    pub conductance: T,
}

impl<T: RealField + Copy + FromPrimitive> Network<T> {
    /// Create a new network
    pub fn new(graph: NetworkGraph<T>, fluid: ConstantPropertyFluid<T>) -> Self {
        Self {
            graph,
            fluid,
            boundary_conditions: HashMap::new(),
            pressures: HashMap::new(),
            flow_rates: HashMap::new(),
            properties: HashMap::new(),
        }
    }

    /// Get the fluid properties
    pub fn fluid(&self) -> &Fluid<T> {
        &self.fluid
    }

    /// Get edges with their properties
    pub fn edges_with_properties(&self) -> Vec<EdgeWithProperties<'_, T>> {
        self.graph
            .edge_references()
            .filter_map(|edge_ref| {
                let edge_idx = edge_ref.id();
                let edge_data = edge_ref.weight();
                let (from, to) = (edge_ref.source(), edge_ref.target());
                self.properties
                    .get(&edge_idx)
                    .map(|props| EdgeWithProperties {
                        id: edge_data.id.clone(),
                        flow_rate: *self.flow_rates.get(&edge_idx).unwrap_or(&T::zero()),
                        nodes: (from, to),
                        properties: props,
                    })
            })
            .collect()
    }

    /// Get all nodes
    pub fn nodes(&self) -> impl Iterator<Item = &Node<T>> {
        self.graph.node_weights()
    }

    /// Get node pressures
    pub fn pressures(&self) -> &HashMap<NodeIndex, T> {
        &self.pressures
    }

    /// Get edge flow rates
    pub fn flow_rates(&self) -> &HashMap<EdgeIndex, T> {
        &self.flow_rates
    }

    /// Set a Dirichlet pressure boundary condition
    pub fn set_pressure(&mut self, node: NodeIndex, pressure: T) {
        self.boundary_conditions
            .insert(node, BoundaryCondition::Dirichlet { value: pressure });
        self.pressures.insert(node, pressure);
    }

    /// Set a Neumann (flow) boundary condition
    pub fn set_neumann_flow(&mut self, node: NodeIndex, flow_rate: T) {
        self.boundary_conditions.insert(
            node,
            BoundaryCondition::Neumann {
                gradient: flow_rate,
            },
        );
    }

    /// Explicitly set a boundary condition
    pub fn set_boundary_condition(&mut self, node: NodeIndex, condition: BoundaryCondition<T>) {
        if let BoundaryCondition::Dirichlet { value } = &condition {
            self.pressures.insert(node, *value);
        }
        self.boundary_conditions.insert(node, condition);
    }

    /// Set flow rate for an edge
    pub fn set_flow_rate(&mut self, edge: EdgeIndex, flow_rate: T) {
        self.flow_rates.insert(edge, flow_rate);
        if let Some(edge_data) = self.graph.edge_weight_mut(edge) {
            edge_data.flow_rate = flow_rate;
        }
    }

    /// Add properties for an edge
    pub fn add_edge_properties(&mut self, edge: EdgeIndex, properties: EdgeProperties<T>) {
        self.properties.insert(edge, properties);
    }

    /// Get the number of nodes in the network
    pub fn node_count(&self) -> usize {
        self.graph.node_count()
    }

    /// Get the number of edges in the network
    pub fn edge_count(&self) -> usize {
        self.graph.edge_count()
    }

    /// Get characteristic length of the network
    pub fn characteristic_length(&self) -> T {
        // Calculate based on average edge length
        if self.properties.is_empty() {
            T::one()
        } else {
            let total_length: T = self
                .properties
                .values()
                .map(|p| p.length)
                .fold(T::zero(), |a, b| a + b);
            total_length / T::from_usize(self.properties.len()).unwrap_or(T::one())
        }
    }

    /// Update network state from solution vector
    pub fn update_from_solution(&mut self, solution: &DVector<T>) -> Result<()> {
        self.pressures.clear();
        for (idx, &value) in solution.iter().enumerate() {
            let node_idx = NodeIndex::new(idx);
            self.pressures.insert(node_idx, value);
        }

        self.flow_rates.clear();
        let epsilon = T::default_epsilon();
        let edge_indices: Vec<_> = self.graph.edge_indices().collect();

        for edge_idx in edge_indices {
            let (from, to) = self
                .graph
                .edge_endpoints(edge_idx)
                .ok_or_else(|| Error::InvalidConfiguration("Missing edge endpoints".into()))?;
            let resistance = self
                .graph
                .edge_weight(edge_idx)
                .map(|edge| edge.resistance)
                .ok_or_else(|| Error::InvalidConfiguration("Missing edge resistance".into()))?;

            if resistance.abs() < epsilon {
                return Err(Error::InvalidConfiguration(format!(
                    "Edge {} has zero resistance, cannot infer flow",
                    edge_idx.index()
                )));
            }

            let p_from = solution[from.index()];
            let p_to = solution[to.index()];
            let flow = (p_from - p_to) / resistance;

            self.flow_rates.insert(edge_idx, flow);

            if let Some(edge) = self.graph.edge_weight_mut(edge_idx) {
                edge.flow_rate = flow;
            }
        }

        Ok(())
    }

    /// Get boundary conditions for the network
    pub fn boundary_conditions(&self) -> &HashMap<NodeIndex, BoundaryCondition<T>> {
        &self.boundary_conditions
    }

    /// Process edges in parallel
    pub fn edges_parallel(&self) -> impl Iterator<Item = ParallelEdge<T>> + '_ {
        self.graph.edge_references().map(move |edge_ref| {
            let (from, to) = (edge_ref.source(), edge_ref.target());
            let edge_data = edge_ref.weight();

            let q = edge_data.flow_rate.abs();
            let two = T::from_f64(2.0).unwrap_or_else(|| T::one());
            let r_eff = edge_data.resistance + two * edge_data.quad_coeff * q;

            ParallelEdge {
                nodes: (from.index(), to.index()),
                conductance: r_eff.recip(),
            }
        })
    }
}
