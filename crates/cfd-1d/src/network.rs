//! Network topology for 1D CFD simulations

use cfd_core::{Error, Fluid, Result};
use nalgebra::{DVector, RealField};
use num_traits::FromPrimitive;
use petgraph::graph::{EdgeIndex, NodeIndex};
use petgraph::visit::EdgeRef;
use petgraph::{Directed, Graph};
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
pub enum EdgeType {
    Pipe,
    Valve,
    Pump,
/// Node in the network
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct Node<T: RealField + Copy> {
    pub id: String,
    pub node_type: NodeType,
    pub position: (T, T),
impl<T: RealField + Copy> Node<T> {
    #[must_use]
    pub fn new(id: String, node_type: NodeType) -> Self {
        Self {
            id,
            node_type,
            position: (T::zero(), T::zero()),
        }
    }
/// Edge in the network (channel/pipe)
pub struct Edge<T: RealField + Copy> {
    pub edge_type: EdgeType,
    pub properties: ChannelProperties<T>,
/// Channel properties
pub struct ChannelProperties<T: RealField + Copy> {
    pub resistance: T,
    pub length: T,
    pub area: T,
    pub hydraulic_diameter: Option<T>,
impl<T: RealField + Copy> ChannelProperties<T> {
    pub fn new(resistance: T, length: T, area: T) -> Self {
            resistance,
            length,
            area,
            hydraulic_diameter: None,
/// Node properties (for compatibility)
pub type NodeProperties<T> = Node<T>;
/// Network metadata
pub struct NetworkMetadata {
    pub name: String,
    pub description: String,
/// Network builder
pub struct NetworkBuilder<T: RealField + Copy> {
    network: Network<T>,
impl<T: RealField + Copy + FromPrimitive + Copy> NetworkBuilder<T> {
    pub fn new(fluid: Fluid<T>) -> Self {
            network: Network::new(fluid),
    pub fn add_node(mut self, node: Node<T>) -> Self {
        self.network.add_node(node);
        self
    pub fn add_edge(
        mut self,
        from: &str,
        to: &str,
        properties: ChannelProperties<T>,
    ) -> Result<Self> {
        self.network.add_edge(from, to, properties)?;
        Ok(self)
    pub fn build(self) -> Network<T> {
        self.network
/// Main network structure
#[derive(Debug, Clone)]
pub struct Network<T: RealField + Copy> {
    graph: NetworkGraph<Node<T>, ChannelProperties<T>>,
    node_indices: HashMap<String, NodeIndex>,
    edge_indices: HashMap<String, EdgeIndex>,
    boundary_conditions: HashMap<NodeIndex, BoundaryCondition<T>>,
    fluid: Fluid<T>,
    pressures: DVector<T>,
    flow_rates: DVector<T>,
impl<T: RealField + Copy + FromPrimitive + Copy> Network<T> {
            graph: Graph::new(),
            node_indices: HashMap::new(),
            edge_indices: HashMap::new(),
            boundary_conditions: HashMap::new(),
            fluid,
            pressures: DVector::zeros(0),
            flow_rates: DVector::zeros(0),
    pub fn add_node(&mut self, node: Node<T>) -> NodeIndex {
        let id = node.id.clone();
        let idx = self.graph.add_node(node);
        self.node_indices.insert(id, idx);
        idx
        &mut self,
    ) -> Result<EdgeIndex> {
        let from_idx = self
            .node_indices
            .get(from)
            .ok_or_else(|| Error::InvalidConfiguration(format!("Node {from} not found")))?;
        let to_idx = self
            .get(to)
            .ok_or_else(|| Error::InvalidConfiguration(format!("Node {to} not found")))?;
        let edge_id = format!("{from}_{to}");
        let idx = self.graph.add_edge(*from_idx, *to_idx, properties);
        self.edge_indices.insert(edge_id, idx);
        Ok(idx)
    pub fn set_boundary_condition(
        node: &str,
        condition: BoundaryCondition<T>,
    ) -> Result<()> {
        let idx = self
            .get(node)
            .ok_or_else(|| Error::InvalidConfiguration(format!("Node {node} not found")))?;
        self.boundary_conditions.insert(*idx, condition);
        Ok(())
    pub fn node_count(&self) -> usize {
        self.graph.node_count()
    pub fn edge_count(&self) -> usize {
        self.graph.edge_count()
    pub fn characteristic_length(&self) -> T {
        if self.edge_count() == 0 {
            T::one()
        } else {
            let total_length: T = self
                .graph
                .edge_weights()
                .map(|e| e.length)
                .fold(T::zero(), |acc, l| acc + l);
            total_length / T::from_usize(self.edge_count()).unwrap_or_else(T::one)
    pub fn fluid(&self) -> &Fluid<T> {
        &self.fluid
    pub fn pressures(&self) -> &DVector<T> {
        &self.pressures
    pub fn flow_rates(&self) -> &DVector<T> {
        &self.flow_rates
    pub fn update_pressures(&mut self, pressures: &DVector<T>) {
        self.pressures = pressures.clone();
    pub fn update_flow_rates(&mut self, flow_rates: &DVector<T>) {
        self.flow_rates = flow_rates.clone();
    pub fn update_from_solution(&mut self, solution: DVector<T>) -> Result<()> {
        if solution.len() != self.node_count() {
            return Err(Error::InvalidConfiguration(
                "Solution dimension mismatch".to_string(),
            ));
        self.pressures = solution;
        self.calculate_flow_rates()?;
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
        self.flow_rates = DVector::from_vec(flows);
    pub fn boundary_conditions(&self) -> impl Iterator<Item = (usize, &BoundaryCondition<T>)> {
        self.boundary_conditions
            .iter()
            .map(|(idx, bc)| (idx.index(), bc))
    pub fn edges_parallel(&self) -> impl rayon::iter::ParallelIterator<Item = EdgeData<T>> + '_ {
        use rayon::prelude::*;
        self.graph
            .edge_references()
            .par_bridge()
            .map(|edge| EdgeData {
                nodes: (edge.source().index(), edge.target().index()),
                conductance: T::one() / edge.weight().resistance,
            })
    pub fn graph(&self) -> &NetworkGraph<Node<T>, ChannelProperties<T>> {
        &self.graph
    /// Get iterator over edges (returns `EdgeData` for solver use)
    pub fn edges(&self) -> impl Iterator<Item = EdgeData<T>> + '_ {
        self.graph.edge_references().map(|edge| EdgeData {
            nodes: (edge.source().index(), edge.target().index()),
            conductance: T::one() / edge.weight().resistance,
        })
    /// Get iterator over edges with full properties (for analysis)
    pub fn edges_with_properties(&self) -> impl Iterator<Item = EdgeProperties<T>> + '_ {
        self.edge_indices.iter().filter_map(move |(id, &idx)| {
            self.graph.edge_endpoints(idx).map(|(s, t)| {
                let edge = &self.graph[idx];
                EdgeProperties {
                    id: id.clone(),
                    nodes: (s.index(), t.index()),
                    properties: edge.clone(),
                    flow_rate: if idx.index() < self.flow_rates.len() {
                        Some(self.flow_rates[idx.index()])
                    } else {
                        None
                    },
                }
    /// Get iterator over nodes
    pub fn nodes(&self) -> impl Iterator<Item = &Node<T>> + '_ {
        self.graph.node_weights()
    /// Get node index by ID
    pub fn get_node_index(&self, id: &str) -> Option<usize> {
        self.node_indices.get(id).map(|idx| idx.index())
    /// Get node index as `NodeIndex` by ID
    pub fn get_node_index_raw(&self, id: &str) -> Option<NodeIndex> {
        self.node_indices.get(id).copied()
    /// Get edge ID by index
    pub fn get_edge_id_by_index(&self, idx: petgraph::graph::EdgeIndex) -> Option<String> {
        self.edge_indices
            .find(|(_, &edge_idx)| edge_idx == idx)
            .map(|(id, _)| id.clone())
    /// Get edges connected to a node
    pub fn node_edges(&self, node_id: &str) -> Result<Vec<EdgeData<T>>> {
        let node_idx = self
            .get(node_id)
            .ok_or_else(|| Error::InvalidConfiguration(format!("Node {node_id} not found")))?;
        Ok(self
            .graph
            .edges(*node_idx)
            .collect())
    pub fn node_indices(&self) -> &HashMap<String, NodeIndex> {
        &self.node_indices
pub struct EdgeData<T: RealField + Copy> {
    pub nodes: (usize, usize),
    pub conductance: T,
/// Edge with full properties for analysis
pub struct EdgeProperties<T: RealField + Copy> {
    pub flow_rate: Option<T>,
