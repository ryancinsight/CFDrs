//! Network wrapper with convenience methods

use super::{NetworkGraph, Node};
use crate::channel::ChannelGeometry;
use cfd_core::{
    error::{Error, Result},
    physics::boundary::BoundaryCondition,
    physics::fluid::{ConstantPropertyFluid, FluidTrait},
};
use nalgebra::{DVector, RealField};
use num_traits::FromPrimitive;
use petgraph::graph::{EdgeIndex, NodeIndex};
use petgraph::visit::EdgeRef;
use std::collections::HashMap;

/// Extended network with fluid properties and convenience methods
#[derive(Clone, Debug)]
pub struct Network<T: RealField + Copy, F: FluidTrait<T> = ConstantPropertyFluid<T>> {
    /// The underlying graph
    pub graph: NetworkGraph<T>,
    /// Fluid properties for the network
    pub fluid: F,
    /// Prescribed boundary conditions
    boundary_conditions: HashMap<NodeIndex, BoundaryCondition<T>>,
    /// Node pressures
    pub pressures: HashMap<NodeIndex, T>,
    /// Edge flow rates
    pub flow_rates: HashMap<EdgeIndex, T>,
    /// Edge properties
    pub properties: HashMap<EdgeIndex, EdgeProperties<T>>,
    /// Residuals from the last solver run
    pub residuals: Vec<T>,
    /// Last used solver method
    pub last_solver_method: Option<crate::solver::LinearSolverMethod>,
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

impl<T: RealField + Copy + FromPrimitive, F: FluidTrait<T>> Network<T, F> {
    /// Create a new network
    pub fn new(graph: NetworkGraph<T>, fluid: F) -> Self {
        Self {
            graph,
            fluid,
            boundary_conditions: HashMap::new(),
            pressures: HashMap::new(),
            flow_rates: HashMap::new(),
            properties: HashMap::new(),
            residuals: Vec::new(),
            last_solver_method: None,
        }
    }

    /// Get the fluid properties
    pub fn fluid(&self) -> &F {
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

    /// Get residuals from the last solver run
    pub fn residuals(&self) -> &Vec<T> {
        &self.residuals
    }

    /// Set a Dirichlet pressure boundary condition
    pub fn set_pressure(&mut self, node: NodeIndex, pressure: T) {
        self.boundary_conditions.insert(
            node,
            BoundaryCondition::Dirichlet {
                value: pressure,
                component_values: None,
            },
        );
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
        if let BoundaryCondition::Dirichlet { value, .. } = &condition {
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

            let (resistance, quad_coeff) = self
                .graph
                .edge_weight(edge_idx)
                .map(|edge| (edge.resistance, edge.quad_coeff))
                .ok_or_else(|| Error::InvalidConfiguration("Missing edge data".into()))?;

            // Invariant checks: physical coefficients must be non-negative
            if resistance < T::zero() {
                return Err(Error::InvalidConfiguration(format!(
                    "Edge {} has negative resistance: {}",
                    edge_idx.index(),
                    resistance
                )));
            }
            if quad_coeff < T::zero() {
                return Err(Error::InvalidConfiguration(format!(
                    "Edge {} has negative quadratic coefficient: {}",
                    edge_idx.index(),
                    quad_coeff
                )));
            }

            if resistance.abs() < epsilon && quad_coeff.abs() < epsilon {
                return Err(Error::InvalidConfiguration(format!(
                    "Edge {} has zero resistance and zero quadratic coefficient, cannot infer flow",
                    edge_idx.index()
                )));
            }

            let p_from = solution[from.index()];
            let p_to = solution[to.index()];
            let dp = p_from - p_to;
            let dp_abs = dp.abs();
            let sign = if dp >= T::zero() { T::one() } else { -T::one() };

            // Calculate flow rate solving the quadratic equation:
            // |dP| = R*|Q| + k*|Q|^2
            // k*|Q|^2 + R*|Q| - |dP| = 0
            // |Q| = (-R + sqrt(R^2 + 4*k*|dP|)) / (2*k)

            let flow = if quad_coeff.abs() < epsilon {
                // Linear case: Q = dP / R
                dp / resistance
            } else {
                // Quadratic case
                let two = T::from_f64(2.0).expect("Field must support 2.0");
                let four = T::from_f64(4.0).expect("Field must support 4.0");

                let discriminant = resistance * resistance + four * quad_coeff * dp_abs;
                let sqrt_disc = discriminant.sqrt();

                // Numerator: -R + sqrt(Delta)
                // Denominator: 2k
                let q_mag = (sqrt_disc - resistance) / (two * quad_coeff);
                sign * q_mag
            };

            self.flow_rates.insert(edge_idx, flow);

            if let Some(edge) = self.graph.edge_weight_mut(edge_idx) {
                edge.flow_rate = flow;
            }
        }

        // After updating flow rates, re-calculate resistances and quadratic coefficients
        // if they are flow-rate dependent.
        self.update_resistances()?;

        Ok(())
    }

    /// Re-calculate resistances and quadratic coefficients for all edges based on current flow rates.
    pub fn update_resistances(&mut self) -> Result<()> {
        let edge_indices: Vec<_> = self.graph.edge_indices().collect();
        let calculator = crate::resistance::ResistanceCalculator::new();

        for edge_idx in edge_indices {
            let flow_rate = *self.flow_rates.get(&edge_idx).unwrap_or(&T::zero());

            // Get geometry and other properties from self.properties
            if let Some(props) = self.properties.get(&edge_idx) {
                if let Some(geometry) = &props.geometry {
                    // Map channel::ChannelGeometry to resistance::ChannelGeometry
                    let res_geometry = match &geometry.cross_section {
                        crate::channel::CrossSection::Circular { diameter } => {
                            crate::resistance::ChannelGeometry::Circular {
                                diameter: *diameter,
                                length: geometry.length,
                            }
                        }
                        crate::channel::CrossSection::Rectangular { width, height } => {
                            crate::resistance::ChannelGeometry::Rectangular {
                                width: *width,
                                height: *height,
                                length: geometry.length,
                            }
                        }
                        crate::channel::CrossSection::Elliptical {
                            major_axis,
                            minor_axis,
                        } => crate::resistance::ChannelGeometry::Elliptical {
                            major_axis: *major_axis,
                            minor_axis: *minor_axis,
                            length: geometry.length,
                        },
                        crate::channel::CrossSection::Trapezoidal {
                            top_width,
                            bottom_width,
                            height,
                        } => crate::resistance::ChannelGeometry::Trapezoidal {
                            top_width: *top_width,
                            bottom_width: *bottom_width,
                            height: *height,
                            length: geometry.length,
                        },
                        crate::channel::CrossSection::Custom {
                            area,
                            hydraulic_diameter,
                        } => crate::resistance::ChannelGeometry::Custom {
                            area: *area,
                            hydraulic_diameter: *hydraulic_diameter,
                            length: geometry.length,
                        },
                    };

                    // Prepare flow conditions
                    // Initialize with zero velocity, but clear it to allow calculator to derive it from flow rate
                    let mut conditions = crate::resistance::FlowConditions::new(T::zero());
                    conditions.velocity = None;
                    conditions.flow_rate = Some(flow_rate.abs());

                    // Re-calculate coefficients
                    let (r, k) = calculator.calculate_coefficients_auto(
                        &res_geometry,
                        &self.fluid,
                        &conditions,
                    )?;

                    // Update the edge in the graph
                    if let Some(edge) = self.graph.edge_weight_mut(edge_idx) {
                        edge.resistance = r;
                        edge.quad_coeff = k;
                    }
                }
            }
        }
        Ok(())
    }

    /// Validate network edge coefficients (resistance and quadratic coefficients).
    ///
    /// Returns an error if any edge has negative resistance or quadratic coefficient.
    pub fn validate_coefficients(&self) -> Result<()> {
        let eps = T::default_epsilon();
        for edge_ref in self.graph.edge_references() {
            let idx = edge_ref.id();
            let w = edge_ref.weight();
            let r = w.resistance;
            let k = w.quad_coeff;
            if r < T::zero() {
                return Err(Error::InvalidConfiguration(format!(
                    "Edge {} has negative resistance: {}",
                    idx.index(),
                    r
                )));
            }
            if k < T::zero() {
                return Err(Error::InvalidConfiguration(format!(
                    "Edge {} has negative quadratic coefficient: {}",
                    idx.index(),
                    k
                )));
            }
            if r.abs() < eps && k.abs() < eps {
                return Err(Error::InvalidConfiguration(format!(
                    "Edge {} has zero resistance and zero quadratic coefficient",
                    idx.index()
                )));
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
            // Derivation: For quadratic loss ΔP = R·Q + k·Q|Q|, the local linearization around Q_k is:
            // ΔP ≈ (R + 2k|Q_k|)·Q.
            // Thus, effective resistance R_eff = R + 2k|Q_k|.
            // We use 2.0 explicitly; failure to represent 2.0 is a critical system failure.
            let two = T::from_f64(2.0).expect("Field must support value 2.0");
            let r_eff = edge_data.resistance + two * edge_data.quad_coeff * q;
            let eps = T::default_epsilon();
            // Conductance is 1/R_eff. We enforce R_eff > ε to avoid division by zero.
            // If R_eff is effectively zero or invalid, we return zero conductance (infinite resistance).
            let conductance = if r_eff > eps && r_eff.is_finite() {
                r_eff.recip()
            } else {
                T::zero()
            };

            ParallelEdge {
                nodes: (from.index(), to.index()),
                conductance,
            }
        })
    }
}
