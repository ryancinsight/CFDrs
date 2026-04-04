//! Network wrapper with convenience methods

use super::{NetworkGraph, Node};
use crate::domain::channel::ChannelGeometry;
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
    /// Node pressures indexed by `NodeIndex::index()`.
    pub pressures: Vec<T>,
    /// Edge flow rates indexed by `EdgeIndex::index()`.
    pub flow_rates: Vec<T>,
    /// Edge properties
    pub properties: HashMap<EdgeIndex, EdgeProperties<T>>,
    /// Residuals from the last solver run
    pub residuals: Vec<T>,
    /// Last used solver method
    pub last_solver_method: Option<crate::solver::core::LinearSolverMethod>,
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
    /// Whether the hydraulic coefficients must be recomputed as flow changes.
    pub resistance_update_policy: ResistanceUpdatePolicy,
    /// Additional properties
    pub properties: HashMap<String, T>,
}

/// Update policy for per-edge hydraulic coefficients during Picard iteration.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum ResistanceUpdatePolicy {
    /// Resistance law is flow-invariant for the selected reduced-order model.
    FlowInvariant,
    /// Resistance law depends on the current flow state and must be refreshed.
    FlowDependent,
}

use crate::domain::channel::{ChannelType, CrossSection, SurfaceProperties, Wettability};
use cfd_schematics::domain::model::{ChannelSpec, CrossSectionSpec};

impl<T: RealField + Copy + FromPrimitive> From<&ChannelSpec> for EdgeProperties<T> {
    /// Convert a `ChannelSpec` from `cfd-schematics` into solver-layer `EdgeProperties`.
    ///
    /// This is the canonical bridge between the schematic domain model and the
    /// 1D solver. It eliminates the need for examples to import `cfd_1d::channel`
    /// types directly.
    fn from(spec: &ChannelSpec) -> Self {
        let length =
            T::from_f64(spec.length_m).expect("Mathematical constant conversion compromised");
        let resistance =
            T::from_f64(spec.resistance).expect("Mathematical constant conversion compromised");

        let (cross_section, area, hydraulic_diameter) = match spec.cross_section {
            CrossSectionSpec::Circular { diameter_m } => {
                let d =
                    T::from_f64(diameter_m).expect("Mathematical constant conversion compromised");
                let a = T::from_f64(std::f64::consts::PI * (diameter_m / 2.0).powi(2))
                    .unwrap_or(T::zero());
                (CrossSection::Circular { diameter: d }, a, Some(d))
            }
            CrossSectionSpec::Rectangular { width_m, height_m } => {
                let w = T::from_f64(width_m).expect("Mathematical constant conversion compromised");
                let h =
                    T::from_f64(height_m).expect("Mathematical constant conversion compromised");
                let a = T::from_f64(width_m * height_m)
                    .expect("Mathematical constant conversion compromised");
                let dh = T::from_f64(2.0 * width_m * height_m / (width_m + height_m))
                    .unwrap_or(T::zero());
                (
                    CrossSection::Rectangular {
                        width: w,
                        height: h,
                    },
                    a,
                    Some(dh),
                )
            }
        };

        let geometry = Some(ChannelGeometry {
            channel_type: ChannelType::Straight,
            length,
            cross_section,
            surface: SurfaceProperties {
                roughness: T::zero(),
                contact_angle: None,
                surface_energy: None,
                wettability: Wettability::Hydrophilic,
            },
            variations: Vec::new(),
        });

        let component_type = match spec.kind {
            cfd_schematics::domain::model::EdgeKind::Pipe => super::ComponentType::Pipe,
            cfd_schematics::domain::model::EdgeKind::Valve => super::ComponentType::Valve,
            cfd_schematics::domain::model::EdgeKind::Pump => super::ComponentType::Pump,
        };

        Self {
            id: spec.id.as_str().to_string(),
            component_type,
            length,
            area,
            hydraulic_diameter,
            resistance,
            geometry,
            resistance_update_policy: ResistanceUpdatePolicy::FlowDependent,
            properties: HashMap::new(),
        }
    }
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
        let n_nodes = graph.node_count();
        let n_edges = graph.edge_count();
        Self {
            graph,
            fluid,
            boundary_conditions: HashMap::new(),
            pressures: vec![T::zero(); n_nodes],
            flow_rates: vec![T::zero(); n_edges],
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
                        flow_rate: self.flow_rates.get(edge_idx.index()).copied().unwrap_or(T::zero()),
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

    /// Get node pressures as a slice indexed by `NodeIndex::index()`.
    pub fn pressures(&self) -> &[T] {
        &self.pressures
    }

    /// Get edge flow rates as a slice indexed by `EdgeIndex::index()`.
    pub fn flow_rates(&self) -> &[T] {
        &self.flow_rates
    }

    /// Get residuals from the last solver run
    pub fn residuals(&self) -> &[T] {
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
        let idx = node.index();
        if idx >= self.pressures.len() {
            self.pressures.resize(idx + 1, T::zero());
        }
        self.pressures[idx] = pressure;
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
            let idx = node.index();
            if idx >= self.pressures.len() {
                self.pressures.resize(idx + 1, T::zero());
            }
            self.pressures[idx] = *value;
        }
        self.boundary_conditions.insert(node, condition);
    }

    /// Set flow rate for an edge
    pub fn set_flow_rate(&mut self, edge: EdgeIndex, flow_rate: T) {
        let idx = edge.index();
        if idx >= self.flow_rates.len() {
            self.flow_rates.resize(idx + 1, T::zero());
        }
        self.flow_rates[idx] = flow_rate;
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
        // Direct copy: pressures[i] = solution[i]
        let n = solution.len();
        self.pressures.resize(n, T::zero());
        for i in 0..n {
            self.pressures[i] = solution[i];
        }

        let epsilon = T::default_epsilon();
        let n_edges = self.graph.edge_count();
        self.flow_rates.resize(n_edges, T::zero());
        let edge_indices: Vec<_> = self.graph.edge_indices().collect();

        for edge_idx in edge_indices {
            let (from, to) = self
                .graph
                .edge_endpoints(edge_idx)
                .ok_or_else(|| Error::InvalidConfiguration("Missing edge endpoints".into()))?;

            let (resistance, quad_coeff, previous_flow_rate) = self
                .graph
                .edge_weight(edge_idx)
                .map(|edge| (edge.resistance, edge.quad_coeff, edge.flow_rate.abs()))
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
            // Using Picard iteration (secant modulus) avoids the need for a residual
            // offset vector. R_eff = R + k|Q_prev|.
            let r_eff = resistance + quad_coeff * previous_flow_rate;
            if !r_eff.is_finite() || r_eff <= epsilon {
                return Err(Error::InvalidConfiguration(format!(
                    "Edge {} has invalid effective resistance after linearization: {}",
                    edge_idx.index(),
                    r_eff
                )));
            }
            let flow = dp / r_eff;

            self.flow_rates[edge_idx.index()] = flow;

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
        let calculator = crate::physics::resistance::ResistanceCalculator::new();
        let epsilon = T::default_epsilon();
        let is_blood_like = self.fluid.name().to_ascii_lowercase().contains("blood");
        let t_std =
            T::from_f64(cfd_core::physics::constants::physics::thermo::T_STANDARD)
                .unwrap_or_else(T::one);
        let p_std = T::from_f64(cfd_core::physics::constants::physics::thermo::P_ATM)
            .unwrap_or_else(T::zero);
        let state = self.fluid.properties_at(t_std, p_std)?;
        let density = state.density;
        let viscosity = state.dynamic_viscosity;
        let durst_re_floor = T::from_f64(10.0).unwrap_or_else(T::one);
        let durst_limit = T::from_f64(50.0).unwrap_or_else(T::one);
        let durst_offset = T::from_f64(2.28).unwrap_or_else(T::one);
        let sixty_four = T::from_f64(64.0).unwrap_or_else(T::one);
        let eight = T::from_f64(8.0).unwrap_or_else(T::one);
        let microchannel_limit = T::from_f64(300.0e-6).unwrap_or_else(T::one);
        let fl_reference_viscosity = T::from_f64(3.2).unwrap_or_else(T::one);
        let fl_floor = T::from_f64(0.5).unwrap_or_else(T::zero);
        let tiny = T::from_f64(1.0e-30).unwrap_or(epsilon);

        for edge_idx in edge_indices {
            let flow_rate = self.flow_rates.get(edge_idx.index()).copied().unwrap_or(T::zero());

            // Get geometry and other properties from self.properties
            if let Some(props) = self.properties.get(&edge_idx) {
                if matches!(
                    props.resistance_update_policy,
                    ResistanceUpdatePolicy::FlowInvariant
                ) {
                    continue;
                }
                if let Some(geometry) = &props.geometry {
                    // Map channel::ChannelGeometry to resistance::ChannelGeometry
                    let res_geometry = match &geometry.cross_section {
                        crate::domain::channel::CrossSection::Circular { diameter } => {
                            crate::physics::resistance::ChannelGeometry::Circular {
                                diameter: *diameter,
                                length: geometry.length,
                            }
                        }
                        crate::domain::channel::CrossSection::Rectangular { width, height } => {
                            crate::physics::resistance::ChannelGeometry::Rectangular {
                                width: *width,
                                height: *height,
                                length: geometry.length,
                            }
                        }
                        crate::domain::channel::CrossSection::Elliptical {
                            major_axis,
                            minor_axis,
                        } => crate::physics::resistance::ChannelGeometry::Elliptical {
                            major_axis: *major_axis,
                            minor_axis: *minor_axis,
                            length: geometry.length,
                        },
                        crate::domain::channel::CrossSection::Trapezoidal {
                            top_width,
                            bottom_width,
                            height,
                        } => crate::physics::resistance::ChannelGeometry::Trapezoidal {
                            top_width: *top_width,
                            bottom_width: *bottom_width,
                            height: *height,
                            length: geometry.length,
                        },
                        crate::domain::channel::CrossSection::Custom {
                            area,
                            hydraulic_diameter,
                        } => crate::physics::resistance::ChannelGeometry::Custom {
                            area: *area,
                            hydraulic_diameter: *hydraulic_diameter,
                            length: geometry.length,
                        },
                    };

                    // Prepare flow conditions
                    // Explicitly use from_flow_rate if available, or initialize and set flow rate.
                    // Using standard ambient temperature (293.15 K) as default instead of zero.
                    let mut conditions = crate::physics::resistance::FlowConditions::new(T::zero());
                    // Force velocity to None to ensure calculator derives it from flow rate
                    conditions.velocity = None;
                    conditions.flow_rate = Some(flow_rate.abs());
                    conditions.temperature = t_std;

                    if let Some(d_h) = props.hydraulic_diameter {
                        if d_h > epsilon && props.area > epsilon {
                            let velocity = flow_rate.abs() / props.area;
                            let shear_rate = eight * velocity / d_h;
                            let apparent_viscosity = self.fluid.viscosity_at_shear(
                                shear_rate,
                                conditions.temperature,
                                conditions.pressure,
                            )?;

                            conditions.shear_rate = Some(shear_rate);
                            if apparent_viscosity > epsilon {
                                conditions.reynolds_number =
                                    Some(density * velocity * d_h / apparent_viscosity);
                            }
                        }
                    }

                    // Re-calculate coefficients
                    let (mut r, k) = calculator.calculate_coefficients_auto(
                        &res_geometry,
                        &self.fluid,
                        &conditions,
                    )?;

                    if let Some(d_h) = props.hydraulic_diameter {
                        if d_h > epsilon && props.area > epsilon {
                            let l_over_dh = props.length / d_h;
                            if l_over_dh < durst_limit && viscosity > epsilon {
                                // Preserve the short-channel developing-flow correction during
                                // Picard updates instead of dropping back to fully developed
                                // resistance after the first flow recomputation.
                                let velocity = flow_rate.abs() / props.area;
                                let reynolds = (density * velocity * d_h / viscosity).max(durst_re_floor);
                                let k_entrance =
                                    durst_offset + sixty_four / (reynolds * l_over_dh).max(tiny);
                                let multiplier =
                                    T::one() + k_entrance / (sixty_four * l_over_dh).max(T::one());
                                r *= multiplier;
                            }

                            if is_blood_like && d_h < microchannel_limit {
                                let fl = cfd_core::physics::fluid::blood::FahraeuasLindqvist::new(
                                    d_h, T::from_f64(0.45).unwrap_or_else(T::zero),
                                );
                                if fl.is_significant() {
                                    let reduction =
                                        (fl.relative_viscosity() / fl_reference_viscosity)
                                            .clamp(fl_floor, T::one());
                                    r *= reduction;
                                }
                            }
                        }
                    }

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
            // Using Picard iteration (secant modulus):
            // ΔP ≈ (R + k|Q_k|)·Q.
            // Thus, effective resistance R_eff = R + k|Q_k|.
            let r_eff = edge_data.resistance + edge_data.quad_coeff * q;
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
