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
use petgraph::Direction;
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

/// Edge-property key for discharge/feed hematocrit.
pub const EDGE_PROPERTY_HEMATOCRIT: &str = "feed_hematocrit";
/// Edge-property key for locally propagated daughter/tube hematocrit.
pub const EDGE_PROPERTY_LOCAL_HEMATOCRIT: &str = "local_hematocrit";
/// Edge-property key for plasma dynamic viscosity in Pa·s.
pub const EDGE_PROPERTY_PLASMA_VISCOSITY_PA_S: &str = "plasma_viscosity_pa_s";
/// Edge-property key for a transiently assembled local apparent viscosity in Pa·s.
pub const EDGE_PROPERTY_LOCAL_APPARENT_VISCOSITY_PA_S: &str = "local_apparent_viscosity_pa_s";

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

struct ResistanceUpdateContext<T: RealField + Copy> {
    calculator: crate::physics::resistance::ResistanceCalculator<T>,
    epsilon: T,
    is_blood_like: bool,
    t_std: T,
    density: T,
    viscosity: T,
    durst_limit: T,
    durst_offset: T,
    sixty_four: T,
    eight: T,
    microchannel_limit: T,
    default_hematocrit: T,
    default_plasma_viscosity: T,
    tiny: T,
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
                        flow_rate: self
                            .flow_rates
                            .get(edge_idx.index())
                            .copied()
                            .unwrap_or(T::zero()),
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
        let update_context = self.resistance_update_context()?;

        for &edge_idx in &edge_indices {
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

        self.propagate_blood_hematocrit(&update_context)?;
        for &edge_idx in &edge_indices {
            let flow = self
                .flow_rates
                .get(edge_idx.index())
                .copied()
                .unwrap_or(T::zero());
            self.update_single_edge_resistance(edge_idx, flow, &update_context)?;
        }

        Ok(())
    }

    /// Build `FlowConditions` for an edge given its geometric properties and current flow rate.
    fn build_edge_flow_conditions(
        &self,
        props: &EdgeProperties<T>,
        flow_rate: T,
        t_std: T,
        density: T,
        eight: T,
        epsilon: T,
    ) -> Result<(crate::physics::resistance::FlowConditions<T>, Option<T>)> {
        let mut conditions =
            crate::physics::resistance::FlowConditions::from_flow_rate(flow_rate.abs());
        conditions.temperature = t_std;
        let mut apparent_viscosity = None;
        if let Some(d_h) = props.hydraulic_diameter {
            if d_h > epsilon && props.area > epsilon {
                let velocity = flow_rate.abs() / props.area;
                let shear_rate = eight * velocity / d_h;
                let apparent_viscosity_local = self.fluid.viscosity_at_shear(
                    shear_rate,
                    conditions.temperature,
                    conditions.pressure,
                )?;
                conditions.shear_rate = Some(shear_rate);
                apparent_viscosity = Some(apparent_viscosity_local);
                if apparent_viscosity_local > epsilon {
                    conditions.reynolds_number =
                        Some(density * velocity * d_h / apparent_viscosity_local);
                }
            }
        }
        Ok((conditions, apparent_viscosity))
    }

    fn resistance_update_context(&self) -> Result<ResistanceUpdateContext<T>> {
        let epsilon = T::default_epsilon();
        let t_std = T::from_f64(cfd_core::physics::constants::physics::thermo::T_STANDARD)
            .unwrap_or_else(T::one);
        let p_std = T::from_f64(cfd_core::physics::constants::physics::thermo::P_ATM)
            .unwrap_or_else(T::zero);
        let state = self.fluid.properties_at(t_std, p_std)?;
        Ok(ResistanceUpdateContext {
            calculator: crate::physics::resistance::ResistanceCalculator::new(),
            epsilon,
            is_blood_like: self.fluid.name().to_ascii_lowercase().contains("blood"),
            t_std,
            density: state.density,
            viscosity: state.dynamic_viscosity,
            durst_limit: T::from_f64(50.0).unwrap_or_else(T::one),
            durst_offset: T::from_f64(2.28).unwrap_or_else(T::one),
            sixty_four: T::from_f64(64.0).unwrap_or_else(T::one),
            eight: T::from_f64(8.0).unwrap_or_else(T::one),
            microchannel_limit: T::from_f64(300.0e-6).unwrap_or_else(T::one),
            default_hematocrit: T::from_f64(0.45).unwrap_or_else(T::zero),
            default_plasma_viscosity: state.dynamic_viscosity
                / T::from_f64(3.2).unwrap_or_else(T::one),
            tiny: T::from_f64(1.0e-30).unwrap_or(epsilon),
        })
    }

    fn update_single_edge_resistance(
        &mut self,
        edge_idx: EdgeIndex,
        flow_rate: T,
        context: &ResistanceUpdateContext<T>,
    ) -> Result<()> {
        let (resistance, quad_coeff) = {
            let Some(props) = self.properties.get(&edge_idx) else {
                return Ok(());
            };
            if matches!(
                props.resistance_update_policy,
                ResistanceUpdatePolicy::FlowInvariant
            ) {
                return Ok(());
            }
            let Some(geometry) = &props.geometry else {
                return Ok(());
            };

            let res_geometry = channel_to_res_geometry(&geometry.cross_section, geometry.length);
            let (conditions, apparent_viscosity) = self.build_edge_flow_conditions(
                props,
                flow_rate,
                context.t_std,
                context.density,
                context.eight,
                context.epsilon,
            )?;

            let rectangular_laminar_max = T::from_f64(2300.0).unwrap_or_else(T::one);
            let is_rectangular = matches!(
                &res_geometry,
                crate::physics::resistance::ChannelGeometry::Rectangular { .. }
            );
            let use_rectangular_darcy_surrogate = is_rectangular
                && conditions
                    .reynolds_number
                    .is_some_and(|re| re >= rectangular_laminar_max);

            let (mut resistance, quad_coeff) = if use_rectangular_darcy_surrogate {
                use crate::physics::resistance::models::ResistanceModel;

                let d_h = props.hydraulic_diameter.ok_or_else(|| {
                    Error::InvalidConfiguration(
                        "Rectangular branch requires hydraulic diameter for explicit Darcy-Weisbach surrogate".to_string(),
                    )
                })?;
                if props.area <= context.epsilon {
                    return Err(Error::InvalidConfiguration(
                        "Rectangular branch requires positive area for explicit Darcy-Weisbach surrogate".to_string(),
                    ));
                }
                crate::physics::resistance::models::DarcyWeisbachModel::new(
                    d_h,
                    props.area,
                    props.length,
                    T::zero(),
                )
                .calculate_coefficients(&self.fluid, &conditions)?
            } else {
                context.calculator.calculate_coefficients_auto(
                    &res_geometry,
                    &self.fluid,
                    &conditions,
                )?
            };

            if let Some(d_h) = props.hydraulic_diameter {
                apply_resistance_corrections(
                    &mut resistance,
                    d_h,
                    props.area,
                    props.length,
                    flow_rate,
                    context.density,
                    context.viscosity,
                    context.is_blood_like,
                    context.epsilon,
                    context.durst_limit,
                    context.durst_offset,
                    context.sixty_four,
                    context.tiny,
                    context.microchannel_limit,
                    apparent_viscosity,
                    props
                        .properties
                        .get(EDGE_PROPERTY_LOCAL_HEMATOCRIT)
                        .or_else(|| props.properties.get(EDGE_PROPERTY_HEMATOCRIT))
                        .copied(),
                    props
                        .properties
                        .get(EDGE_PROPERTY_LOCAL_APPARENT_VISCOSITY_PA_S)
                        .copied(),
                    props
                        .properties
                        .get(EDGE_PROPERTY_PLASMA_VISCOSITY_PA_S)
                        .copied(),
                    context.default_hematocrit,
                    context.default_plasma_viscosity,
                );
            }

            (resistance, quad_coeff)
        };

        if let Some(edge) = self.graph.edge_weight_mut(edge_idx) {
            edge.resistance = resistance;
            edge.quad_coeff = quad_coeff;
        }
        Ok(())
    }

    /// Re-calculate resistances and quadratic coefficients for all edges based on current flow rates.
    pub fn update_resistances(&mut self) -> Result<()> {
        let edge_indices: Vec<_> = self.graph.edge_indices().collect();
        let context = self.resistance_update_context()?;
        self.propagate_blood_hematocrit(&context)?;

        for edge_idx in edge_indices {
            let flow_rate = self
                .flow_rates
                .get(edge_idx.index())
                .copied()
                .unwrap_or(T::zero());
            self.update_single_edge_resistance(edge_idx, flow_rate, &context)?;
        }
        Ok(())
    }

    fn propagate_blood_hematocrit(&mut self, context: &ResistanceUpdateContext<T>) -> Result<()> {
        if !context.is_blood_like {
            return Ok(());
        }

        let node_indices: Vec<_> = self.graph.node_indices().collect();
        let tolerance = T::from_f64(1.0e-9).unwrap_or(context.epsilon);
        let max_sweeps = self.graph.node_count().clamp(2, 12);
        let mut node_hematocrit = vec![context.default_hematocrit; self.graph.node_count()];
        let mut next_node_hematocrit = vec![context.default_hematocrit; self.graph.node_count()];
        let mut inflows = Vec::with_capacity(8);
        let mut outflows = Vec::with_capacity(8);
        let mut updates = Vec::with_capacity(self.graph.edge_count());

        for props in self.properties.values_mut() {
            let seed = props
                .properties
                .get(EDGE_PROPERTY_HEMATOCRIT)
                .copied()
                .unwrap_or(context.default_hematocrit);
            props
                .properties
                .entry(EDGE_PROPERTY_LOCAL_HEMATOCRIT.to_string())
                .or_insert(seed);
        }

        for node_idx in &node_indices {
            node_hematocrit[node_idx.index()] = self.node_hematocrit_estimate(
                *node_idx,
                &node_hematocrit,
                context.default_hematocrit,
                &mut inflows,
                &mut outflows,
            );
        }

        for _ in 0..max_sweeps {
            updates.clear();
            let mut max_change = T::zero();

            next_node_hematocrit.copy_from_slice(&node_hematocrit);

            for node_idx in &node_indices {
                let estimate = self.node_hematocrit_estimate(
                    *node_idx,
                    &node_hematocrit,
                    context.default_hematocrit,
                    &mut inflows,
                    &mut outflows,
                );
                max_change = max_change.max((node_hematocrit[node_idx.index()] - estimate).abs());
                next_node_hematocrit[node_idx.index()] = estimate;
            }

            for node_idx in &node_indices {
                self.node_edge_fluxes_into(*node_idx, true, &mut inflows);
                self.node_edge_fluxes_into(*node_idx, false, &mut outflows);
                if outflows.is_empty() {
                    continue;
                }
                let node_hct = next_node_hematocrit[node_idx.index()];

                if inflows.len() == 1 && outflows.len() == 2 {
                    if let Some((edge_a, hematocrit_a, edge_b, hematocrit_b)) = self
                        .bifurcation_hematocrit_split(
                            inflows[0].0,
                            outflows[0],
                            outflows[1],
                            node_hct,
                            context.epsilon,
                        )
                    {
                        max_change = max_change.max(
                            (self.edge_hematocrit(edge_a, context.default_hematocrit)
                                - hematocrit_a)
                                .abs(),
                        );
                        max_change = max_change.max(
                            (self.edge_hematocrit(edge_b, context.default_hematocrit)
                                - hematocrit_b)
                                .abs(),
                        );
                        updates.push((edge_a, hematocrit_a));
                        updates.push((edge_b, hematocrit_b));
                        continue;
                    }
                }

                for &(edge_idx, _) in &outflows {
                    max_change = max_change.max(
                        (self.edge_hematocrit(edge_idx, context.default_hematocrit) - node_hct)
                            .abs(),
                    );
                    updates.push((edge_idx, node_hct));
                }
            }

            for (edge_idx, hematocrit) in updates.drain(..) {
                if let Some(props) = self.properties.get_mut(&edge_idx) {
                    props.properties.insert(
                        EDGE_PROPERTY_LOCAL_HEMATOCRIT.to_string(),
                        hematocrit.clamp(T::zero(), T::one()),
                    );
                }
            }
            std::mem::swap(&mut node_hematocrit, &mut next_node_hematocrit);

            if max_change <= tolerance {
                break;
            }
        }

        Ok(())
    }

    fn node_hematocrit_estimate(
        &self,
        node_idx: NodeIndex,
        node_hematocrit: &[T],
        default_hematocrit: T,
        inflows: &mut Vec<(EdgeIndex, T)>,
        outflows: &mut Vec<(EdgeIndex, T)>,
    ) -> T {
        self.node_edge_fluxes_into(node_idx, true, inflows);
        let total_in = inflows.iter().fold(T::zero(), |acc, (_, q)| acc + *q);
        if total_in > T::default_epsilon() {
            return inflows.iter().fold(T::zero(), |acc, (edge_idx, q)| {
                acc + *q * self.edge_hematocrit(*edge_idx, default_hematocrit)
            }) / total_in;
        }

        self.node_edge_fluxes_into(node_idx, false, outflows);
        let total_out = outflows.iter().fold(T::zero(), |acc, (_, q)| acc + *q);
        if total_out > T::default_epsilon() {
            let seeded =
                outflows
                    .iter()
                    .fold((T::zero(), T::zero()), |(acc_h, acc_q), (edge_idx, q)| {
                        let seeded_h = self
                            .properties
                            .get(edge_idx)
                            .and_then(|props| {
                                props.properties.get(EDGE_PROPERTY_HEMATOCRIT).copied()
                            })
                            .unwrap_or_else(|| self.edge_hematocrit(*edge_idx, default_hematocrit));
                        (acc_h + *q * seeded_h, acc_q + *q)
                    });
            if seeded.1 > T::default_epsilon() {
                return seeded.0 / seeded.1;
            }
        }

        node_hematocrit
            .get(node_idx.index())
            .copied()
            .unwrap_or(default_hematocrit)
    }

    fn node_edge_fluxes_into(
        &self,
        node_idx: NodeIndex,
        inflow: bool,
        fluxes: &mut Vec<(EdgeIndex, T)>,
    ) {
        fluxes.clear();

        for edge_ref in self.graph.edges_directed(node_idx, Direction::Incoming) {
            let edge_idx = edge_ref.id();
            let flow = self
                .flow_rates
                .get(edge_idx.index())
                .copied()
                .unwrap_or_else(T::zero);
            if inflow && flow > T::zero() {
                fluxes.push((edge_idx, flow));
            } else if !inflow && flow < T::zero() {
                fluxes.push((edge_idx, -flow));
            }
        }

        for edge_ref in self.graph.edges_directed(node_idx, Direction::Outgoing) {
            let edge_idx = edge_ref.id();
            let flow = self
                .flow_rates
                .get(edge_idx.index())
                .copied()
                .unwrap_or_else(T::zero);
            if inflow && flow < T::zero() {
                fluxes.push((edge_idx, -flow));
            } else if !inflow && flow > T::zero() {
                fluxes.push((edge_idx, flow));
            }
        }
    }

    fn edge_hematocrit(&self, edge_idx: EdgeIndex, default_hematocrit: T) -> T {
        self.properties
            .get(&edge_idx)
            .and_then(|props| {
                props
                    .properties
                    .get(EDGE_PROPERTY_LOCAL_HEMATOCRIT)
                    .or_else(|| props.properties.get(EDGE_PROPERTY_HEMATOCRIT))
                    .copied()
            })
            .unwrap_or(default_hematocrit)
    }

    fn bifurcation_hematocrit_split(
        &self,
        parent_edge: EdgeIndex,
        daughter_a: (EdgeIndex, T),
        daughter_b: (EdgeIndex, T),
        feed_hematocrit: T,
        epsilon: T,
    ) -> Option<(EdgeIndex, T, EdgeIndex, T)> {
        let parent_diameter = self
            .properties
            .get(&parent_edge)?
            .hydraulic_diameter
            .and_then(|d| nalgebra::try_convert::<T, f64>(d))?;
        let daughter_a_diameter = self
            .properties
            .get(&daughter_a.0)?
            .hydraulic_diameter
            .and_then(|d| nalgebra::try_convert::<T, f64>(d))?;
        let daughter_b_diameter = self
            .properties
            .get(&daughter_b.0)?
            .hydraulic_diameter
            .and_then(|d| nalgebra::try_convert::<T, f64>(d))?;
        let q_a = nalgebra::try_convert::<T, f64>(daughter_a.1)?;
        let q_b = nalgebra::try_convert::<T, f64>(daughter_b.1)?;
        let h_feed = nalgebra::try_convert::<T, f64>(feed_hematocrit)?;

        let total_q = q_a + q_b;
        if total_q <= nalgebra::try_convert::<T, f64>(epsilon)? {
            return None;
        }

        let result_a = crate::physics::cell_separation::plasma_skimming::pries_phase_separation(
            h_feed.clamp(0.0, 1.0),
            (q_a / total_q).clamp(0.0, 1.0),
            daughter_a_diameter * 1.0e6,
            daughter_b_diameter * 1.0e6,
            parent_diameter * 1.0e6,
        )
        .ok()?;

        let flow_fraction_b: f64 = (q_b / total_q).clamp(0.0, 1.0);
        let hematocrit_b = if flow_fraction_b > 1.0e-15_f64 {
            ((1.0 - result_a.cell_fraction) * h_feed / flow_fraction_b).clamp(0.0, 1.0)
        } else {
            0.0
        };

        Some((
            daughter_a.0,
            T::from_f64(result_a.daughter_hematocrit)?,
            daughter_b.0,
            T::from_f64(hematocrit_b)?,
        ))
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

/// Map a `CrossSection` + length to the resistance-crate `ChannelGeometry`.
fn channel_to_res_geometry<T: nalgebra::RealField + Copy>(
    cross_section: &crate::domain::channel::CrossSection<T>,
    length: T,
) -> crate::physics::resistance::ChannelGeometry<T> {
    match cross_section {
        crate::domain::channel::CrossSection::Circular { diameter } => {
            crate::physics::resistance::ChannelGeometry::Circular {
                diameter: *diameter,
                length,
            }
        }
        crate::domain::channel::CrossSection::Rectangular { width, height } => {
            crate::physics::resistance::ChannelGeometry::Rectangular {
                width: *width,
                height: *height,
                length,
            }
        }
        crate::domain::channel::CrossSection::Elliptical {
            major_axis,
            minor_axis,
        } => crate::physics::resistance::ChannelGeometry::Elliptical {
            major_axis: *major_axis,
            minor_axis: *minor_axis,
            length,
        },
        crate::domain::channel::CrossSection::Trapezoidal {
            top_width,
            bottom_width,
            height,
        } => crate::physics::resistance::ChannelGeometry::Trapezoidal {
            top_width: *top_width,
            bottom_width: *bottom_width,
            height: *height,
            length,
        },
        crate::domain::channel::CrossSection::Custom {
            area,
            hydraulic_diameter,
        } => crate::physics::resistance::ChannelGeometry::Custom {
            area: *area,
            hydraulic_diameter: *hydraulic_diameter,
            length,
        },
    }
}

/// Apply short-channel (Durst) and Fåhræus–Lindqvist blood viscosity corrections to `r`.
///
/// The Durst branch uses the published entrance-loss relation directly; no
/// Reynolds-number floor is introduced in the physical model. A tiny
/// denominator guard remains only to avoid division by zero at exact zero flow.
#[allow(clippy::too_many_arguments)]
fn apply_resistance_corrections<T: nalgebra::RealField + Copy>(
    r: &mut T,
    d_h: T,
    area: T,
    length: T,
    flow_rate: T,
    density: T,
    viscosity: T,
    is_blood_like: bool,
    epsilon: T,
    durst_limit: T,
    durst_offset: T,
    sixty_four: T,
    tiny: T,
    microchannel_limit: T,
    apparent_viscosity: Option<T>,
    hematocrit: Option<T>,
    apparent_viscosity_override: Option<T>,
    plasma_viscosity: Option<T>,
    default_hematocrit: T,
    default_plasma_viscosity: T,
) {
    if d_h > epsilon && area > epsilon {
        let l_over_dh = length / d_h;
        if l_over_dh < durst_limit && viscosity > epsilon {
            let velocity = flow_rate.abs() / area;
            let reynolds = density * velocity * d_h / viscosity;
            let k_entrance = durst_offset + sixty_four / (reynolds * l_over_dh).max(tiny);
            let multiplier = T::one() + k_entrance / (sixty_four * l_over_dh).max(T::one());
            *r *= multiplier;
        }

        if is_blood_like && d_h < microchannel_limit {
            let local_mu = apparent_viscosity.unwrap_or(viscosity);
            if local_mu > epsilon {
                let target_mu = apparent_viscosity_override.or_else(|| {
                    blood_microchannel_apparent_viscosity(
                        d_h,
                        flow_rate,
                        area,
                        hematocrit.unwrap_or(default_hematocrit),
                        plasma_viscosity.unwrap_or(default_plasma_viscosity),
                    )
                });
                if let Some(target_mu) = target_mu {
                    let ratio = target_mu / local_mu;
                    if ratio.is_finite() && ratio > epsilon {
                        *r *= ratio;
                    }
                }
            }
        }
    }
}

/// Estimate microchannel blood apparent viscosity from local shear, diameter,
/// hematocrit, and plasma viscosity using Secomb diameter correction together
/// with Quemada low-shear aggregation.
pub fn blood_microchannel_apparent_viscosity<T: nalgebra::RealField + Copy + FromPrimitive>(
    d_h: T,
    flow_rate: T,
    area: T,
    hematocrit: T,
    plasma_viscosity: T,
) -> Option<T> {
    let d_h_m = nalgebra::try_convert::<T, f64>(d_h)?;
    let flow_rate_m3_s = nalgebra::try_convert::<T, f64>(flow_rate.abs())?;
    let area_m2 = nalgebra::try_convert::<T, f64>(area)?;
    let hematocrit = nalgebra::try_convert::<T, f64>(hematocrit)?;
    let plasma_viscosity = nalgebra::try_convert::<T, f64>(plasma_viscosity)?;
    if d_h_m <= 0.0 || area_m2 <= 0.0 || plasma_viscosity <= 0.0 {
        return None;
    }

    let velocity = flow_rate_m3_s / area_m2;
    let shear_rate = (8.0 * velocity / d_h_m).abs();
    let diameter_um = d_h_m * 1.0e6;
    let secomb = crate::physics::cell_separation::fahraeus_lindqvist::secomb_network_viscosity(
        diameter_um,
        hematocrit.clamp(0.0, 0.95),
        plasma_viscosity,
    );
    let quemada = crate::physics::cell_separation::rouleaux_aggregation::checked_quemada_viscosity(
        shear_rate,
        hematocrit.clamp(0.0, 0.95),
        plasma_viscosity,
    )
    .unwrap_or(secomb);

    let target = secomb.max(quemada);
    T::from_f64(target)
}

#[cfg(test)]
mod tests {
    use super::*;
    use approx::assert_relative_eq;
    use crate::physics::cell_separation::fahraeus_lindqvist::secomb_network_viscosity;

    #[test]
    fn blood_microchannel_apparent_viscosity_falls_back_to_secomb_when_quemada_invalid() {
        let d_h = 50.0e-6;
        let area = std::f64::consts::PI * (d_h * 0.5) * (d_h * 0.5);
        let hematocrit = 0.95;
        let plasma_viscosity = 1.2e-3;

        let result = blood_microchannel_apparent_viscosity(
            d_h,
            0.0,
            area,
            hematocrit,
            plasma_viscosity,
        )
        .expect("microchannel apparent viscosity");

        let expected = secomb_network_viscosity(d_h * 1.0e6, hematocrit, plasma_viscosity);
        assert!((result - expected).abs() < 1e-15 * expected.max(1.0));
    }

    #[test]
    fn short_channel_correction_uses_exact_low_reynolds_durst_multiplier() {
        let mut resistance = 1.0_f64;
        let d_h: f64 = 100.0e-6;
        let area = std::f64::consts::PI * (d_h * 0.5_f64).powi(2);
        let length = 5.0_f64 * d_h;
        let flow_rate: f64 = 1.0e-13;
        let density: f64 = 1_000.0;
        let viscosity: f64 = 1.0e-3;
        let epsilon = f64::EPSILON;
        let tiny = 1.0e-30_f64;
        let l_over_dh = length / d_h;
        let reynolds = density * (flow_rate.abs() / area) * d_h / viscosity;
        let expected_multiplier =
            1.0 + (2.28 + 64.0 / (reynolds * l_over_dh).max(tiny)) / (64.0 * l_over_dh).max(1.0);

        apply_resistance_corrections(
            &mut resistance,
            d_h,
            area,
            length,
            flow_rate,
            density,
            viscosity,
            false,
            epsilon,
            50.0,
            2.28,
            64.0,
            tiny,
            300.0e-6,
            None,
            None,
            None,
            None,
            0.45,
            viscosity / 3.2,
        );

        assert_relative_eq!(resistance, expected_multiplier, max_relative = 1e-12);
    }
}
