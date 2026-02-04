//! Resistance analysis for network components

use super::traits::NetworkAnalyzer;
use crate::analysis::error::ResistanceCalculationError;
use crate::analysis::ResistanceAnalysis;
use crate::network::{Network, NetworkGraphExt};
use cfd_core::error::Result;
use cfd_core::physics::constants::physics::thermo::{P_ATM, T_STANDARD};
use nalgebra::RealField;
use num_traits::{Float, FromPrimitive};
use petgraph::algo::all_simple_paths;
use petgraph::visit::EdgeRef;
use std::iter::Sum;

/// Resistance analyzer for network components
pub struct ResistanceAnalyzer<T: RealField + Copy> {
    _phantom: std::marker::PhantomData<T>,
}

impl<T: RealField + Copy> Default for ResistanceAnalyzer<T> {
    fn default() -> Self {
        Self::new()
    }
}

impl<T: RealField + Copy> ResistanceAnalyzer<T> {
    /// Create new resistance analyzer
    #[must_use]
    pub fn new() -> Self {
        Self {
            _phantom: std::marker::PhantomData,
        }
    }
}

impl<T: RealField + Copy + FromPrimitive + Float + Sum> NetworkAnalyzer<T>
    for ResistanceAnalyzer<T>
{
    type Result = ResistanceAnalysis<T>;

    fn analyze(&mut self, network: &Network<T>) -> Result<ResistanceAnalysis<T>> {
        let mut analysis = ResistanceAnalysis::new();
        self.populate_edge_resistances(network, &mut analysis)?;

        for path in self.critical_paths(network, &analysis) {
            analysis.add_critical_path(path);
        }

        Ok(analysis)
    }

    fn name(&self) -> &'static str {
        "ResistanceAnalyzer"
    }
}

impl<T: RealField + Copy + FromPrimitive + Float + Sum> ResistanceAnalyzer<T> {
    fn populate_edge_resistances(
        &self,
        network: &Network<T>,
        analysis: &mut ResistanceAnalysis<T>,
    ) -> Result<()> {
        let fluid = network.fluid();

        for edge in network.edges_with_properties() {
            let flow_rate = if edge.flow_rate == T::zero() {
                None
            } else {
                Some(edge.flow_rate)
            };

            let resistance = self
                .calculate_resistance(edge.properties, fluid, flow_rate)
                .map_err(|e| {
                    cfd_core::error::Error::InvalidInput(format!(
                        "Failed to analyze resistance for edge '{}': {}",
                        edge.id, e
                    ))
                })?;

            analysis.add_resistance(edge.id.clone(), resistance);

            let component_type = edge.properties.component_type;
            analysis.add_resistance_by_type(component_type.as_str().to_string(), resistance);
        }

        Ok(())
    }

    fn critical_paths(
        &self,
        network: &Network<T>,
        analysis: &ResistanceAnalysis<T>,
    ) -> Vec<Vec<String>> {
        let inlet_nodes = network.graph.inlet_nodes();
        let outlet_nodes = network.graph.outlet_nodes();
        if inlet_nodes.is_empty() || outlet_nodes.is_empty() {
            return Vec::new();
        }

        let mut best_resistance: Option<T> = None;
        let mut best_paths: Vec<Vec<String>> = Vec::new();
        let epsilon = T::default_epsilon();

        for inlet in inlet_nodes {
            for outlet in &outlet_nodes {
                for node_path in
                    all_simple_paths::<Vec<_>, _>(&network.graph, inlet, *outlet, 0, None)
                {
                    let mut edge_candidates: Vec<Vec<_>> = Vec::new();
                    let mut valid = true;

                    for window in node_path.windows(2) {
                        let from = window[0];
                        let to = window[1];
                        let edges: Vec<_> = network
                            .graph
                            .edges_connecting(from, to)
                            .map(|edge| edge.id())
                            .collect();
                        if edges.is_empty() {
                            valid = false;
                            break;
                        }
                        edge_candidates.push(edges);
                    }

                    if !valid {
                        continue;
                    }

                    let mut edge_paths: Vec<Vec<_>> = vec![Vec::new()];
                    for edges in edge_candidates {
                        let mut next_paths = Vec::new();
                        for path in &edge_paths {
                            for edge_idx in &edges {
                                let mut extended = path.clone();
                                extended.push(*edge_idx);
                                next_paths.push(extended);
                            }
                        }
                        edge_paths = next_paths;
                    }

                    for edge_path in edge_paths {
                        let mut resistance_sum = T::zero();
                        let mut edge_ids = Vec::with_capacity(edge_path.len());
                        let mut edge_valid = true;

                        for edge_idx in edge_path {
                            if let Some(edge) = network.graph.edge_weight(edge_idx) {
                                let resistance = analysis
                                    .resistances
                                    .get(&edge.id)
                                    .copied()
                                    .unwrap_or(edge.resistance);
                                resistance_sum += resistance;
                                edge_ids.push(edge.id.clone());
                            } else {
                                edge_valid = false;
                                break;
                            }
                        }

                        if !edge_valid || edge_ids.is_empty() {
                            continue;
                        }

                        match best_resistance {
                            None => {
                                best_resistance = Some(resistance_sum);
                                best_paths = vec![edge_ids];
                            }
                            Some(best) => {
                                let diff = num_traits::Float::abs(resistance_sum - best);
                                let scale = num_traits::Float::abs(best) + T::one();
                                if resistance_sum > best {
                                    best_resistance = Some(resistance_sum);
                                    best_paths = vec![edge_ids];
                                } else if diff <= epsilon * scale {
                                    best_paths.push(edge_ids);
                                }
                            }
                        }
                    }
                }
            }
        }

        best_paths
    }
}

impl<T: RealField + Copy + FromPrimitive + Float> ResistanceAnalyzer<T> {
    fn calculate_resistance(
        &self,
        properties: &crate::network::EdgeProperties<T>,
        fluid: &cfd_core::physics::fluid::Fluid<T>,
        flow_rate: Option<T>,
    ) -> std::result::Result<T, ResistanceCalculationError> {
        use crate::resistance::{FlowConditions, HagenPoiseuilleModel, ResistanceModel};

        // Require hydraulic diameter - no silent fallbacks
        let hydraulic_diameter = properties
            .hydraulic_diameter
            .ok_or(ResistanceCalculationError::MissingHydraulicDiameter)?;

        // Create resistance model with validated parameters
        let model = HagenPoiseuilleModel::new(hydraulic_diameter, properties.length);

        const REF_TEMPERATURE_KEY: &str = "reference_temperature";
        const REF_PRESSURE_KEY: &str = "reference_pressure";

        let temperature = properties
            .properties
            .get(REF_TEMPERATURE_KEY)
            .copied()
            .unwrap_or_else(|| T::from_f64(T_STANDARD).unwrap_or_else(T::one));
        let pressure = properties
            .properties
            .get(REF_PRESSURE_KEY)
            .copied()
            .unwrap_or_else(|| T::from_f64(P_ATM).unwrap_or_else(T::one));

        let conditions = FlowConditions {
            reynolds_number: flow_rate.map(|q| {
                let velocity = q / properties.area;
                fluid.density * velocity * hydraulic_diameter / fluid.viscosity
            }),
            velocity: flow_rate.map(|q| q / properties.area),
            flow_rate,
            temperature,
            pressure,
        };

        // Calculate resistance and propagate any errors
        model
            .calculate_resistance(fluid, &conditions)
            .map_err(|e| ResistanceCalculationError::ModelError(e.to_string()))
    }
}
