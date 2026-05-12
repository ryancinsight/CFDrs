//! Flow analysis for network components

use super::traits::NetworkAnalyzer;
use crate::domain::channel::FlowRegime;
use crate::domain::network::Network;
use crate::solver::analysis::FlowAnalysis;
use cfd_core::error::Result;
use nalgebra::RealField;
use num_traits::{Float, FromPrimitive};
use petgraph::visit::EdgeRef;
use std::iter::Sum;

/// Flow analyzer for network components
pub struct FlowAnalyzer<T: RealField + Copy> {
    _phantom: std::marker::PhantomData<T>,
}

impl<T: RealField + Copy> Default for FlowAnalyzer<T> {
    fn default() -> Self {
        Self::new()
    }
}

impl<T: RealField + Copy> FlowAnalyzer<T> {
    /// Create new flow analyzer
    #[must_use]
    pub fn new() -> Self {
        Self {
            _phantom: std::marker::PhantomData,
        }
    }
}

impl<T: RealField + Copy + FromPrimitive + Float + Sum> NetworkAnalyzer<T> for FlowAnalyzer<T> {
    type Result = FlowAnalysis<T>;

    fn analyze(&mut self, network: &Network<T>) -> Result<FlowAnalysis<T>> {
        let mut analysis = FlowAnalysis::new();

        // Analyze flow in each edge
        for edge in network.edges_with_properties() {
            let flow_rate = edge.flow_rate;
            if flow_rate != T::zero() {
                analysis.add_component_flow(edge.id.clone(), flow_rate);

                let area = edge.properties.area;
                let velocity = flow_rate / area;
                let velocity_mag = Float::abs(velocity);
                analysis.add_velocity(edge.id.clone(), velocity_mag);

                let hydraulic_diameter = edge.properties.hydraulic_diameter.unwrap_or_else(|| {
                    (T::one() + T::one() + T::one() + T::one()) * area
                        / (T::pi() * Float::sqrt(area))
                });

                let reynolds = network.fluid().density * velocity_mag * hydraulic_diameter
                    / network.fluid().viscosity;
                analysis.add_reynolds_number(edge.id.clone(), reynolds);

                if hydraulic_diameter > T::zero() {
                    let eight =
                        T::from_f64(8.0).expect("Mathematical constant conversion compromised");
                    let shear_rate = eight * velocity_mag / hydraulic_diameter;
                    let shear_stress = network.fluid().viscosity * shear_rate;
                    analysis.add_wall_shear_rate(edge.id.clone(), shear_rate);
                    analysis.add_wall_shear_stress(edge.id.clone(), shear_stress);
                }

                // Determine flow regime
                let regime = self.determine_flow_regime(network, &edge, flow_rate);
                analysis.add_flow_regime(edge.id.clone(), regime);
            }
        }

        // Calculate total system flow
        let total_flow = self.calculate_total_flow(network);
        analysis.set_total_flow(total_flow);

        Ok(analysis)
    }

    fn name(&self) -> &'static str {
        "FlowAnalyzer"
    }
}

impl<T: RealField + Copy + FromPrimitive + Float> FlowAnalyzer<T> {
    fn determine_flow_regime(
        &self,
        network: &Network<T>,
        edge: &crate::domain::network::EdgeWithProperties<T>,
        flow_rate: T,
    ) -> FlowRegime {
        let fluid = network.fluid();
        let properties = edge.properties;

        // Calculate Reynolds number
        let area = properties.area;
        let hydraulic_diameter = properties.hydraulic_diameter.unwrap_or_else(|| {
            (T::one() + T::one() + T::one() + T::one()) * area / (T::pi() * Float::sqrt(area))
        });

        let velocity = flow_rate / area;
        let reynolds = fluid.density * Float::abs(velocity) * hydraulic_diameter / fluid.viscosity;

        FlowRegime::from_reynolds_number(reynolds)
    }

    fn calculate_total_flow(&self, network: &Network<T>) -> T {
        // Sum outflow from all outlet nodes
        use petgraph::graph::NodeIndex;
        let mut total = T::zero();
        for (idx, node) in network.nodes().enumerate() {
            if matches!(node.node_type, crate::domain::network::NodeType::Outlet) {
                // Sum flow rates of edges connected to this outlet node
                let node_idx = NodeIndex::new(idx);
                for edge_ref in network.graph.edges(node_idx) {
                    let edge_idx = edge_ref.id();
                    if let Some(&flow) = network.flow_rates().get(edge_idx.index()) {
                        total += Float::abs(flow);
                    }
                }
            }
        }
        total
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::domain::network::{
        ComponentType, EdgeProperties, Network, NetworkBuilder, ResistanceUpdatePolicy,
    };
    use cfd_core::physics::fluid::database::water_20c;
    use std::collections::HashMap;

    fn reverse_flow_network(flow_rate: f64) -> Network<f64> {
        let mut builder = NetworkBuilder::new();
        let inlet = builder.add_inlet("inlet".to_string());
        let outlet = builder.add_outlet("outlet".to_string());
        let edge = builder.connect_with_pipe(inlet, outlet, "pipe".to_string());
        let graph = builder
            .build()
            .expect("two-node network with inlet and outlet must validate");
        let mut network = Network::new(graph, water_20c::<f64>().expect("water properties exist"));
        network.add_edge_properties(
            edge,
            EdgeProperties {
                id: "pipe".to_string(),
                component_type: ComponentType::Pipe,
                length: 0.1,
                area: 1.0e-4,
                hydraulic_diameter: Some(0.01),
                resistance: 1.0,
                geometry: None,
                resistance_update_policy: ResistanceUpdatePolicy::FlowInvariant,
                properties: HashMap::new(),
            },
        );
        network.set_flow_rate(edge, flow_rate);
        network
    }

    #[test]
    fn reverse_flow_regime_uses_reynolds_magnitude() -> Result<()> {
        let mut analyzer = FlowAnalyzer::<f64>::new();
        let analysis = analyzer.analyze(&reverse_flow_network(-3.0e-5))?;

        let reynolds = *analysis
            .reynolds_numbers
            .get("pipe")
            .expect("pipe Reynolds number must be recorded");
        let regime = analysis
            .flow_regimes
            .get("pipe")
            .expect("pipe flow regime must be recorded");

        assert!(
            reynolds > 2300.0 && reynolds < 4000.0,
            "test network must produce transitional Reynolds magnitude, got {reynolds}"
        );
        assert_eq!(regime, &FlowRegime::Transitional);
        Ok(())
    }

    #[test]
    fn forward_and_reverse_flow_share_scalar_diagnostics() -> Result<()> {
        let mut analyzer = FlowAnalyzer::<f64>::new();
        let forward = analyzer.analyze(&reverse_flow_network(3.0e-5))?;
        let reverse = analyzer.analyze(&reverse_flow_network(-3.0e-5))?;

        assert_eq!(
            forward.flow_regimes.get("pipe"),
            reverse.flow_regimes.get("pipe")
        );
        assert_eq!(
            forward.reynolds_numbers.get("pipe"),
            reverse.reynolds_numbers.get("pipe")
        );
        assert_eq!(
            forward.wall_shear_rates.get("pipe"),
            reverse.wall_shear_rates.get("pipe")
        );
        Ok(())
    }
}
