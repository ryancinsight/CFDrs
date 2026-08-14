//! Flow analysis for network components

use super::traits::NetworkAnalyzer;
use crate::domain::channel::FlowRegime;
use crate::domain::network::Network;
use crate::solver::analysis::FlowAnalysis;
use aequitas::systems::si::quantities::{
    Dimensionless, Pressure, ReciprocalTime, Velocity, VolumetricFlowRate,
};
use cfd_core::conversion::{SafeFromF64, SafeFromUsize};
use cfd_core::error::Result;
use cfd_core::CfdScalar;
use eunomia::NumericElement;
use petgraph::visit::EdgeRef;
use std::iter::Sum;

/// Flow analyzer for network components
pub struct FlowAnalyzer<T: CfdScalar + Copy> {
    _phantom: std::marker::PhantomData<T>,
}

impl<T: CfdScalar + Copy> Default for FlowAnalyzer<T> {
    fn default() -> Self {
        Self::new()
    }
}

impl<T: CfdScalar + Copy> FlowAnalyzer<T> {
    /// Create new flow analyzer
    #[must_use]
    pub fn new() -> Self {
        Self {
            _phantom: std::marker::PhantomData,
        }
    }
}

impl<T: CfdScalar + Copy + SafeFromF64 + SafeFromUsize + Sum> NetworkAnalyzer<T>
    for FlowAnalyzer<T>
{
    type Result = FlowAnalysis<T>;

    fn analyze(&mut self, network: &Network<T>) -> Result<FlowAnalysis<T>> {
        let mut analysis = FlowAnalysis::new();

        // Analyze flow in each edge
        for edge in network.edges_with_properties() {
            let flow_rate = edge.flow_rate.into_base();
            if flow_rate != T::ZERO {
                analysis
                    .add_component_flow(edge.id.clone(), VolumetricFlowRate::from_base(flow_rate));

                let area = edge.properties.area.into_base();
                let velocity = flow_rate / area;
                let velocity_mag = <T as NumericElement>::abs(velocity);
                analysis.add_velocity(edge.id.clone(), Velocity::from_base(velocity_mag));

                let hydraulic_diameter = edge.properties.hydraulic_diameter.map_or_else(
                    || {
                        (T::ONE + T::ONE + T::ONE + T::ONE) * area
                            / (T::pi() * <T as NumericElement>::sqrt(area))
                    },
                    |diameter| diameter.into_base(),
                );

                let reynolds =
                    network.fluid().density.into_base() * velocity_mag * hydraulic_diameter
                        / network.fluid().viscosity.into_base();
                analysis.add_reynolds_number(edge.id.clone(), Dimensionless::from_base(reynolds));

                if hydraulic_diameter > T::ZERO {
                    let eight = T::from_f64_or_one(8.0);
                    let shear_rate = eight * velocity_mag / hydraulic_diameter;
                    let shear_stress = network.fluid().viscosity.into_base() * shear_rate;
                    analysis.add_wall_shear_rate(
                        edge.id.clone(),
                        ReciprocalTime::from_base(shear_rate),
                    );
                    analysis
                        .add_wall_shear_stress(edge.id.clone(), Pressure::from_base(shear_stress));
                }

                // Determine flow regime
                let regime = self.determine_flow_regime(network, &edge, flow_rate);
                analysis.add_flow_regime(edge.id.clone(), regime);
            }
        }

        // Calculate total system flow
        let total_flow = self.calculate_total_flow(network);
        analysis.set_total_flow(VolumetricFlowRate::from_base(total_flow));

        Ok(analysis)
    }

    fn name(&self) -> &'static str {
        "FlowAnalyzer"
    }
}

impl<T: CfdScalar + Copy + SafeFromF64> FlowAnalyzer<T> {
    fn determine_flow_regime(
        &self,
        network: &Network<T>,
        edge: &crate::domain::network::EdgeWithProperties<T>,
        flow_rate: T,
    ) -> FlowRegime {
        let fluid = network.fluid();
        let properties = edge.properties;

        // Calculate Reynolds number
        let area = properties.area.into_base();
        let hydraulic_diameter = properties.hydraulic_diameter.map_or_else(
            || {
                (T::ONE + T::ONE + T::ONE + T::ONE) * area
                    / (T::pi() * <T as NumericElement>::sqrt(area))
            },
            |diameter| diameter.into_base(),
        );

        let velocity = flow_rate / area;
        let reynolds =
            fluid.density.into_base() * <T as NumericElement>::abs(velocity) * hydraulic_diameter
                / fluid.viscosity.into_base();

        FlowRegime::from_reynolds_number(reynolds)
    }

    fn calculate_total_flow(&self, network: &Network<T>) -> T {
        // Sum outflow from all outlet nodes
        use petgraph::graph::NodeIndex;
        let mut total = T::ZERO;
        for (idx, node) in network.nodes().enumerate() {
            if matches!(node.node_type, crate::domain::network::NodeType::Outlet) {
                // Sum flow rates of edges connected to this outlet node
                let node_idx = NodeIndex::new(idx);
                for edge_ref in network.graph.edges(node_idx) {
                    let edge_idx = edge_ref.id();
                    if let Some(&flow) = network.flow_rates().get(edge_idx.index()) {
                        total += <T as NumericElement>::abs(flow.into_base());
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
    use aequitas::systems::si::quantities::HydraulicResistance;
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
                length: aequitas::systems::si::quantities::Length::from_base(0.1),
                area: aequitas::systems::si::quantities::Area::from_base(1.0e-4),
                hydraulic_diameter: Some(aequitas::systems::si::quantities::Length::from_base(
                    0.01,
                )),
                resistance: HydraulicResistance::from_base(1.0),
                geometry: None,
                resistance_update_policy: ResistanceUpdatePolicy::FlowInvariant,
                properties: HashMap::new(),
            },
        );
        network.set_flow_rate(edge, VolumetricFlowRate::from_base(flow_rate));
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
            reynolds.into_base() > 2300.0 && reynolds.into_base() < 4000.0,
            "test network must produce transitional Reynolds magnitude, got {:?}",
            reynolds.into_base()
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
