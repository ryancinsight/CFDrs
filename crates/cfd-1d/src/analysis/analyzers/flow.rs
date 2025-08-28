//! Flow analysis for network components

use super::traits::NetworkAnalyzer;
use crate::analysis::FlowAnalysis;
use crate::channel::FlowRegime;
use crate::network::Network;
use cfd_core::Result;
use nalgebra::RealField;
use num_traits::{Float, FromPrimitive};
use std::iter::Sum;

/// Flow analyzer for network components
pub struct FlowAnalyzer<T: RealField + Copy> {
    _phantom: std::marker::PhantomData<T>,
}

impl<T: RealField + Copy> FlowAnalyzer<T> {
    /// Create new flow analyzer
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
            if let Some(flow_rate) = edge.flow_rate {
                analysis.add_component_flow(edge.id.clone(), flow_rate);

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

    fn name(&self) -> &str {
        "FlowAnalyzer"
    }
}

impl<T: RealField + Copy + FromPrimitive + Float> FlowAnalyzer<T> {
    fn determine_flow_regime(
        &self,
        network: &Network<T>,
        edge: &crate::network::EdgeWithProperties<T>,
        flow_rate: T,
    ) -> FlowRegime {
        let fluid = network.fluid();
        let properties = edge.properties;

        // Calculate Reynolds number
        let area = properties.area;
        let hydraulic_diameter = properties.hydraulic_diameter.unwrap_or_else(|| {
            T::from_f64(4.0).unwrap_or_else(T::one) * area
                / (T::from_f64(std::f64::consts::PI).unwrap_or_else(T::one) * Float::sqrt(area))
        });

        let velocity = flow_rate / area;
        let reynolds = fluid.reynolds_number(velocity, hydraulic_diameter);

        FlowRegime::from_reynolds_number(reynolds)
    }

    fn calculate_total_flow(&self, network: &Network<T>) -> T {
        // Sum outflow from all outlet nodes
        let mut total = T::zero();
        for (idx, node) in network.nodes().enumerate() {
            if matches!(node.node_type, crate::network::NodeType::Outlet) {
                if let Some(flow) = network.flow_rates().get(idx) {
                    total = total + flow.abs();
                }
            }
        }
        total
    }
}
