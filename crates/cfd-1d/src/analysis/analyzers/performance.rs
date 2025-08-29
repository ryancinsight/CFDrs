//! Performance analysis for network components

use super::traits::NetworkAnalyzer;
use crate::analysis::PerformanceMetrics;
use crate::network::Network;
use cfd_core::error::Result;
use nalgebra::RealField;
use num_traits::{Float, FromPrimitive};
use petgraph::visit::EdgeRef;
use std::iter::Sum;

/// Performance analyzer for network components
pub struct PerformanceAnalyzer<T: RealField + Copy> {
    _phantom: std::marker::PhantomData<T>,
}

impl<T: RealField + Copy> PerformanceAnalyzer<T> {
    /// Create new performance analyzer
    pub fn new() -> Self {
        Self {
            _phantom: std::marker::PhantomData,
        }
    }
}

impl<T: RealField + Copy + FromPrimitive + Float + Sum> NetworkAnalyzer<T>
    for PerformanceAnalyzer<T>
{
    type Result = PerformanceMetrics<T>;

    fn analyze(&mut self, network: &Network<T>) -> Result<PerformanceMetrics<T>> {
        let mut metrics = PerformanceMetrics::new();

        // Calculate total pressure drop
        let pressure_drop = self.calculate_total_pressure_drop(network);
        metrics.set_total_pressure_drop(pressure_drop);

        // Calculate total flow rate
        let flow_rate = self.calculate_total_flow_rate(network);
        metrics.set_total_flow_rate(flow_rate);

        // Calculate network efficiency
        let efficiency = self.calculate_efficiency(network, pressure_drop, flow_rate);
        metrics.set_efficiency(efficiency);

        // Calculate power consumption
        let power = pressure_drop * flow_rate;
        metrics.set_power_consumption(power);

        Ok(metrics)
    }

    fn name(&self) -> &str {
        "PerformanceAnalyzer"
    }
}

impl<T: RealField + Copy + FromPrimitive + Float + Sum> PerformanceAnalyzer<T> {
    fn calculate_total_pressure_drop(&self, network: &Network<T>) -> T {
        let pressures = network.pressures();
        if pressures.is_empty() {
            return T::zero();
        }

        // Find max and min pressures
        let max_pressure = pressures
            .values()
            .copied()
            .fold(T::zero(), |a, b| if a > b { a } else { b });
        let min_pressure =
            pressures
                .values()
                .copied()
                .fold(max_pressure, |a, b| if a < b { a } else { b });

        max_pressure - min_pressure
    }

    fn calculate_total_flow_rate(&self, network: &Network<T>) -> T {
        // Sum absolute flow rates at outlets
        use petgraph::graph::NodeIndex;
        let mut total = T::zero();
        for (idx, node) in network.nodes().enumerate() {
            if matches!(node.node_type, crate::network::NodeType::Outlet) {
                let node_idx = NodeIndex::new(idx);
                for edge_ref in network.graph.edges(node_idx) {
                    let edge_idx = edge_ref.id();
                    if let Some(&flow) = network.flow_rates().get(&edge_idx) {
                        total = total + Float::abs(flow);
                    }
                }
            }
        }

        total
    }

    fn calculate_efficiency(&self, network: &Network<T>, pressure_drop: T, flow_rate: T) -> T {
        if pressure_drop <= T::zero() || flow_rate <= T::zero() {
            return T::zero();
        }

        // Calculate theoretical minimum power
        let theoretical_power = pressure_drop * flow_rate;

        // Calculate actual power (including losses)
        let actual_power = self.calculate_actual_power(network);

        if actual_power > T::zero() {
            theoretical_power / actual_power
        } else {
            T::one()
        }
    }

    fn calculate_actual_power(&self, network: &Network<T>) -> T {
        // Sum power losses in all components
        let mut total_power = T::zero();

        for edge in network.edges_with_properties() {
            let flow_rate = edge.flow_rate;
            if flow_rate != T::zero() {
                let (from_idx, to_idx) = edge.nodes;
                let pressures = network.pressures();

                if let (Some(&p_from), Some(&p_to)) =
                    (pressures.get(&from_idx), pressures.get(&to_idx))
                {
                    let pressure_drop = Float::abs(p_from - p_to);
                    total_power = total_power + pressure_drop * Float::abs(flow_rate);
                }
            }
        }

        total_power
    }
}
