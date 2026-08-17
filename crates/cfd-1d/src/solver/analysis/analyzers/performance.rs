//! Performance analysis for network components

use super::traits::NetworkAnalyzer;
use crate::domain::network::Network;
use crate::solver::analysis::PerformanceMetrics;
use aequitas::systems::si::quantities::{Dimensionless, Power, Pressure, VolumetricFlowRate};
use cfd_core::conversion::SafeFromUsize;
use cfd_core::error::Result;
use cfd_core::CfdScalar;
use eunomia::NumericElement;
use petgraph::visit::EdgeRef;
use std::iter::Sum;

/// Performance analyzer for network components
pub struct PerformanceAnalyzer<T: CfdScalar + Copy> {
    _phantom: std::marker::PhantomData<T>,
}

impl<T: CfdScalar + Copy> Default for PerformanceAnalyzer<T> {
    fn default() -> Self {
        Self::new()
    }
}

impl<T: CfdScalar + Copy> PerformanceAnalyzer<T> {
    /// Create new performance analyzer
    #[must_use]
    pub fn new() -> Self {
        Self {
            _phantom: std::marker::PhantomData,
        }
    }
}

impl<T: CfdScalar + Copy + NumericElement + SafeFromUsize + Sum> NetworkAnalyzer<T>
    for PerformanceAnalyzer<T>
{
    type Result = PerformanceMetrics<T>;

    fn analyze(&mut self, network: &Network<T>) -> Result<PerformanceMetrics<T>> {
        let mut metrics = PerformanceMetrics::new();

        // Calculate total pressure drop
        let pressure_drop = self.calculate_total_pressure_drop(network);
        metrics.set_total_pressure_drop(Pressure::from_base(pressure_drop));

        // Calculate total flow rate
        let flow_rate = self.calculate_total_flow_rate(network);
        metrics.set_total_flow_rate(VolumetricFlowRate::from_base(flow_rate));

        // Calculate network efficiency
        let efficiency = self.calculate_efficiency(network, pressure_drop, flow_rate);
        metrics.set_efficiency(Dimensionless::from_base(efficiency));

        // Calculate power consumption
        let power = pressure_drop * flow_rate;
        metrics.set_power_consumption(Power::from_base(power));

        Ok(metrics)
    }

    fn name(&self) -> &'static str {
        "PerformanceAnalyzer"
    }
}

impl<T: CfdScalar + Copy + NumericElement + SafeFromUsize + Sum> PerformanceAnalyzer<T> {
    fn calculate_total_pressure_drop(&self, network: &Network<T>) -> T {
        let pressures = network.pressures();
        if pressures.is_empty() {
            return T::ZERO;
        }

        // Find max and min pressures
        let max_pressure = pressures
            .iter()
            .map(|pressure| pressure.into_base())
            .fold(T::ZERO, |a, b| if a > b { a } else { b });
        let min_pressure = pressures
            .iter()
            .map(|pressure| pressure.into_base())
            .fold(max_pressure, |a, b| if a < b { a } else { b });

        max_pressure - min_pressure
    }

    fn calculate_total_flow_rate(&self, network: &Network<T>) -> T {
        // Sum absolute flow rates at outlets
        use petgraph::graph::NodeIndex;
        let mut total = T::ZERO;
        for (idx, node) in network.nodes().enumerate() {
            if matches!(node.node_type, crate::domain::network::NodeType::Outlet) {
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

    fn calculate_efficiency(&self, network: &Network<T>, pressure_drop: T, flow_rate: T) -> T {
        if pressure_drop <= T::ZERO || flow_rate <= T::ZERO {
            return T::ZERO;
        }

        // Calculate theoretical minimum power
        let theoretical_power = pressure_drop * flow_rate;

        // Calculate actual power (including losses)
        let actual_power = self.calculate_actual_power(network);

        if actual_power > T::ZERO {
            theoretical_power / actual_power
        } else {
            T::ONE
        }
    }

    fn calculate_actual_power(&self, network: &Network<T>) -> T {
        // Sum power losses in all components
        let mut total_power = T::ZERO;

        for edge in network.edges_with_properties() {
            let flow_rate = edge.flow_rate.into_base();
            if flow_rate != T::ZERO {
                let (from_idx, to_idx) = edge.nodes;
                let pressures = network.pressures();

                if let (Some(&p_from), Some(&p_to)) = (
                    pressures.get(from_idx.index()),
                    pressures.get(to_idx.index()),
                ) {
                    let pressure_drop =
                        <T as NumericElement>::abs(p_from.into_base() - p_to.into_base());
                    total_power += pressure_drop * <T as NumericElement>::abs(flow_rate);
                }
            }
        }

        total_power
    }
}
