//! Pressure analysis for network components

use super::traits::NetworkAnalyzer;
use crate::analysis::PressureAnalysis;
use crate::network::Network;
use cfd_core::error::Result;
use nalgebra::RealField;
use num_traits::FromPrimitive;
use std::iter::Sum;

/// Pressure analyzer for network components
pub struct PressureAnalyzer<T: RealField + Copy> {
    _phantom: std::marker::PhantomData<T>,
}

impl<T: RealField + Copy> Default for PressureAnalyzer<T> {
    fn default() -> Self {
        Self::new()
    }
}

impl<T: RealField + Copy> PressureAnalyzer<T> {
    /// Create new pressure analyzer
    #[must_use]
    pub fn new() -> Self {
        Self {
            _phantom: std::marker::PhantomData,
        }
    }
}

impl<T: RealField + Copy + FromPrimitive + Sum> NetworkAnalyzer<T> for PressureAnalyzer<T> {
    type Result = PressureAnalysis<T>;

    fn analyze(&mut self, network: &Network<T>) -> Result<PressureAnalysis<T>> {
        let mut analysis = PressureAnalysis::new();

        // Collect node pressures
        use petgraph::graph::NodeIndex;
        let pressures = network.pressures();
        for (idx, node) in network.nodes().enumerate() {
            let node_idx = NodeIndex::new(idx);
            if let Some(&pressure) = pressures.get(&node_idx) {
                analysis.add_pressure(node.id.clone(), pressure);
            }
        }

        // Collect pressure drops and gradients
        for edge in network.edges_with_properties() {
            let (from_idx, to_idx) = edge.nodes;
            if let (Some(&p_from), Some(&p_to)) = (pressures.get(&from_idx), pressures.get(&to_idx))
            {
                let pressure_drop = p_from - p_to;
                analysis.add_pressure_drop(edge.id.clone(), pressure_drop);

                // Calculate pressure gradient
                let length = edge.properties.length;
                if length > T::zero() {
                    let gradient = pressure_drop / length;
                    analysis.add_pressure_gradient(edge.id.clone(), gradient);
                }
            }
        }

        Ok(analysis)
    }

    fn name(&self) -> &'static str {
        "PressureAnalyzer"
    }
}
